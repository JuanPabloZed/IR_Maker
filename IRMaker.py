from numpy import array,log,exp,sin,pi,int16,flip,column_stack,int32,float32,argmax,max as maax,asarray,pad,blackman,log10,min as miin,transpose,size as siize,shape,floor
from numpy.fft import rfft,rfftfreq

from scipy.io.wavfile import read,write as wriite
from scipy.signal import convolve, spectrogram
from soundfile import write

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtMultimedia import QSound
import pyqtgraph as pg
from pathlib import Path
from qt_material import apply_stylesheet

from funcs import npow2, smooth, normalize

import sys
from os import mkdir,path

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = path.abspath(".")

    return path.join(base_path, relative_path)

class Ui_MainWIndow(QtWidgets.QMainWindow):
    def __init__(self,parent=None):
        super(Ui_MainWIndow,self).__init__(parent)
        self.setWindowTitle('IR Maker')
        self.setWindowIcon(QtGui.QIcon(resource_path("add\irmaker.png")))
        self.setFixedSize(1050,900)
                
        self.about_button = QtWidgets.QPushButton(self)
        self.about_button.setGeometry(QtCore.QRect(960, 5, 80, 22))
        self.about_button.setObjectName("about_label")
        self.about_button.setFlat(True)
        self.about_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        font = QtGui.QFont()
        font.setUnderline(True)
        self.about_button.setFont(font)
        self.about_button.setText("About")
        self.about_button.clicked.connect(lambda: self.aboutDial())

        self.layoutWidget = QtWidgets.QWidget(self)
        self.layoutWidget.setGeometry(QtCore.QRect(10, 10, 821, 221))
        self.layoutWidget.setObjectName("layoutWidget")

        self.boxes_layout = QtWidgets.QHBoxLayout(self.layoutWidget)
        self.boxes_layout.setContentsMargins(0, 0, 0, 0)
        self.boxes_layout.setObjectName("boxes_layout")
        self.sweep_box = QtWidgets.QGroupBox(self.layoutWidget)
        self.sweep_box.setObjectName("sweep_box")
        self.sweep_box.setTitle("Sweep")
        
        self.browsesweep_button = QtWidgets.QPushButton(self.sweep_box,clicked = lambda: self.openFileSweep())
        self.browsesweep_button.setGeometry(QtCore.QRect(9, 40, 250, 28))
        self.browsesweep_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.browsesweep_button.setObjectName("browsesweep_button")
        self.browsesweep_button.setText("Browse sweep file")
        
        self.beg_freq = QtWidgets.QLineEdit(self.sweep_box)
        self.beg_freq.setGeometry(QtCore.QRect(153, 85, 89, 22))
        self.beg_freq.setText("20")
        self.beg_freq.setObjectName("beg_freq")
        self.beg_freq.setPlaceholderText("in Hz")
        self.beg_freq.setStyleSheet("color: #8bc34a")
        self.beg_freq.textChanged.connect(lambda: self.check_all())

        self.end_freq = QtWidgets.QLineEdit(self.sweep_box)
        self.end_freq.setGeometry(QtCore.QRect(153, 110, 89, 22))
        self.end_freq.setText("20000")
        self.end_freq.setStyleSheet("color: #8bc34a")
        self.end_freq.setObjectName("end_freq")
        self.end_freq.setPlaceholderText("in Hz")
        self.end_freq.textChanged.connect(lambda: self.check_all())
        
        self.begfreq_label = QtWidgets.QLabel(self.sweep_box)
        self.begfreq_label.setGeometry(QtCore.QRect(15, 87, 131, 21))
        self.begfreq_label.setAlignment(QtCore.Qt.AlignRight)
        self.begfreq_label.setObjectName("begfreq_label")
        self.begfreq_label.setText("Beginning freq. (Hz)")
        
        self.endfreq_label = QtWidgets.QLabel(self.sweep_box)
        self.endfreq_label.setGeometry(QtCore.QRect(34, 112, 101, 21))
        self.endfreq_label.setAlignment(QtCore.Qt.AlignLeft)
        self.endfreq_label.setObjectName("endfreq_label")
        self.endfreq_label.setText("Ending freq. (Hz)")
        
        self.sweepgen_button = QtWidgets.QPushButton(self.sweep_box,clicked=lambda: self.sweep())
        self.sweepgen_button.setGeometry(QtCore.QRect(19, 150, 229, 61))
        self.sweepgen_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.sweepgen_button.setObjectName("sweepgen_button")
        self.sweepgen_button.setText("Sweep generator")
        
        self.boxes_layout.addWidget(self.sweep_box)

        self.response_box = QtWidgets.QGroupBox(self.layoutWidget)
        self.response_box.setObjectName("response_box")
        self.response_box.setTitle("Response")
        
        self.loadfolder_button = QtWidgets.QPushButton(self.response_box,clicked=lambda: self.openResponseFolder())
        self.loadfolder_button.setGeometry(QtCore.QRect(9, 40, 251, 28))
        self.loadfolder_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.loadfolder_button.setObjectName("loadfolder_button")
        self.loadfolder_button.setText("Select recording(s) folder")
        
        self.files_list = QtWidgets.QListView(self.response_box)
        self.files_list.setGeometry(QtCore.QRect(10, 75, 251, 136))
        self.files_list.setStatusTip("")
        self.files_list.setWhatsThis("")
        self.files_list.setHorizontalScrollMode(QtWidgets.QAbstractItemView.ScrollPerPixel)
        self.files_list.setViewMode(QtWidgets.QListView.ListMode)
        self.files_list.setSelectionRectVisible(True)
        self.files_list.setObjectName("files_list")
        self.fileModel = QtWidgets.QFileSystemModel()
        self.files_list.setModel(self.fileModel)
        self.files_list.selectionModel().selectionChanged.connect(lambda: self.selectInList())
        self.noselec = 1

        self.boxes_layout.addWidget(self.response_box)

        self.output_box = QtWidgets.QGroupBox(self.layoutWidget)
        self.output_box.setObjectName("output_box")
        self.output_box.setTitle("Output")

        self.mpt_checkbox = QtWidgets.QCheckBox(self.output_box)
        # self.mpt_checkbox.setGeometry(QtCore.QRect(10, 190, 151, 20))
        # self.mpt_checkbox.setObjectName("mpt_checkbox")
        # self.mpt_checkbox.setText("MP transform")
        # self.mpt_checkbox.stateChanged.connect(lambda: self.err_mess())

        self.bitdepth_combo = QtWidgets.QComboBox(self.output_box)
        self.bitdepth_combo.setGeometry(QtCore.QRect(138 , 165, 73, 22))
        self.bitdepth_combo.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.bitdepth_combo.setObjectName("bitdepth_combo")
        self.bitdepth_combo.setStyleSheet("color: #8bc34a")
        self.bitdepth_combo.addItem("")
        self.bitdepth_combo.addItem("")
        self.bitdepth_combo.addItem("")
        self.bitdepth_combo.setItemText(0, "16 bit")
        self.bitdepth_combo.setItemText(1, "24 bit")
        self.bitdepth_combo.setItemText(2, "32 bit") 
        self.bitdepth_combo.setCurrentIndex(1)     

        self.bitdepth_label = QtWidgets.QLabel(self.output_box)
        self.bitdepth_label.setGeometry(QtCore.QRect(70, 166, 51, 21))
        self.bitdepth_label.setObjectName("bitdepth_label")
        self.bitdepth_label.setText("Bit depth")
        
        self.autosave_radio = QtWidgets.QRadioButton(self.output_box, clicked=lambda: self.browseout_button.setEnabled(False))
        self.autosave_radio.setGeometry(QtCore.QRect(10, 25, 201, 51))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.autosave_radio.sizePolicy().hasHeightForWidth())
        self.autosave_radio.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setItalic(False)
        self.autosave_radio.setFont(font)
        self.autosave_radio.setText("Automatic output (creates \"IR\"\n"
                                    "folder in recordings folder)")
        self.autosave_radio.setIconSize(QtCore.QSize(20, 20))
        self.autosave_radio.setChecked(True)
        self.autosave_radio.setObjectName("autosave_radio")
        self.autosave_radio.toggled.connect(lambda: self.check_all())

        self.customsave_radio = QtWidgets.QRadioButton(self.output_box,clicked=lambda: self.browseout_button.setEnabled(True))
        self.customsave_radio.setGeometry(QtCore.QRect(10, 70, 161, 20))
        self.customsave_radio.setObjectName("customsave_radio")
        self.customsave_radio.setText("Custom output")
        
        self.browseout_button = QtWidgets.QPushButton(self.output_box,clicked=lambda: self.customSaveOut())
        self.browseout_button.setGeometry(QtCore.QRect(61, 97, 151, 28))
        self.browseout_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.browseout_button.setAutoFillBackground(False)
        self.browseout_button.setObjectName("browseout_button")
        self.browseout_button.setText("Browse output")
        self.browseout_button.setEnabled(False)
        
        self.srate_label = QtWidgets.QLabel(self.output_box)
        self.srate_label.setGeometry(QtCore.QRect(28, 135, 111, 21))
        self.srate_label.setObjectName("srate_label")
        self.srate_label.setText("Sample rate (Hz)")
        
        self.srate = QtWidgets.QLineEdit(self.output_box)
        self.srate.setGeometry(QtCore.QRect(138, 135, 106, 22))
        self.srate.setObjectName("srate")
        self.srate.setStyleSheet("color: #8bc34a")
        self.srate.setPlaceholderText("in Hz")
        self.srate_ss = self.srate.styleSheet()
        self.srate.textChanged.connect(lambda: self.check_all())

        self.boxes_layout.addWidget(self.output_box)

        self.createir_button = QtWidgets.QPushButton(self,clicked=lambda: self.programme())
        self.createir_button.setObjectName("createir_button")
        self.createir_button.setGeometry(840,30,201,81)
        self.createir_button.setText("Create IR")
        self.createir_button.setEnabled(False)
        
        self.playir_button = QtWidgets.QPushButton(self,clicked = lambda: self.do_nothing())
        self.playir_button.setGeometry(840, 128, 201, 81)
        self.playir_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.playir_button.setObjectName("playir_button")
        self.playir_button.setText("Play IR")
        self.playir_button.setEnabled(False)
        
        self.spectrogram_plot = pg.PlotWidget(self)
        self.spectrogram_plot.setEnabled(True)
        self.spectrogram_plot.setGeometry(QtCore.QRect(10, 560, 1030, 330))
        self.spectrogram_plot.setObjectName("spectrogram_plot")
        self.spectrogram_plot.setBackground('k')
        self.spectrogram_plot.setVisible(False)
        self.spectrogram_plot.setTitle('Spectrogram')
        self.spectrogram_plot.setLabel('bottom','Time',units='s')
        self.spectrogram_plot.setLabel('left','Frequency',units='Hz')
      
        self.ir_plot = pg.PlotWidget(self)
        self.ir_plot.setGeometry(QtCore.QRect(10, 240, 1030, 315))
        self.ir_plot.setObjectName("ir_plot")
        self.ir_plot.setBackground('k')
        self.ir_plot.setTitle('Waveform')
        self.ir_plot.showGrid(x=True,y=True)
        self.ir_plot.setLabel('bottom', 'Time',units='s')
        self.ir_plot.setLabel('left', 'Amplitude')
     
        self.spectral_plot = pg.PlotWidget(self)
        self.spectral_plot.setGeometry(QtCore.QRect(10, 560, 1030, 330))
        self.spectral_plot.setObjectName("spectral_plot")
        self.spectral_plot.setBackground('k')
        self.spectral_plot.setTitle('Spectrum')
        self.spectral_plot.setLogMode(x=True)
        self.spectral_plot.showGrid(x=True,y=True)
        self.spectral_plot.setLabel('bottom','Frequency',units='Hz')
        self.spectral_plot.setLabel('left','Amplitude',units='dB')
    
        self.recordfile_path = ''
        self.sweep_path = ''
        self.save_name = '' # for autosave
        self.save_path = '' # for custom save

    def check_all(self):
        if self.customsave_radio.isChecked() == False:
        # if ONE thing not done, disable create button
            if self.browsesweep_button.text() == "Browse sweep file" or self.browsesweep_button.text() == '' \
                or self.beg_freq.text() == '' or self.beg_freq.text() == '0'\
                or self.end_freq.text() == '' or self.end_freq.text() == '0'\
                or self.noselec == 1\
                or self.srate.text() == '' or self.srate.text() == '0':

                self.createir_button.setEnabled(False)
            # if ALL things done, enable create button
            elif self.browsesweep_button.text() != "Browse sweep file" and self.browsesweep_button.text() != '' \
                and self.beg_freq.text() != '' and self.beg_freq.text() != '0'\
                and self.end_freq.text() != '' and self.end_freq.text() != '0'\
                and self.noselec != 1\
                and self.srate.text() != '' and self.srate.text() != '0':

                self.createir_button.setEnabled(True)

        elif self.customsave_radio.isChecked() == True:
        # if ONE thing not done, disable create button
            if self.browsesweep_button.text() == "Browse sweep file" or self.browsesweep_button.text() == '' \
                or self.beg_freq.text() == '' or self.beg_freq.text() == '0'\
                or self.end_freq.text() == '' or self.end_freq.text() == '0'\
                or self.noselec == 1\
                or self.browseout_button.text() == "Browse output" or self.browseout_button.text() == ''\
                or self.srate.text() == '' or self.srate.text() == '0':
            
                self.createir_button.setEnabled(False)
            # if ALL things done, enable create button
            elif self.browsesweep_button.text() != "Browse sweep file" and self.browsesweep_button.text() != '' \
                and self.beg_freq.text() != '' and self.beg_freq.text() != '0'\
                and self.end_freq.text() != '' and self.end_freq.text() != '0'\
                and self.noselec != 1\
                and self.browseout_button.text() != "Browse output" and self.browseout_button.text() != ''\
                and self.srate.text() != '' and self.srate.text() != '0':

                self.createir_button.setEnabled(True)

    # functions section
    def programme(self):
        # compute the IR
        ### ERRORS MANAGEMENT ###


        ir = self.deconvolver()
        # custom path or auto path for saving IR file
        self.outpath = ''
        if self.autosave_radio.isChecked() == True:
            irdir_path = self.folderpath + '/IR'
            # if autosave file already created, don't create it again
            if path.isdir(irdir_path) == False:
                mkdir(irdir_path)
            self.outpath = irdir_path +'/'+ self.save_name
        elif self.customsave_radio.isChecked() == True:
            self.outpath = self.save_path

        # save the IR file depending on the desired bit depth
        if self.bitdepth_combo.currentIndex() == 0: # 16 bits
            write(self.outpath,ir,int(self.srate.text()),'PCM_16')
        elif self.bitdepth_combo.currentIndex() == 1: # 24 bits
            write(self.outpath,ir,int(self.srate.text()),'PCM_24')
        elif self.bitdepth_combo.currentIndex() == 2: # 32 bits
            write(self.outpath,ir,int(self.srate.text()),'PCM_32')
        
        # plot data
        self.temporalgraph_fct(ir)
        self.spectralgraph_fct(ir)

        self.playir_button.clicked.disconnect()
        self.playir_button.clicked.connect(lambda: self.playIR())
        self.playir_button.setEnabled(True)

        return

    def do_nothing(self):
        return

    def aboutDial(self):
        dialog = abtDial(self)
        dialog.show()

    def playIR(self):
        datapath = self.outpath
        QSound.play(datapath)
        return
    
    def deconvolver(self):
        # load files
        _, outfile = read(self.recordfile_path, mmap=False)
        print(shape(outfile))
        if outfile.ndim == 2:
            if shape(outfile)[1] > 2:
                QtWidgets.QMessageBox.critical(self,"Error","Multichannel files not supported.\nPlease select a response file in mono or stereo only.")
                return
        # sweep read
        ress, ess1 = read(self.sweep_path)
        f1 = int(self.beg_freq.text())
        f2 = int(self.end_freq.text())
        # if sweep in stereo format
        if ess1.ndim == 2 :
            ess = ess1[:,0] # converting to mono
        elif ess1.ndim == 1 : 
            ess = ess1  # don't touch anything
        
        # Sweep parameters for theoretical inverse convolution
        R = log(f2/f1)   # sweep rate
        T = len(ess)/ress   # duration (s)
    
        # Theoretical inverse convolution of the sweep
        t = array([x/ress for x in range(len(ess))])
        k = exp(t*R / T)
        invess = flip(ess)/k
        # sweep normalization
        norminvess = normalize(invess)
        # if recording in stereo
        if outfile.ndim == 2:
            # removing the last zeros of the recording
            i = -1
            while outfile[i , 1] == 0:
                i -= 1
            outfileL = outfile[0:len(outfile)+i , 0]
            outfileR = outfile[0:len(outfile)+i , 1]
    
            # normalization of the recording
            normoutfileL = normalize(outfileL)
            normoutfileR = normalize(outfileR)
            
            # "deconvolving"
            irL = convolve(normoutfileL, norminvess)
            irR = convolve(normoutfileR, norminvess)
            # temporal shift to put the begining of the IR to temporal origin
                # left channel
            i=0
            while abs(irL[i]) < 10 and abs(irR[i]) < 10:
                i += 1
            irL = irL[i:len(irL)]
                # right channel
            irR = irR[i:len(irR)]
            # merge in one stereo wav file
            ir = column_stack((irL, irR))
            # normalization of the IR for int16 encoding
            normir = normalize(ir)*2**31-1
            normir = normir.astype(int32)
            normirL = normir[:,0]
            normirR = normir[:,1]
            i = -1
            while abs(normirL[i]) <= 300000:
                i -= 1
            normirL = normirL[0:i]
            normirR = normirR[0:i]
            normir = column_stack((normirL, normirR))
            
            return (normir) 
            
        #if recording in mono
        elif outfile.ndim == 1:
            # removing the last zeros of the recording
            i = -1
            while abs(outfile[i]) == 0:
                i -= 1
            outfile = outfile[0:len(outfile)+i]
    
            # normalization of the recording
            normoutfile = normalize(outfile)
    
            # "deconvolving"
            ir = convolve(normoutfile, norminvess)
            # temporal shift, same as stereo mode
            i=0
            while abs(ir[i]) < 10:
                i += 1
            ir = ir[i:len(ir)]
            # normalization of the IR 
            normir = normalize(ir)
            normir = normir.astype(float32)
            i = -1
            while abs(normir[i]) <= 0.00009:
                i -= 1
            normir = normir[0:i]
            index = argmax(normir)
            if index > 10:
                normir = normir[index-10:]
            
            # normir = smooth(normir,10)
                
          
            return (normir)

    def temporalgraph_fct(self,data):
        # if mono IR
        if data.ndim == 1:
            self.ir_plot.clear()
            normdata_mono = data/maax(abs(data))
            t = [x/int(self.srate.text()) for x in range(len(normdata_mono))]
            pen = pg.mkPen(color = '#8bc34a')
            self.ir_plot.plot(t[0:len(normdata_mono)],normdata_mono,pen=pen)
            # y ticks
            aymono = self.ir_plot.getAxis('left')
            yticksmono = [(0,'Mono')]
            aymono.setTicks([yticksmono])
            self.ir_plot.setXRange(0,t[-1])
            self.ir_plot.showGrid(x=True, y=True)
            return
        
        # if stereo IR
        if data.ndim == 2:
            dataL = data[:,0]
            dataR = data[:,1]
            normdataL = dataL/maax(abs(data))
            normdataR = dataR/maax(abs(data))
            t = [x/int(self.srate.text()) for x in range(len(data[:,0]))]
            # plot
            self.ir_plot.clear()
            # left
            pen = pg.mkPen(color = '#bf0000')
            self.ir_plot.plot(t[0:len(dataL)],normdataL + 2.2,pen=pen)
            # right
            pen = pg.mkPen(color = '#8bc34a')
            self.ir_plot.plot(t[0:len(dataR)],normdataR,pen=pen)
            
            # y axis ticks
            aystereo = self.ir_plot.getAxis('left')
            yticksstereo = [(0,'Right'),(2.2,'Left')]
            aystereo.setTicks([yticksstereo])
            self.ir_plot.setXRange(0,t[-1])
            self.ir_plot.showGrid(x=True, y=True)
            return

    def spectralgraph_fct(self,data):
        print(f'dim ir = {data.ndim}')
        # ir mono ir, show FFT
        if data.ndim == 1:
            print('fft')
            self.spectrogram_plot.setVisible(False)
            self.spectral_plot.setVisible(True)
            # get signal & compute FFT
            npoutfile = asarray(data)/maax(data)
            pad_length = npow2(npow2(len(npoutfile)))
            padded_npoutfile = pad(npoutfile,(0,pad_length-len(npoutfile)),'constant',constant_values=(0,0))
            h_panned_npoutfile = padded_npoutfile*blackman(pad_length)
            # fft_outfile = np.fft.rfft(h_panned_npoutfile)
            fft_outfile = rfft(h_panned_npoutfile)
            fft_toplot = smooth(20*log10(abs(fft_outfile)),45)-maax(20*log10(abs(fft_outfile)))
            f = rfftfreq(pad_length,1/int(self.srate.text()))
            i=0
            j=-1
            while f[i]<float(self.beg_freq.text()):
                i += 1
            while f[j]>float(self.end_freq.text()):
                j -= 1
            f = f[i:j]
            fft_toplot = fft_toplot[i:j]

            pen = pg.mkPen(color = '#5e8332')
            brush = pg.mkBrush(color='#8bc34a')

            self.spectral_plot.clear()
            self.spectral_plot.plot(f,fft_toplot,fillLevel=1.15*miin(fft_toplot),brush=brush,pen=pen)
            self.spectral_plot.setXRange(log10(f[0]),log10(f[-1]))
            return
        
        # if stereo IR, show spectrogram
        if data.ndim == 2:
            self.spectral_plot.setVisible(False)
            self.spectrogram_plot.setVisible(True)
            # spectro de la moyenne des deux canaux
            IR = (asarray(data[:,0]) + asarray(data[:,1]))/2
            f,t,Sxx = spectrogram(IR,fs=int(self.srate.text()),nfft=len(IR)//50,nperseg=len(IR)//400,scaling='spectrum')
            Sxx = 20*log10(transpose(Sxx))
            img = pg.ImageItem()
            img.setImage(Sxx)
            tr = pg.Qt.QtGui.QTransform()
            tr.scale(t[-1] / siize(Sxx, axis=0), f[-1] / siize(Sxx, axis=1))  
            img.setTransform(tr)
            self.spectrogram_plot.setLimits(xMin=0, xMax=t[-1], yMin=f[0], yMax=f[-1])
            hist = pg.HistogramLUTItem()
            hist.setImageItem(img)
            hist.setLevels(miin(Sxx), maax(Sxx))
            hist.gradient.restoreState(
            {'mode': 'rgb',
            'ticks': [(1.0, (253, 231, 36, 255)),
                    (0.85, (94, 201, 97, 255)),
                    (0.65, (32, 144, 140, 255)),
                    (0.47, (58, 82, 139, 255)),
                    (0.0, (68, 1, 84, 255))]})
            self.spectrogram_plot.clear()
            self.spectrogram_plot.addItem(img)
            self.spectrogram_plot.showGrid(x=True,y=True)
            self.spectrogram_plot.setYRange(int(self.beg_freq.text()),int(self.end_freq.text()))
            return
        
    def selectInList(self):
        self.noselec=0
        for index in self.files_list.selectedIndexes():
            recordName = self.files_list.model().fileName(index)
            self.recordfile_path = self.files_list.model().filePath(index)
            print(self.recordfile_path)
        self.save_name = Path(recordName).stem + ' - IR.wav'
        self.check_all()
        return

    def openFileSweep(self):
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,"Select sweep file", "","*wav")
        test = read(fileName)[1] 
        if test.ndim == 2 :
            _,c = shape(test)
            while c > 2 :
                QtWidgets.QMessageBox.critical(self,"Error","File contains too much channels.\nPlease select sweep file only in mono or stereo.")
                fileName,_ = QtWidgets.QFileDialog.getOpenFileName(self,"Please choose a mono or stereo file only", "",'*.wav')
                test = read(fileName)[1] 
                if test.ndim == 1:
                    c=1
                else :
                    _,c = shape(test)
        self.sweep_path=fileName
        
        self.browsesweep_button.setText(Path(fileName).name)
        self.check_all()
        return

    def openResponseFolder(self):
        self.folderpath = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select recordings folder')
        self.files_list.setRootIndex(self.fileModel.setRootPath(self.folderpath))
        self.fileModel.setFilter(QtCore.QDir.NoDotAndDotDot |  QtCore.QDir.Files)
        self.fileModel.setNameFilters(['*.wav'])
        self.fileModel.setNameFilterDisables(False)
        self.loadfolder_button.setText(Path(self.folderpath).name)
        return
    
    def customSaveOut(self):
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self,"Select saving location","","*.wav")
        if fileName[-4:] != '.wav':
                fileName = fileName + '.wav'
        self.save_path=fileName
        self.browseout_button.setText(Path(fileName).name)
        self.check_all()
        return
    
    def sweep(self):
        dialog=Ui_SweepGenerator(self)
        dialog.show()
        return
 
class Ui_SweepGenerator(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(Ui_SweepGenerator,self).__init__(parent)
        self.setWindowTitle('Sweep generator')
        self.setWindowIcon(QtGui.QIcon(resource_path("add\irmaker.png")))
        self.setFixedSize(620,630)
       
        self.spectro_plot = pg.PlotWidget(self)
        self.spectro_plot.setEnabled(True)
        self.spectro_plot.setGeometry(QtCore.QRect(10, 400, 600, 220))
        self.spectro_plot.setObjectName("spectro_plot")
        self.spectro_plot.setBackground('k')
        self.spectro_plot.setTitle('Spectrogram')
        self.spectro_plot.setLabel('left','Frequency',units='Hz')
        self.spectro_plot.setLabel('bottom','Time',units='s')
      
        self.sweep_plot = pg.PlotWidget(self)
        self.sweep_plot.setGeometry(QtCore.QRect(10, 240, 600, 150))
        self.sweep_plot.setObjectName("sweep_plot")
        self.sweep_plot.setBackground('k')
        self.sweep_plot.setTitle('Waveform')
        self.sweep_plot.setLabel('left','Amplitude')
        self.sweep_plot.setLabel('bottom','Time',units="s")
      
        self.params_box = QtWidgets.QGroupBox(self)
        self.params_box.setGeometry(QtCore.QRect(40, 10, 221, 211))
        self.params_box.setObjectName("params_box")
        self.params_box.setTitle("Sweep parameters")
      
        self.label = QtWidgets.QLabel(self.params_box)
        self.label.setGeometry(QtCore.QRect(10, 35, 121, 21))
        self.label.setObjectName("label")
        self.label.setText("Beginning freq. (Hz)")
      
        self.label_2 = QtWidgets.QLabel(self.params_box)
        self.label_2.setGeometry(QtCore.QRect(10, 65, 121, 21))
        self.label_2.setObjectName("label_2")
        self.label_2.setText("Ending freq. (Hz)")
      
        self.beg_freq = QtWidgets.QLineEdit(self.params_box)
        self.beg_freq.setGeometry(QtCore.QRect(130, 35, 81, 22))
        self.beg_freq.setObjectName("beg_freq")
        self.beg_freq.setText("20")
        self.beg_freq.setStyleSheet('color : #8bc34a')
        self.beg_freq.setPlaceholderText("in Hz")
        self.beg_freq.textChanged.connect(lambda: self.check_all())
      
        self.end_freq = QtWidgets.QLineEdit(self.params_box)
        self.end_freq.setGeometry(QtCore.QRect(130, 65, 81, 22))
        self.end_freq.setObjectName("end_freq")
        self.end_freq.setText("20000")
        self.end_freq.setStyleSheet('color : #8bc34a')
        self.end_freq.setPlaceholderText("in Hz")
        self.end_freq.textChanged.connect(lambda: self.check_all())
      
        self.duration = QtWidgets.QLineEdit(self.params_box)
        self.duration.setGeometry(QtCore.QRect(130, 95, 81, 22))
        self.duration.setObjectName("duration")
        self.duration.setText("12")
        self.duration.setStyleSheet('color : #8bc34a')
        self.duration.setPlaceholderText("in sec")
        self.duration.textChanged.connect(lambda: self.check_all())
      
        self.label_3 = QtWidgets.QLabel(self.params_box)
        self.label_3.setGeometry(QtCore.QRect(10, 95, 121, 21))
        self.label_3.setObjectName("label_3")
        self.label_3.setText("Duration (sec)")
        
        self.srate = QtWidgets.QComboBox(self.params_box)
        self.srate.setGeometry(QtCore.QRect(130, 125, 81, 22))
        self.srate.setObjectName("srate")
        self.srate.addItem('44100')
        self.srate.addItem('48000')
        self.srate.addItem('88200')
        self.srate.addItem('96000')
        self.srate.addItem('176400')
        self.srate.addItem('192000')
        self.srate.setStyleSheet('color : #8bc34a')
        self.srate.setPlaceholderText("in Hz")
        self.srate.currentIndexChanged.connect(lambda: self.check_all())
        
        self.label_4 = QtWidgets.QLabel(self.params_box)
        self.label_4.setGeometry(QtCore.QRect(10, 125, 121, 21))
        self.label_4.setObjectName("label_4")
        self.label_4.setText("Sample rate (Hz)")
        
        self.save_button = QtWidgets.QPushButton(self.params_box,clicked=lambda: self.saveFile())
        self.save_button.setGeometry(QtCore.QRect(20, 170, 181, 28))
        self.save_button.setObjectName("save_button")
        self.save_button.setText("Browse...")
        
        self.label_5 = QtWidgets.QLabel(self.params_box)
        self.label_5.setGeometry(QtCore.QRect(0, 145, 221, 28))
        self.label_5.setObjectName("label_5")
        self.label_5.setText("Saving location")
        self.label_5.setAlignment(QtCore.Qt.AlignCenter)

        self.gen_button = QtWidgets.QPushButton(self,clicked=lambda: self.sweep())
        self.gen_button.setGeometry(QtCore.QRect(300, 50, 270, 51))
        self.gen_button.setObjectName("gen_button")
        self.gen_button.setText("Generate sweep")
        self.gen_button.setEnabled(False)
        
        self.play_button = QtWidgets.QPushButton(self,clicked=lambda: self.do_nothing())
        self.play_button.setGeometry(QtCore.QRect(300, 120, 270, 51))
        self.play_button.setObjectName("play_button")
        self.play_button.setVisible(False)
        self.play_button.setText("Play sweep")        
        
        self.save_loc = ''        
    
    def check_all(self):
        if self.beg_freq.text() == '' or self.beg_freq.text() == '0'\
            or self.end_freq.text() == '' or self.end_freq.text() == '0'\
            or self.duration.text() == '' or self.duration.text() == '0'\
            or self.save_button.text() == 'Browse...' or self.save_button.text() == '':
            self.gen_button.setEnabled(False)
        
        elif self.beg_freq.text() != '' and self.beg_freq.text() != '0'\
            and self.end_freq.text() != '' and self.end_freq.text() != '0'\
            and self.duration.text() != '' and self.duration.text() != '0'\
            and self.save_button.text() != 'Browse...' and self.save_button.text() != '':
            self.gen_button.setEnabled(True)

    def sweep(self):
        self.data = self.generate_sweep()
        wriite(self.save_loc, int(self.srate.currentText()), self.data)
        self.waveform_disp()
        self.spectro_disp()
        self.play_button.setVisible(True)
        self.play_button.clicked.disconnect()
        self.play_button.clicked.connect(lambda: self.play())
        return
    
    def generate_sweep(self):
        f1 = float(self.beg_freq.text())
        f2 = float(self.end_freq.text())
        T = float(self.duration.text())
        sr = int(self.srate.currentText())
        R = log(f2/f1)   # sweep rate
        # time = np.array([i/sr for i in range(T*sr)]) # time array for the graph
        sweep = array([sin(2*pi*f1*T/R*(exp(t/sr*R/T)-1)) for t in range(int(floor(T*sr)))])*32767
        sweep = sweep.astype(int16)
        endfade = int(0.01*int(sr)*T)
        for i in range(endfade):
            sweep[i] *= i/endfade
            sweep[-1*(i+1)] *= i/endfade

        return sweep
    
    def waveform_disp(self):
        self.sweep_plot.clear()
        normdata = self.data/maax(abs(self.data))
        t = [x/int(self.srate.currentText()) for x in range(len(normdata))]
        pen = pg.mkPen(color = '#8bc34a')
        self.sweep_plot.plot(t[0:len(normdata)],normdata,pen=pen)
        # y ticks
        aymono = self.sweep_plot.getAxis('left')
        yticksmono = [(0,'Mono')]
        aymono.setTicks([yticksmono])
        self.sweep_plot.setXRange(0,t[-1])
        self.sweep_plot.showGrid(x=True, y=True)
        
    def spectro_disp(self):
        IR = asarray(self.data)
        f,t,Sxx = spectrogram(IR,fs=int(self.srate.currentText()),nfft=len(IR)//50,nperseg=len(IR)//400,scaling='spectrum')
        Sxx = 20*log10(transpose(Sxx))
        img = pg.ImageItem()
        img.setImage(Sxx)
        tr = pg.Qt.QtGui.QTransform()
        tr.scale(t[-1] / siize(Sxx, axis=0), f[-1] / siize(Sxx, axis=1))  
        img.setTransform(tr)
        self.spectro_plot.setLimits(xMin=0, xMax=t[-1], yMin=f[0], yMax=f[-1])
        hist = pg.HistogramLUTItem()
        hist.setImageItem(img)
        hist.setLevels(miin(Sxx), maax(Sxx))
        hist.gradient.restoreState(
        {'mode': 'rgb',
        'ticks': [(1.0, (253, 231, 36, 255)),
                (0.85, (94, 201, 97, 255)),
                (0.65, (32, 144, 140, 255)),
                (0.47, (58, 82, 139, 255)),
                (0.0, (68, 1, 84, 255))]})
        self.spectro_plot.clear()
        self.spectro_plot.addItem(img)
        self.spectro_plot.showGrid(x=True,y=True)
        self.spectro_plot.setYRange(int(self.beg_freq.text()),int(self.end_freq.text()))
        return
    
    def do_nothing(self):
        return

    def play(self):
        data = self.save_loc
        QSound.play(data)
        return

    def saveFile(self):
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self,"Select saving location","","*.wav")
        if fileName[-4:] != '.wav':
                fileName = fileName + '.wav'
        self.save_loc=fileName
        self.save_button.setText(Path(fileName).name)
        self.check_all()
        return

class abtDial(QtWidgets.QMainWindow):
    def __init__(self,parent=None):
        super(abtDial,self).__init__(parent)
        self.setWindowTitle('About')
        self.setWindowIcon(QtGui.QIcon(resource_path('add\irmaker.png')))
        self.setFixedSize(400,300)
        self.label = QtWidgets.QLabel(self)
        self.label.setGeometry(QtCore.QRect(10, 260, 380, 16))
        self.label.setAlignment(QtCore.Qt.AlignHCenter)
        self.label.setObjectName("label")
        self.label.setText('Version 1.0.0')

        self.label_2 = QtWidgets.QLabel(self)
        self.label_2.setGeometry(QtCore.QRect(140, 10, 120, 120))
        self.label_2.setText("")
        self.label_2.setPixmap(QtGui.QPixmap(resource_path("add\irmaker.png")))
        self.label_2.setScaledContents(True)
        self.label_2.setObjectName("label_2")

        self.label_3 = QtWidgets.QLabel(self)
        self.label_3.setGeometry(QtCore.QRect(170, 130, 61, 16))
        self.label_3.setObjectName("label_3")
        self.label_3.setText('IR Maker™')

        self.label_5 = QtWidgets.QLabel(self)
        self.label_5.setGeometry(QtCore.QRect(35, 222, 89, 16))
        self.label_5.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_5.setObjectName("label_5")
        self.label_5.setText('Git repository : ')

        self.label_6 = QtWidgets.QLabel(self)
        self.label_6.setGeometry(QtCore.QRect(123, 222, 248, 16))
        font = QtGui.QFont()
        font.setUnderline(True)
        self.label_6.setFont(font)
        self.label_6.setTextFormat(QtCore.Qt.MarkdownText)
        self.label_6.setText('https://github.com/JuanPabloZed/IR_Maker')
        self.label_6.setOpenExternalLinks(True)
        self.label_6.setObjectName("label_6")

        self.label_7 = QtWidgets.QLabel(self)
        self.label_7.setGeometry(QtCore.QRect(10, 170, 380, 34))
        self.label_7.setAlignment(QtCore.Qt.AlignCenter)
        self.label_7.setObjectName("label_7")
        self.label_7.setText("Copyright (c) 2023 Nathan Zwahlen, Benjamin Quiédeville\n& Hugo Perrier")


def main():
    app = QtWidgets.QApplication(sys.argv)
    main_window = Ui_MainWIndow()
    apply_stylesheet(app,theme='dark_lightgreen.xml')
    main_window.show()
    app.exec_()

if __name__ == '__main__':
    main()