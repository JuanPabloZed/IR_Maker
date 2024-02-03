from numpy.fft import rfft,rfftfreq
from numpy import (array, log, exp, sin, pi, int16, flip, column_stack,
                   int32, float32, argmax, max as maax, asarray,pad,
                   hamming, log10, min as miin, transpose, 
                   size as siize, shape, floor)

from scipy.io.wavfile import read,write as wriite
from scipy.signal import convolve, spectrogram
from soundfile import write

from PyQt6 import uic
from PyQt6.QtGui import QPixmap,QFileSystemModel
from PyQt6.QtCore import QDir
# from PyQt6.QtCore import Qt
from PyQt6.QtCore import QUrl
from PyQt6.QtMultimedia import QAudioOutput, QMediaPlayer
from PyQt6.QtWidgets import (QMainWindow, QDialog ,QPushButton, QApplication, 
                             QMessageBox, QLineEdit, QLabel, QComboBox, 
                             QCheckBox, QFileDialog, 
                             QRadioButton, QGroupBox, QListView)

from pyqtgraph.Qt.QtGui import QTransform
from pyqtgraph import PlotWidget, mkPen, mkBrush, ImageItem, HistogramLUTItem
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

class Ui_MainWIndow(QMainWindow):
    def __init__(self,parent=None):
        super(Ui_MainWIndow,self).__init__(parent)
        
        # load ui
        uic.loadUi(resource_path("add6\main_window.ui"),self)

        # show elements
        self.about_button = self.findChild(QPushButton,"about_button")
        self.about_button.clicked.connect(lambda: self.aboutDial())
        
        self.sweep_box  = self.findChild(QGroupBox,"sweep_box")

        self.browsesweep_button = self.findChild(QPushButton,"browsesweep_button")
        self.browsesweep_button.clicked.connect(lambda: self.openFileSweep())

        self.beg_freq = self.findChild(QLineEdit,"beg_freq")
        self.beg_freq.setStyleSheet("color: #8bc34a")
        self.beg_freq.textChanged.connect(lambda: self.check_all())

        self.end_freq = self.findChild(QLineEdit,"end_freq")
        self.end_freq.setStyleSheet("color: #8bc34a")
        self.end_freq.textChanged.connect(lambda: self.check_all())

        self.begfreq_label = self.findChild(QLabel,"begfreq_label")

        self.endfreq_label = self.findChild(QLabel,"endfreq_label")
        self.end_freq.setStyleSheet("color: #8bc34a")

        self.sweepgen_button = self.findChild(QPushButton,"sweepgen_button")
        self.sweepgen_button.clicked.connect(lambda: self.sweep())

        self.response_box = self.findChild(QGroupBox,"response_box")

        self.loadfolder_button = self.findChild(QPushButton,"loadfolder_button")
        self.loadfolder_button.clicked.connect(lambda: self.openResponseFolder())

        self.files_list = self.findChild(QListView,"files_list")
        self.fileModel = QFileSystemModel()
        self.files_list.setModel(self.fileModel)
        self.files_list.selectionModel().selectionChanged.connect(lambda: self.selectInList())
        self.noselec = 1

        self.output_box = self.findChild(QGroupBox,"output_box")

        self.mpt_checkbox = self.findChild(QCheckBox,"mpt_checkbox")
        self.mpt_checkbox.stateChanged.connect(lambda: self.err_mess())
        self.mpt_checkbox.setText("MP Transform (coming soon)")
        self.mpt_checkbox.setEnabled(False)

        self.bitdepth_combo = self.findChild(QComboBox,"bitdepth_combo")
        self.bitdepth_combo.setStyleSheet("color: #8bc34a")

        self.bitdepth_label = self.findChild(QLabel,"bitdepth_label")

        self.autosave_radio = self.findChild(QRadioButton,"autosave_radio")
        self.autosave_radio.setChecked(True)
        self.autosave_radio.clicked.connect(lambda: self.browseout_button.setEnabled(False))
        self.autosave_radio.toggled.connect(lambda: self.check_all())

        self.customsave_radio = self.findChild(QRadioButton,"customsave_radio")
        self.customsave_radio.clicked.connect(lambda: self.browseout_button.setEnabled(True))

        self.browseout_button = self.findChild(QPushButton,"browseout_button")
        self.browseout_button.clicked.connect(lambda: self.customSaveOut())
        self.browseout_button.setEnabled(False)

        self.srate_label = self.findChild(QLabel,"srate_label")

        self.srate = self.findChild(QLineEdit,"srate")
        self.srate.setStyleSheet("color: #8bc34a")
        self.srate.textChanged.connect(lambda: self.check_all())

        self.createir_button = self.findChild(QPushButton,"createir_button")
        self.createir_button.setEnabled(False)
        self.createir_button.clicked.connect(lambda: self.programme())

        self.playir_button = self.findChild(QPushButton,"playir_button")
        self.playir_button.clicked.connect(lambda: self.do_nothing())
        self.playir_button.setEnabled(False)

        self.spectrogram_plot = self.findChild(PlotWidget,"spectrogram_plot")
        self.spectrogram_plot.setVisible(False)
        self.spectrogram_plot.setBackground('k')
        self.spectrogram_plot.setTitle('Spectrogram')
        self.spectrogram_plot.setLabel('bottom','Time',units='s')
        self.spectrogram_plot.setLabel('left','Frequency',units='Hz')

        self.ir_plot = self.findChild(PlotWidget,"ir_plot")
        self.ir_plot.setBackground('k')
        self.ir_plot.setTitle('Waveform')
        self.ir_plot.showGrid(x=True,y=True)
        self.ir_plot.setLabel('bottom', 'Time',units='s')
        self.ir_plot.setLabel('left', 'Amplitude')
        self.ir_plot.setMouseEnabled(y=False)

        self.spectral_plot = self.findChild(PlotWidget,"spectral_plot")
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
        self.player = QMediaPlayer()
        self.audio_output = QAudioOutput()
        self.player.setAudioOutput(self.audio_output)
        self.player.setSource(QUrl.fromLocalFile(datapath))
        self.player.play()
        return
    
    def deconvolver(self):
        # load files
        _, outfile = read(self.recordfile_path, mmap=False)
        if outfile.ndim == 2:
            if shape(outfile)[1] > 2:
                QMessageBox.critical(self,"Error","Multichannel files not supported.\nPlease select a response file in mono or stereo only.")
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
            while abs(irL[i]) < 10 and abs(irR[i]) < 10 and i<len(irL)-1:
                i += 1
            crop = True
            if i==len(irL)-1: #* if went all the way until end of the signal, leave it as is and don't crop anything
                i=0
                crop = False
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
            while abs(normirL[i]) <= 300000 and abs(i) < len(normirL)-1 and crop:
                i -= 1
            normirL = normirL[0:i]
            normirR = normirR[0:i]
            normir = column_stack((normirL, normirR))
            
            return (normir) 
            
        #if recording in mono
        elif outfile.ndim == 1:
            # removing the last zeros of the recording
            i = -1
            while abs(outfile[i]) == 0 and abs(i) < len(outfile)-1:
                i -= 1  
            outfile = outfile[0:len(outfile)+i]
    
            # normalization of the recording
            normoutfile = normalize(outfile)
    
            # "deconvolving"
            ir = convolve(normoutfile, norminvess)
            # temporal shift, same as stereo mode
            i=0
            while abs(ir[i]) < 10 and i < len(ir)-1:
                i += 1
            crop = True 
            if i == len(ir)-1: #* if went all the way until end of the signal, leave it as is and don't crop anything
                    i = 0
                    crop = False

            ir = ir[i:len(ir)]
            # normalization of the IR 
            normir = normalize(ir)
            normir = normir.astype(float32)
            i = -1
            while abs(normir[i]) <= 0.00009 and abs(i) < len(normir)-1 and crop:
                i -= 1
            normir = normir[0:i]
            
            index = argmax(normir)
            if index > 10 and crop:
                normir = normir[index-10:]

            return (normir)

    def temporalgraph_fct(self,data):
        # if mono IR
        if data.ndim == 1:
            self.ir_plot.clear()
            normdata_mono = data/maax(abs(data))
            t = [x/int(self.srate.text()) for x in range(len(normdata_mono))]
            pen = mkPen(color = '#8bc34a')
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
            pen = mkPen(color = '#bf0000')
            self.ir_plot.plot(t[0:len(dataL)],normdataL + 2.2,pen=pen)
            # right
            pen = mkPen(color = '#8bc34a')
            self.ir_plot.plot(t[0:len(dataR)],normdataR,pen=pen)
            
            # y axis ticks
            aystereo = self.ir_plot.getAxis('left')
            yticksstereo = [(0,'Right'),(2.2,'Left')]
            aystereo.setTicks([yticksstereo])
            self.ir_plot.setXRange(0,t[-1])
            self.ir_plot.showGrid(x=True, y=True)
            return

    def spectralgraph_fct(self,data):
        # ir mono ir, show FFT
        if data.ndim == 1:
            self.spectrogram_plot.setVisible(False)
            self.spectral_plot.setVisible(True)
            # get signal & compute FFT
            npoutfile = asarray(data)/maax(data)
            pad_length = npow2(npow2(len(npoutfile)))
            padded_npoutfile = pad(npoutfile,(0,pad_length-len(npoutfile)),'constant',constant_values=(0,0))
            h_panned_npoutfile = padded_npoutfile*hamming(pad_length)
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

            pen = mkPen(color = '#5e8332')
            brush = mkBrush(color='#8bc34aa8')

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
            img = ImageItem()
            img.setImage(Sxx)
            tr = QTransform()
            tr.scale(t[-1] / siize(Sxx, axis=0), f[-1] / siize(Sxx, axis=1))  
            img.setTransform(tr)
            self.spectrogram_plot.setLimits(xMin=0, xMax=t[-1], yMin=f[0], yMax=f[-1])
            hist = HistogramLUTItem()
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
        self.save_name = Path(recordName).stem + ' - IR.wav'
        self.check_all()
        return

    def openFileSweep(self):
        fileName, _ = QFileDialog.getOpenFileName(self,"Select sweep file", "","*wav")
        test = read(fileName)[1] 
        if test.ndim == 2 :
            _,c = shape(test)
            while c > 2 :
                QMessageBox.critical(self,"Error","File contains too much channels.\nPlease select sweep file only in mono or stereo.")
                fileName,_ = QFileDialog.getOpenFileName(self,"Please choose a mono or stereo file only", "",'*.wav')
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
        self.folderpath = QFileDialog.getExistingDirectory(self, 'Select recordings folder')
        self.files_list.setRootIndex(self.fileModel.setRootPath(self.folderpath))
        # self.fileModel.setFilter(QDir.NoDotAndDotDot |  QDir.Files)
        self.fileModel.setFilter(QDir.Filter.Files)
        self.fileModel.setNameFilters(['*.wav'])
        self.fileModel.setNameFilterDisables(False)
        self.loadfolder_button.setText(Path(self.folderpath).name)
        return
    
    def customSaveOut(self):
        fileName, _ = QFileDialog.getSaveFileName(self,"Select saving location","","*.wav")
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
 
class Ui_SweepGenerator(QDialog):
    def __init__(self, parent=None):
        super(Ui_SweepGenerator,self).__init__(parent)

        # load ui
        uic.loadUi(resource_path("add6\sweep_window.ui"),self)
        # show elements
        self.spectro_plot = self.findChild(PlotWidget,"spectro_plot")
        self.spectro_plot.setBackground('k')
        self.spectro_plot.setTitle('Spectrogram')
        self.spectro_plot.setLabel('left','Frequency',units='Hz')
        self.spectro_plot.setLabel('bottom','Time',units='s')

        self.sweep_plot = self.findChild(PlotWidget,"sweep_plot")
        self.sweep_plot.setBackground('k')
        self.sweep_plot.setTitle('Waveform')
        self.sweep_plot.setLabel('left','Amplitude')
        self.sweep_plot.setLabel('bottom','Time',units="s")

        self.params_box = self.findChild(QGroupBox, "params_box")

        self.label = self.findChild(QLabel,"label") 

        self.label_2 = self.findChild(QLabel,'label_2')

        self.beg_freq = self.findChild(QLineEdit,"beg_freq") 
        self.beg_freq.setStyleSheet('color : #8bc34a')
        self.beg_freq.textChanged.connect(lambda: self.check_all())

        self.end_freq = self.findChild(QLineEdit,"end_freq")
        self.end_freq.setStyleSheet('color : #8bc34a')
        self.end_freq.textChanged.connect(lambda: self.check_all())

        self.duration = self.findChild(QLineEdit,"duration") 
        self.duration.setStyleSheet('color : #8bc34a')
        self.duration.textChanged.connect(lambda: self.check_all())

        self.label_3 = self.findChild(QLabel,'label_3')

        self.srate = self.findChild(QComboBox,"srate")
        self.srate.setStyleSheet('color : #8bc34a')
        self.srate.currentIndexChanged.connect(lambda: self.check_all())

        self.label_4 = self.findChild(QLabel,"label_4")

        self.save_button = self.findChild(QPushButton,"save_button")
        self.save_button.clicked.connect(lambda: self.saveFile())

        self.label_5 = self.findChild(QLabel,"label_5")

        self.gen_button = self.findChild(QPushButton,"gen_button")
        self.gen_button.clicked.connect(lambda: self.sweep())
        self.gen_button.setEnabled(False)

        self.play_button = self.findChild(QPushButton,"play_button")
        self.play_button.setVisible(False)
        self.play_button.clicked.connect(lambda: self.do_nothing())

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
        endfade = int(0.1*int(self.srate.currentText()))
        for i in range(endfade):
            sweep[i] *= i/endfade
            sweep[-1*(i+1)] *= i/endfade

        return sweep
    
    def waveform_disp(self):
        self.sweep_plot.clear()
        normdata = self.data/maax(abs(self.data))
        t = [x/int(self.srate.currentText()) for x in range(len(normdata))]
        pen = mkPen(color = '#8bc34a')
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
        img = ImageItem()
        img.setImage(Sxx)
        tr = QTransform()
        tr.scale(t[-1] / siize(Sxx, axis=0), f[-1] / siize(Sxx, axis=1))  
        img.setTransform(tr)
        self.spectro_plot.setLimits(xMin=0, xMax=t[-1], yMin=f[0], yMax=f[-1])
        hist = HistogramLUTItem()
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
        self.player = QMediaPlayer()
        self.audio_output = QAudioOutput()
        self.player.setAudioOutput(self.audio_output)
        self.player.setSource(QUrl.fromLocalFile(data))
        # audio_output.setVolume(50)
        self.player.play()        # QSound.play(data)
        return

    def saveFile(self):
        fileName, _ = QFileDialog.getSaveFileName(self,"Select saving location","","*.wav")
        if fileName[-4:] != '.wav':
                fileName = fileName + '.wav'
        self.save_loc=fileName
        self.save_button.setText(Path(fileName).name)
        self.check_all()
        return

class abtDial(QDialog):
    def __init__(self,parent=None):
        super(abtDial,self).__init__(parent)
        #load ui file
        uic.loadUi(resource_path(R"add6\about.ui"),self)

        #show elements
        self.label = self.findChild(QLabel,"label")
        self.label_2 = self.findChild(QLabel,"label_2")
        self.label_2.setPixmap(QPixmap(resource_path("add6\irmaker.png")))
        self.label_3 = self.findChild(QLabel,"label_3")
        self.label_4 = self.findChild(QLabel,"label_4")
        self.label_5 = self.findChild(QLabel,"label_5")
        self.label_6 = self.findChild(QLabel,"label_6")

def main():
    app = QApplication(sys.argv)
    main_window = Ui_MainWIndow()
    apply_stylesheet(app,theme='dark_lightgreen.xml')
    main_window.show()
    app.exec()

if __name__ == '__main__':
    main()