import numpy as np

from scipy.io.wavfile import read
from scipy.signal import convolve, spectrogram
from soundfile import write

from PyQt5.QtWidgets import *
from PyQt5 import QtCore, uic
import pyqtgraph as pg
from pathlib import Path

from Sweep_Window import Sweep_Window
from funcs import next_power_of_2,smooth,normalize

import sys
from os import mkdir,path

# Main Window
class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setWindowTitle('IR Maker')
        # load ui made with qt designer
        uic.loadUi('design_main_window.ui',self)

        # get elements for the engine

        # sweep section
        self.browsesweep_button = self.findChild(QPushButton, 'browseSweep_button')
        self.beg_freq = self.findChild(QLineEdit, 'begfreq_txtedit')
        self.beg_freq.setPlaceholderText('in Hz')
        self.end_freq = self.findChild(QLineEdit, 'endfreq_txtedit')
        self.end_freq.setPlaceholderText('in Hz')
        self.sweepgen_button = self.findChild(QPushButton, 'sweepgen_button')
        self.sweep_path = ''
        # response section
        self.recordfile_path = ''
        self.loadfolder_button = self.findChild(QPushButton,'selectfolder_button')
        self.files_list = self.findChild(QListView,'files_list')
        self.fileModel = QFileSystemModel()
        self.files_list.setModel(self.fileModel)
        # output section
        self.autosave_radio = self.findChild(QRadioButton,'autosave_radio')
        self.customsave_radio = self.findChild(QRadioButton,'customsave_radio')
        self.browseout_button = self.findChild(QPushButton,'browseout_button')
        self.bitdepth_combo = self.findChild(QComboBox,'bitdepth_combo')
        self.srate = self.findChild(QLineEdit,'sr_txtedit')
        self.srate.setPlaceholderText('in Hz')
        self.mpt_checkbox = self.findChild(QCheckBox,'mpt_checkbox')
        # execute section
        self.createir_button = self.findChild(QPushButton,'compute_button')
        self.createir_button.setEnabled(False)
        self.playir_button = self.findChild(QPushButton,'play_ir')
        self.playir_button.setVisible(False)
        # plots
        self.ir_plot = self.findChild(pg.PlotWidget, "temporal_plot")
        self.ir_plot.setBackground('w')
        self.ir_plot.setTitle('Waveform')
        self.ir_plot.setLabel('bottom', 'Time',units='s')
        self.ir_plot.setLabel('left', 'Amplitude')

        self.spectral_plot = self.findChild(pg.PlotWidget,"spectral_plot")
        self.spectral_plot.setBackground('w')
        self.spectral_plot.setTitle('Spectrum')
        self.spectral_plot.setLogMode(x=True)
        self.spectral_plot.showGrid(x=True,y=True)
        self.spectral_plot.setLabel('bottom','Frequency',units='Hz')
        self.spectral_plot.setLabel('left','Amplitude',units='dB')
       
        self.spectrogram_img = self.findChild(pg.PlotWidget,'spectrogram_plot')
        self.spectrogram_img.setBackground('w')
        self.spectrogram_img.setVisible(False)
        self.spectrogram_img.setLabel('bottom','Time',units='s')
        self.spectrogram_img.setLabel('left','Frequency',units='Hz')

        # engine
        # sweep section
        self.browsesweep_button.clicked.connect(lambda: self.openFileSweep())
        self.sweepgen_button.clicked.connect(lambda: self.sweep())
        # response section
        self.loadfolder_button.clicked.connect(lambda: self.openResponseFolder())
        self.files_list.selectionModel().selectionChanged.connect(lambda: self.selectInList())
        # output section
        self.save_path = '' # for custom save
        self.save_name = '' # for autosave
        self.customsave_radio.clicked.connect(lambda: self.browseout_button.setEnabled(True))
        self.autosave_radio.clicked.connect(lambda: self.browseout_button.setEnabled(False))
        self.browseout_button.setEnabled(False)
        self.browseout_button.clicked.connect(lambda: self.customSaveOut())
        
        # IR creation
        self.createir_button.clicked.connect(lambda: self.programme())


        self.show()
    
    
    # functions section
    def programme(self):
        # compute the IR
        ir = self.deconvolver()
        # custom path or auto path for saving IR file
        outpath = ''
        if self.autosave_radio.isChecked() == True:
            irdir_path = self.folderpath + '/IR'
            # if autosave file already created, don't create it again
            if path.isdir(irdir_path) == False:
                mkdir(irdir_path)
            outpath = irdir_path +'/'+ self.save_name
        elif self.customsave_radio.isChecked() == True:
            outpath = self.save_path

        # save the IR file depending on the desired bit depth
        if self.bitdepth_combo.currentIndex() == 0: # 16 bits
            write(outpath,ir,int(self.srate.text()),'PCM_16')
        elif self.bitdepth_combo.currentIndex() == 1: # 24 bits
            write(outpath,ir,int(self.srate.text()),'PCM_24')
        elif self.bitdepth_combo.currentIndex() == 2: # 32 bits
            write(outpath,ir,int(self.srate.text()),'PCM_32')
        
        # plot data
        self.temporalgraph_fct(ir)
        self.spectralgraph_fct(ir)

        self.playir_button.setVisible(True)
        return
    
    def normalize(self,values):
        return values / np.max(values)
    
    def deconvolver(self):
        """
        Calulates the IR of a system in which the sine sweep contained in sweeppath has been emitted.
        Uses theoretical inverse convolution of an Exponential Sine Sweep (ESS) to compute the IR.
        Creates a mono or stereo file depending on the format of the recording of the system response contained in recordpath.
        Saves the IR in a wav file named by targetname.
        -----
        INPUTS :
            - sweeppath (str) : path to the sine sweep wav file (MUST BE ESS)
            - recordpath (str) : path to the system's response to the sine sweep
            - targetname (str) : name of the IR wav file (don't forget the '.wav' at the end)
            - f1 (float) : initial frequency for the sweep
            - f2 (float) : final frequency for the sweep
            - sr (int) : IR sample rate (44100Hz, 48000Hz, 88200Hz, 96000Hz, 192000Hz)
        OUTPUTS : 
            - ir (or irL & irR if stereo output), the ndarray containing the IR
        """
        # load files
        _, outfile = read(self.recordfile_path, mmap=False)
    
        # sweep read
        ress, ess1 = read(self.sweep_path)
        f1 = int(self.beg_freq.text())
        f2 = int(self.end_freq.text())
        # if sweep in stereo format
        if ess1.ndim == 2 :
            ess = ess1[:,0] # converting to mono
        else : 
            ess = ess1  # don't touch anything
    
        # Sweep parameters for theoretical inverse convolution
        R = np.log(f2/f1)   # sweep rate
        T = len(ess)/ress   # duration (s)
    
        # Theoretical inverse convolution of the sweep
        t = np.array([x/ress for x in range(len(ess))])
        k = np.exp(t*R / T)
        invess = np.flip(ess)/k
        # sweep normalization
        norminvess = self.normalize(invess)
        # if recording in stereo
        if outfile.ndim == 2:
            # removing the last zeros of the recording
            i = -1
            while outfile[i , 1] == 0:
                i -= 1
            outfileL = outfile[0:len(outfile)+i , 0]
            outfileR = outfile[0:len(outfile)+i , 1]
    
            # normalization of the recording
            normoutfileL = self.normalize(outfileL)
            normoutfileR = self.normalize(outfileR)
            
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
            ir = np.column_stack((irL, irR))
            # normalization of the IR for int16 encoding
            normir = self.normalize(ir)
            normir = normir.astype(np.float32)
            normirL = normir[:,0]
            normirR = normir[:,1]
            i = -1
            while abs(normirL[i]) <= 0.00015:
                i -= 1
            normirL = normirL[0:i]
            normirR = normirR[0:i]
            normir = np.column_stack((normirL, normirR))
            
            return (normir) 
            
        #if recording in mono
        elif outfile.ndim == 1:
            # removing the last zeros of the recording
            i = -1
            while abs(outfile[i]) == 0:
                i -= 1
            outfile = outfile[0:len(outfile)+i]
    
            # normalization of the recording
            normoutfile = self.normalize(outfile)
    
            # "deconvolving"
            ir = convolve(normoutfile, norminvess)
            # temporal shift, same as stereo mode
            i=0
            while abs(ir[i]) < 10:
                i += 1
            ir = ir[i:len(ir)]
            # normalization of the IR 
            normir = self.normalize(ir)
            normir = normir.astype(np.float32)
            i = -1
            while abs(normir[i]) <= 0.00009:
                i -= 1
            normir = normir[0:i]
            index = np.argmax(normir)
            if index > 20:
                normir = normir[index-20:]
          
            return (normir)

    def temporalgraph_fct(self,data):
        # if mono IR
        if data.ndim == 1:
            self.ir_plot.clear()
            normdata_mono = data/np.max(abs(data))
            t = [x/int(self.srate.text()) for x in range(len(normdata_mono))]
            pen = pg.mkPen(color = 'b')
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
            normdataL = dataL/np.max(abs(data))
            normdataR = dataR/np.max(abs(data))
            t = [x/int(self.srate.text()) for x in range(len(data[:,0]))]
            # plot
            self.ir_plot.clear()
            # left
            pen = pg.mkPen(color = 'r')
            self.ir_plot.plot(t[0:len(dataL)],normdataL + 2.2,pen=pen)
            # right
            pen = pg.mkPen(color = 'b')
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
            self.spectrogram_img.setVisible(False)
            self.spectral_plot.setVisible(True)
            # get signal & compute FFT
            npoutfile = np.asarray(data)/np.max(data)
            pad_length = next_power_of_2(next_power_of_2(len(npoutfile)))
            padded_npoutfile = np.pad(npoutfile,(0,pad_length-len(npoutfile)),'constant',constant_values=(0,0))
            h_panned_npoutfile = padded_npoutfile*np.blackman(pad_length)
            # fft_outfile = np.fft.rfft(h_panned_npoutfile)
            fft_outfile = np.fft.rfft(h_panned_npoutfile)
            fft_toplot = smooth(20*np.log10(abs(fft_outfile)),45)-np.max(20*np.log10(abs(fft_outfile)))
            f = np.fft.rfftfreq(pad_length,1/int(self.srate.text()))
            i=0
            j=-1
            while f[i]<float(self.beg_freq.text()):
                i += 1
            while f[j]>float(self.end_freq.text()):
                j -= 1
            f = f[i:j]
            fft_toplot = fft_toplot[i:j]

            pen = pg.mkPen(color = 'r')
            brush = pg.mkBrush(color=(255,150,150))

            self.spectral_plot.clear()
            self.spectral_plot.plot(f,fft_toplot,fillLevel=1.15*np.min(fft_toplot),brush=brush,pen=pen)
            self.spectral_plot.setXRange(np.log10(f[0]),np.log10(f[-1]))
            return
        
        # if stereo IR, show spectrogram
        if data.ndim == 2:
            self.spectral_plot.setVisible(False)
            self.spectrogram_img.setVisible(True)
            # spectro de la moyenne des deux canaux
            IR = (np.asarray(data[:,0]) + np.asarray(data[:,1]))/2
            f,t,Sxx = spectrogram(IR,fs=int(self.srate.text()),nfft=len(IR)//50,nperseg=len(IR)//400,scaling='spectrum')
            Sxx = 20*np.log10(np.matrix.transpose(Sxx))
            img = pg.ImageItem()
            img.setImage(Sxx)
            tr = pg.Qt.QtGui.QTransform()
            tr.scale(t[-1] / np.size(Sxx, axis=0), f[-1] / np.size(Sxx, axis=1))  
            img.setTransform(tr)
            self.spectrogram_img.setLimits(xMin=0, xMax=t[-1], yMin=f[0], yMax=f[-1])
            hist = pg.HistogramLUTItem()
            hist.setImageItem(img)
            hist.setLevels(np.min(Sxx), np.max(Sxx))
            hist.gradient.restoreState(
            {'mode': 'rgb',
            'ticks': [(1.0, (253, 231, 36, 255)),
                    (0.85, (94, 201, 97, 255)),
                    (0.65, (32, 144, 140, 255)),
                    (0.47, (58, 82, 139, 255)),
                    (0.0, (68, 1, 84, 255))]})
            self.spectrogram_img.clear()
            self.spectrogram_img.addItem(img)
            self.spectrogram_img.showGrid(x=True,y=True)
            self.spectrogram_img.setYRange(int(self.beg_freq.text()),int(self.end_freq.text()))
            return
        
    def selectInList(self):
        self.createir_button.setEnabled(True)
        for index in self.files_list.selectedIndexes():
            recordName = self.files_list.model().fileName(index)
            self.recordfile_path = self.files_list.model().filePath(index)
            print(self.recordfile_path)
        self.save_name = Path(recordName).stem + ' - IR.wav'
        self.createir_button.setText(self.save_name)
        return

    def openFileSweep(self):
        fileName, _ = QFileDialog.getOpenFileName(self,"Select sweep file", "","*wav")
        test = read(fileName)[1] 
        if test.ndim == 2 :
            _,c = np.shape(test)
            while c > 2 :
                warn=QMessageBox()
                warn.setText('File contains too much channels. Please choose a mono or stereo file only.')
                warn.setTitle('Error')
                fileName, _ = QFileDialog.getOpenFileName(self,"Please choose a mono or stereo file only", "",'*.wav')
                test = read(fileName)[1] 
                if test.ndim == 1:
                    c=1
                else :
                    _,c = np.shape(test)
        self.sweep_path=fileName
        
        self.browsesweep_button.setText(Path(fileName).name)
        return

    def openResponseFolder(self):
        self.folderpath = QFileDialog.getExistingDirectory(self, 'Select recordings folder')
        self.files_list.setRootIndex(self.fileModel.setRootPath(self.folderpath))
        self.fileModel.setFilter(QtCore.QDir.NoDotAndDotDot |  QtCore.QDir.Files)
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
        return
    
    def sweep(self):
        dialog=Sweep_Window(self)
        dialog.show()
        return
    



# run app
app = QApplication(sys.argv)
main_window = MainWindow()
app.exec_()