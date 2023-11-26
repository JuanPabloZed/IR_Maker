import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.io.wavfile import write
from scipy.signal import spectrogram
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import pyqtgraph as pg

class Sweep_Window(QMainWindow):
    def __init__(self, parent=None):
        super(Sweep_Window, self).__init__(parent)
        Sweep_Window.resize(self, 1000, 700)
        self.setWindowTitle('Sweep generator')
        # Saving file
        labelsave = QLabel("Select saving location", self)
        labelsave.setAlignment(Qt.AlignCenter)
        labelsave.setGeometry(30, 30, 180, 30)
        self.file_save = QPushButton("Saving location",self)
        self.file_save.setGeometry(30, 60, 180, 30)
        self.save_data=' '
        self.file_save.clicked.connect(lambda: self.saveSweepDialog())
        # Freq de début
        labelfreqbeg = QLabel("Begining Frequency (Hz)", self)
        labelfreqbeg.setAlignment(Qt.AlignCenter)
        labelfreqbeg.setGeometry(220, 30, 180, 30)
        self.begin_freq = QLineEdit(self)
        self.begin_freq.setMaxLength(20)
        self.begin_freq.setPlaceholderText("Enter value")
        self.begin_freq.setGeometry(220, 60, 180, 30)
        # Freq de fin
        labelfreqend = QLabel("Ending frequency (Hz)", self)
        labelfreqend.setAlignment(Qt.AlignCenter)
        labelfreqend.setGeometry(410, 30, 180, 30)
        self.end_freq = QLineEdit(self)
        self.end_freq.setMaxLength(20)
        self.end_freq.setPlaceholderText("Enter value")
        self.end_freq.setGeometry(410, 60, 180, 30)
        # Duration of sweep
        sweeplabel=QLabel("Duration (s)",self)
        sweeplabel.setAlignment(Qt.AlignCenter)
        sweeplabel.setGeometry(600, 30, 180, 30)
        self.T = QLineEdit(self)
        self.T.setPlaceholderText("Enter value")
        self.T.setGeometry(600, 60, 180, 30)
        # Sampling rate
        labelsr = QLabel("Sampling frequency (Hz)", self)
        labelsr.setAlignment(Qt.AlignCenter)
        labelsr.setGeometry(790, 30, 180, 30)
        self.sr = QLineEdit(self)
        self.sr.setMaxLength(20)
        self.sr.setPlaceholderText("Enter value")
        self.sr.setGeometry(790, 60, 180, 30)
        
        # Generate sweep
        gen_sweep = QPushButton("Generate ESS",self)
        gen_sweep.setGeometry(220,100,570,60)       
        gen_sweep.clicked.connect(lambda : self.sweep(self.begin_freq.text(), self.end_freq.text(),self.sr.text(),self.T.text(),self.save_data))
        # plot
        self.labelgraph = QLabel("Your ESS",self)
        self.labelgraph.setAlignment(Qt.AlignCenter)
        self.labelgraph.setGeometry(30,220,940,20)
        self.graph=pg.PlotWidget(self)
        self.graph.setGeometry(30, 250, 940, 420)
        self.graph.setLabel('left', 'Amplitude')
        self.graph.setLabel('bottom', 'Time (s)')
        self.graph.setBackground('w')
        # switch spectrogramme/temporel
        self.spectro_button=QPushButton("Spectrogram of your ESS",self)
        self.spectro_button.setGeometry(220,160,570,60)
        self.spectro_button.setVisible(False)
        self.x=0
        self.spectro_button.clicked.connect(lambda : self.Spectrogram(self.x,int(self.sr.text())))
        # spectro
        self.spectro_sweep=pg.PlotWidget(self)
        self.spectro_sweep.setGeometry(30, 250, 940, 420)
        self.spectro_sweep.setLabel('bottom','Time (s)')
        self.spectro_sweep.setLabel('left','Frequency (Hz)')
        self.spectro_sweep.setBackground('w')
        self.spectro_sweep.setVisible(False)

    def saveSweepDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Select saving location","","*.wav", options=options)
        if fileName[-4:] != '.wav':
                fileName = fileName + '.wav'
        self.save_data=fileName
        self.file_save.setText(Path(fileName).name)
        return
    
    def sweepgenerator(self,f1, f2, T, sr, savepath):
        """
        Generates an Exponential Sine Sweep (ESS) and saves it in a .wav file.
        -----
        INPUTS :
            - f1 (float) : initial frequency for the sweep (Hz)
            - f2 (float) : final frequency for the sweep (Hz)
            - T (int) : duration of the sweep (sec)
            - sr (int) : sample rate of the wav file
            - savepath (str) : path where the ESS will be saved
    
        OUTPUTS : 
            None, the function only saves the wav file of the ESS.
        """
        R = np.log(f2/f1)   # sweep rate
        time = np.array([i/sr for i in range(T*sr)]) # time array for the graph
        self.x = np.array([np.sin(2*np.pi*f1*T/R*(np.exp(t/sr*R/T)-1)) for t in range(T*sr)])*32767
        self.x = self.x.astype(np.int16)
        # genertion of the sweep
        
        write(savepath, sr, self.x) # save the file
        self.graph_fct(self.x)
        return
    
    def sweep(self,beg_f,end_f,sr,T,save_data):
        # sweep parameters
        f1 = float(beg_f)
        f2 = float(end_f)
        sr = int(sr)
        T = int(T)
        savepath = save_data
        # sweep generation
        self.sweepgenerator(f1, f2, T, sr, savepath)
        self.spectro_button.setVisible(True)
        return

    def Spectrogram(self,x,fs):
        self.graph.setVisible(False)
        self.spectro_sweep.setVisible(True)
        npoutfile = np.asarray(x)
        f,t,Sxx = spectrogram(npoutfile,fs=fs,nfft=len(npoutfile)//50,nperseg=len(npoutfile)//400)
        Sxx = 20*np.log10(np.matrix.transpose(Sxx))
        img = pg.ImageItem()
        img.setImage(Sxx)
        tr = pg.Qt.QtGui.QTransform()
        tr.scale(t[-1] / np.size(Sxx, axis=0), f[-1] / np.size(Sxx, axis=1))  
        img.setTransform(tr)
        self.spectro_sweep.setLimits(xMin=0, xMax=t[-1], yMin=f[0], yMax=f[-1])
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
        self.spectro_sweep.addItem(img)
        self.spectro_sweep.showGrid(x=True,y=True)
        self.spectro_button.setText('Temporal signal')
        self.spectro_button.clicked.connect(lambda : self.replotSweep())
        return
    
    def replotSweep(self):
        self.spectro_sweep.setVisible(False)
        self.graph.setVisible(True)
        self.spectro_button.setText('Spectrogram')
        self.spectro_button.clicked.connect(lambda : self.replotSpectro())
        return
        
    def replotSpectro(self):
        self.graph.setVisible(False)
        self.spectro_sweep.setVisible(True)
        self.spectro_button.setText('Temporal signal')
        self.spectro_button.clicked.connect(lambda : self.replotSweep())
        return

        
    def graph_fct(self, data):
        self.spectro_sweep.setVisible(False)
        self.graph.setVisible(True)
        t = [x/int(self.sr.text()) for x in range(int(self.T.text())*int(self.sr.text()))]
        self.graph.clear()
        pen = pg.mkPen(color = 'b')
        self.graph.plot(t,data/max(abs(data)),pen=pen)
        #self.graph.setYRange(-1.1*max(abs(data)),1.1*max(abs(data)))
        return

def next_power_of_2(n):
        return 1 << (int(np.log2(n - 1)) + 1)