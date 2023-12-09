import numpy as np

from scipy.io.wavfile import write
from scipy.signal import spectrogram

from PyQt5.QtWidgets import QPushButton,QLineEdit,QMainWindow,QFileDialog
from PyQt5.QtMultimedia import QSound
from PyQt5 import QtCore,uic

import pyqtgraph as pg

from pathlib import Path

class SweepWindow(QMainWindow):
    def __init__(self,parent=None):
        super(SweepWindow,self).__init__(parent)
        self.setWindowTitle('Sweep generator')
        # load ui
        uic.loadUi('c:\IRs\IR_Maker\design_sweep_generator.ui',self)

        # get elements
        # sweep params
        self.beg_freq = self.findChild(QLineEdit,'begfreq_textedit')
        self.beg_freq.setPlaceholderText('in Hz')
        self.end_freq = self.findChild(QLineEdit,'endfreq_textedit')
        self.end_freq.setPlaceholderText('in Hz')
        self.duration = self.findChild(QLineEdit,'dur_textedit')
        self.duration.setPlaceholderText('in sec')
        self.srate = self.findChild(QLineEdit,'srate_textedit')
        self.srate.setPlaceholderText('in Hz')
        # saving location
        self.save_loc = ''
        self.save_button = self.findChild(QPushButton,'savesweep_button')
        # big buttons
        self.gen_button = self.findChild(QPushButton,'gensweep_button')
        self.play_button = self.findChild(QPushButton,'playsweep_button')
        self.play_button.setVisible(False)
        self.play_button.clicked.connect(lambda: self.do_nothing())
        # plots
        self.sweep_plot = self.findChild(pg.PlotWidget,'temporal_plot')
        self.sweep_plot.setBackground('w')
        self.sweep_plot.setTitle('Waveform')
        self.sweep_plot.setLabel('left','Amplitude')
        self.sweep_plot.setLabel('bottom','Time',units="s")

        self.spectro_plot = self.findChild(pg.PlotWidget,'spectrogram_plot')
        self.spectro_plot.setBackground('w')
        self.spectro_plot.setTitle('Spectrogram')
        self.spectro_plot.setLabel('left','Frequency',units='Hz')
        self.spectro_plot.setLabel('bottom','Time',units='s')
        
        # engine
        self.save_button.clicked.connect(lambda: self.saveFile())
        self.gen_button.clicked.connect(lambda: self.sweep())
    
    def sweep(self):
        self.data = self.generate_sweep()
        write(self.save_loc, int(self.srate.text()), self.data)
        self.waveform_disp()
        self.spectro_disp()
        self.play_button.setVisible(True)
        self.play_button.clicked.disconnect()
        self.play_button.clicked.connect(lambda: self.play())
        return
    
    def generate_sweep(self):
        f1 = float(self.beg_freq.text())
        f2 = float(self.end_freq.text())
        T = int(self.duration.text())
        sr = int(self.srate.text())
        R = np.log(f2/f1)   # sweep rate
        # time = np.array([i/sr for i in range(T*sr)]) # time array for the graph
        sweep = np.array([np.sin(2*np.pi*f1*T/R*(np.exp(t/sr*R/T)-1)) for t in range(T*sr)])*32767
        sweep = sweep.astype(np.int16)
        return sweep
    
    def waveform_disp(self):
        self.sweep_plot.clear()
        normdata = self.data/np.max(abs(self.data))
        t = [x/int(self.srate.text()) for x in range(len(normdata))]
        pen = pg.mkPen(color = 'b')
        self.sweep_plot.plot(t[0:len(normdata)],normdata,pen=pen)
        # y ticks
        aymono = self.sweep_plot.getAxis('left')
        yticksmono = [(0,'Mono')]
        aymono.setTicks([yticksmono])
        self.sweep_plot.setXRange(0,t[-1])
        self.sweep_plot.showGrid(x=True, y=True)
        
    def spectro_disp(self):
        IR = np.asarray(self.data)
        f,t,Sxx = spectrogram(IR,fs=int(self.srate.text()),nfft=len(IR)//50,nperseg=len(IR)//400,scaling='spectrum')
        Sxx = 20*np.log10(np.matrix.transpose(Sxx))
        img = pg.ImageItem()
        img.setImage(Sxx)
        tr = pg.Qt.QtGui.QTransform()
        tr.scale(t[-1] / np.size(Sxx, axis=0), f[-1] / np.size(Sxx, axis=1))  
        img.setTransform(tr)
        self.spectro_plot.setLimits(xMin=0, xMax=t[-1], yMin=f[0], yMax=f[-1])
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
        fileName, _ = QFileDialog.getSaveFileName(self,"Select saving location","","*.wav")
        if fileName[-4:] != '.wav':
                fileName = fileName + '.wav'
        self.save_loc=fileName
        self.save_button.setText(Path(fileName).name)
        return
