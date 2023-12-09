from numpy import log,array,sin,pi,exp,int16, max as maax,log10,asarray,size as siize,min as miin,transpose
from scipy.io.wavfile import write
from scipy.signal import spectrogram

from PyQt5 import QtCore, QtWidgets
from PyQt5.QtMultimedia import QSound

import pyqtgraph as pg

from pathlib import Path

class Ui_SweepGenerator(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(Ui_SweepGenerator,self).__init__(parent)
        self.setWindowTitle('Sweep generator')
        Ui_SweepGenerator.resize(self,620,630)
       
        self.spectro_plot = pg.PlotWidget(self)
        self.spectro_plot.setEnabled(True)
        self.spectro_plot.setGeometry(QtCore.QRect(10, 400, 600, 220))
        self.spectro_plot.setObjectName("spectro_plot")
        self.spectro_plot.setBackground('w')
        self.spectro_plot.setTitle('Spectrogram')
        self.spectro_plot.setLabel('left','Frequency',units='Hz')
        self.spectro_plot.setLabel('bottom','Time',units='s')
      
        self.sweep_plot = pg.PlotWidget(self)
        self.sweep_plot.setGeometry(QtCore.QRect(10, 240, 600, 150))
        self.sweep_plot.setObjectName("sweep_plot")
        self.sweep_plot.setBackground('w')
        self.sweep_plot.setTitle('Waveform')
        self.sweep_plot.setLabel('left','Amplitude')
        self.sweep_plot.setLabel('bottom','Time',units="s")
      
        self.params_box = QtWidgets.QGroupBox(self)
        self.params_box.setGeometry(QtCore.QRect(50, 10, 201, 211))
        self.params_box.setObjectName("params_box")
        self.params_box.setTitle("Sweep parameters")
      
        self.label = QtWidgets.QLabel(self.params_box)
        self.label.setGeometry(QtCore.QRect(10, 30, 91, 21))
        self.label.setObjectName("label")
        self.label.setText("Beginning freq.")
      
        self.label_2 = QtWidgets.QLabel(self.params_box)
        self.label_2.setGeometry(QtCore.QRect(10, 60, 91, 21))
        self.label_2.setObjectName("label_2")
        self.label_2.setText("Ending freq.")
      
        self.beg_freq = QtWidgets.QLineEdit(self.params_box)
        self.beg_freq.setGeometry(QtCore.QRect(100, 30, 91, 22))
        self.beg_freq.setObjectName("beg_freq")
        self.beg_freq.setPlaceholderText("in Hz")
      
        self.end_freq = QtWidgets.QLineEdit(self.params_box)
        self.end_freq.setGeometry(QtCore.QRect(100, 60, 91, 22))
        self.end_freq.setObjectName("end_freq")
        self.end_freq.setPlaceholderText("in Hz")
      
        self.duration = QtWidgets.QLineEdit(self.params_box)
        self.duration.setGeometry(QtCore.QRect(100, 90, 91, 22))
        self.duration.setObjectName("duration")
        self.duration.setPlaceholderText("in sec")
      
        self.label_3 = QtWidgets.QLabel(self.params_box)
        self.label_3.setGeometry(QtCore.QRect(10, 90, 91, 21))
        self.label_3.setObjectName("label_3")
        self.label_3.setText("Duration")
        
        self.srate = QtWidgets.QLineEdit(self.params_box)
        self.srate.setGeometry(QtCore.QRect(100, 120, 91, 22))
        self.srate.setObjectName("srate")
        self.srate.setPlaceholderText("in Hz")
        
        self.label_4 = QtWidgets.QLabel(self.params_box)
        self.label_4.setGeometry(QtCore.QRect(10, 120, 91, 21))
        self.label_4.setObjectName("label_4")
        self.label_4.setText("Sampling freq.")
        
        self.save_button = QtWidgets.QPushButton(self.params_box,clicked=lambda: self.saveFile())
        self.save_button.setGeometry(QtCore.QRect(100, 170, 93, 28))
        self.save_button.setObjectName("save_button")
        self.save_button.setText("Browse...")
        
        self.label_5 = QtWidgets.QLabel(self.params_box)
        self.label_5.setGeometry(QtCore.QRect(10, 170, 91, 28))
        self.label_5.setObjectName("label_5")
        self.label_5.setText("Saving location")
        
        self.line = QtWidgets.QFrame(self.params_box)
        self.line.setGeometry(QtCore.QRect(30, 140, 141, 41))
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")

        self.gen_button = QtWidgets.QPushButton(self,clicked=lambda: self.sweep())
        self.gen_button.setGeometry(QtCore.QRect(300, 50, 270, 51))
        self.gen_button.setObjectName("gen_button")
        self.gen_button.setText("Generate sweep")
        
        self.play_button = QtWidgets.QPushButton(self,clicked=lambda: self.do_nothing())
        self.play_button.setGeometry(QtCore.QRect(300, 120, 270, 51))
        self.play_button.setObjectName("play_button")
        self.play_button.setVisible(False)
        self.play_button.setText("Play sweep")        
        
        self.save_loc = ''        
    
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
        R = log(f2/f1)   # sweep rate
        # time = np.array([i/sr for i in range(T*sr)]) # time array for the graph
        sweep = array([sin(2*pi*f1*T/R*(exp(t/sr*R/T)-1)) for t in range(T*sr)])*32767
        sweep = sweep.astype(int16)
        return sweep
    
    def waveform_disp(self):
        self.sweep_plot.clear()
        normdata = self.data/maax(abs(self.data))
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
        IR = asarray(self.data)
        f,t,Sxx = spectrogram(IR,fs=int(self.srate.text()),nfft=len(IR)//50,nperseg=len(IR)//400,scaling='spectrum')
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
        return
