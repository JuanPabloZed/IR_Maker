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
        Sweep_Window.resize(self, 1500, 800)
        self.setWindowTitle('Sweep generator')
        # Saving file
        labelsave = QLabel("Select saving location", self)
        labelsave.setAlignment(Qt.AlignCenter)
        labelsave.setGeometry(280, 30, 180, 30)
        self.file_save = QPushButton("Saving location",self)
        self.file_save.setGeometry(280, 60, 180, 30)
        self.save_data=' '
        self.file_save.clicked.connect(lambda: self.saveSweepDialog())
        # Freq de début
        labelfreqbeg = QLabel("Begining Frequency (Hz)", self)
        labelfreqbeg.setAlignment(Qt.AlignCenter)
        labelfreqbeg.setGeometry(470, 30, 180, 30)
        self.begin_freq = QLineEdit(self)
        self.begin_freq.setMaxLength(20)
        self.begin_freq.setPlaceholderText("Enter value")
        self.begin_freq.setGeometry(470, 60, 180, 30)
        # Freq de fin
        labelfreqend = QLabel("Ending frequency (Hz)", self)
        labelfreqend.setAlignment(Qt.AlignCenter)
        labelfreqend.setGeometry(660, 30, 180, 30)
        self.end_freq = QLineEdit(self)
        self.end_freq.setMaxLength(20)
        self.end_freq.setPlaceholderText("Enter value")
        self.end_freq.setGeometry(660, 60, 180, 30)
        # Duration of sweep
        sweeplabel=QLabel("Duration (s)",self)
        sweeplabel.setAlignment(Qt.AlignCenter)
        sweeplabel.setGeometry(850, 30, 180, 30)
        self.T = QLineEdit(self)
        self.T.setPlaceholderText("Enter value")
        self.T.setGeometry(850, 60, 180, 30)
        # Sampling rate
        labelsr = QLabel("Sampling frequency (Hz)", self)
        labelsr.setAlignment(Qt.AlignCenter)
        labelsr.setGeometry(1040, 30, 180, 30)
        self.sr = QLineEdit(self)
        self.sr.setMaxLength(20)
        self.sr.setPlaceholderText("Enter value")
        self.sr.setGeometry(1040, 60, 180, 30)
        
        # Generate sweep
        gen_sweep = QPushButton("Generate ESS",self)
        gen_sweep.setGeometry(280,100,940,60)       
        gen_sweep.clicked.connect(lambda : self.sweep(self.begin_freq.text(), self.end_freq.text(),self.sr.text(),self.T.text(),self.save_data))
        # plot
        self.labelgraph = QLabel("Your ESS",self)
        self.labelgraph.setAlignment(Qt.AlignCenter)
        self.labelgraph.setGeometry(30,220,1440,20)
        self.graph=pg.PlotWidget(self)
        self.graph.setGeometry(30, 250, 1440, 520)
        self.graph.setLabel('left', 'Amplitude')
        self.graph.setLabel('bottom', 'Time (s)')
        self.graph.setBackground('w')
        
        # spectrogramme
        self.Spectro_button=QPushButton("Spectrogram of your ESS",self)
        self.Spectro_button.setGeometry(280,160,940,60)
        self.Spectro_button.setVisible(False)
        self.x=0
        self.Spectro_button.clicked.connect(lambda : self.Spectrogram(self.x,int(self.sr.text())))
        '''
        self.labelspectro = QLabel("Spectrogram of Generated Sweep",self)
        self.labelspectro.setAlignment(Qt.AlignCenter)
        self.labelspectro.setGeometry(30,410,1440,20)
        self.spectro=pg.PlotWidget(self)
        self.spectro.setGeometry(30, 430, 940, 230)
        '''
        
        '''
        pg.setConfigOptions(imageAxisOrder='row-major')# Interpret image data as row-major instead of col-major
        pg.mkQApp()
        self.win = pg.GraphicsLayoutWidget(self)
        self.win.setGeometry(30, 430, 1440, 330)
        self.p1 = self.win.addPlot()# A plot area (ViewBox + axes) for displaying the image
        self.img = pg.ImageItem()# Item for displaying image data
        self.p1.addItem(self.img)
        self.hist = pg.HistogramLUTItem()# Add a histogram with which to control the gradient of the image
        self.hist.setImageItem(self.img)# Link the histogram to the image
        self.win.addItem(self.hist)# If you don't add the histogram to the window, it stays invisible, but I find it useful.
        #win.show()# Show the window
        self.p1.setLabel('bottom', "Time", units='s') # Add labels to the axis
        # If you include the units, Pyqtgraph automatically scales the axis and adjusts the SI prefix (in this case kHz)
        self.p1.setLabel('left', "Frequency", units='Hz')
        '''
    def saveSweepDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Select saving location","","*.wav", options=options)
        self.file_save.setText(Path(fileName).name)
        self.save_data=fileName
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
        '''
        plt.plot(time, x, linewidth = 0.5) # plot the sweep to show the user the result
        # graph enhancement
        plt.title('Your sine sweep')
        plt.xlabel('Time(s)')
        plt.ylabel('Amplitude')
        plt.show()
        '''
        write(savepath + '.wav' , sr, self.x) # save the file
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
        self.Spectro_button.setVisible(True)
        return

    def Spectrogram(self,x,fs):
        npoutfile = np.asarray(x)
        pad_length = next_power_of_2(next_power_of_2(len(npoutfile)))
        padded_npoutfile = np.pad(npoutfile,(0,pad_length-len(npoutfile)),'constant',constant_values=(0,0))
        f,t,Sxx = spectrogram(padded_npoutfile,fs=fs,nfft=len(padded_npoutfile)//50,nperseg=len(padded_npoutfile)//400)
        plt.pcolormesh(t,f,20*np.log10(Sxx))
        cb = plt.colorbar()
        cb.remove()
        plt.ylim(15,22000)
        plt.xlim(0,len(x)/fs)
        plt.ylabel('Frequency [Hz]')
        plt.xlabel('Time [sec]')
        plt.title("Spectrogram of your ESS")
        plt.show()
    

        '''
        self.hist.setLevels(np.min(Sxx), np.max(Sxx))# Fit the min and max levels of the histogram to the data available
        # This gradient is roughly comparable to the gradient used by Matplotlib
        # You can adjust it and then save it using hist.gradient.saveState()
        self.hist.gradient.restoreState(
                {'mode': 'rgb',
                'ticks': [(0.5, (0, 182, 188, 255)),
                        (1.0, (246, 111, 0, 255)),
                        (0.0, (75, 0, 113, 255))]})
        self.img.setImage(Sxx)# Sxx contains the amplitude for each pixel
        #self.img.scale(t[-1]/np.size(Sxx, axis=1),
                #f[-1]/np.size(Sxx, axis=0))# Scale the X and Y Axis to time and frequency (standard is pixels)
        #self.p1.setLimits(xMin=0, xMax=t[-1], yMin=0, yMax=f[-1])# Limit panning/zooming to the spectrogram
        
        return
        '''
    def graph_fct(self, data):
        t = [x/int(self.sr.text()) for x in range(int(self.T.text())*int(self.sr.text()))]
        self.graph.clear()
        pen = pg.mkPen(color = 'b')
        self.graph.plot(t,data/max(abs(data)),pen=pen)
        #self.graph.setYRange(-1.1*max(abs(data)),1.1*max(abs(data)))
        return

def next_power_of_2(n):
        return 1 << (int(np.log2(n - 1)) + 1)