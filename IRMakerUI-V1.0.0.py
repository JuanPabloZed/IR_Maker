# IR Maker®\nDevelopped by Nathan Zwahlen, Benjamin Quiédeville & Hugo Perrier.
# V1.0.0
#Deconvolver maison
import numpy as np
from pathlib import Path
from scipy.io.wavfile import read, write
from scipy.signal import convolve
from scipy import signal
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import pyqtgraph as pg
import sys
from tkinter import filedialog

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setWindowTitle("Ir Maker")
        MainWindow.resize(self, 1000, 730)
        #Creator Info
        labelcreator=QLabel("IR Maker® Developped by Nathan Zwahlen, Benjamin Quiédeville & Hugo Perrier.",self)
        labelcreator.setGeometry(490,700,500,15)
        labelcreator.setAlignment(Qt.AlignRight)
        labelversion=QLabel("V1.0.0",self)
        labelversion.setGeometry(490,715,500,15)
        labelversion.setAlignment(Qt.AlignRight)
        #Choix du fichier du sweep
        labelsweep = QLabel("Select sweep file", self)
        labelsweep.setAlignment(Qt.AlignCenter)
        labelsweep.setGeometry(30, 30, 180, 30)
        self.file_sweep=QPushButton("Select file",self)
        self.file_sweep.setGeometry(30, 60, 180, 30)
        self.sweep_data=' '
        self.file_sweep.clicked.connect(lambda : self.openFileSweep())
        
        #Choix du deuxième fichier 
        labelresp = QLabel("Select response file", self)
        labelresp.setAlignment(Qt.AlignCenter)
        labelresp.setGeometry(30, 100, 180, 30)
        self.file_response = QPushButton("Select file",self)
        self.file_response.setGeometry(30, 130, 180, 30)
        self.response_data=' '
        self.file_response.clicked.connect(lambda: self.openFileResponse())
        #Choix fichier de save
        labelsave = QLabel("Select saving location", self)
        labelsave.setAlignment(Qt.AlignCenter)
        labelsave.setGeometry(30, 170, 180, 30)
        self.file_save = QPushButton("Saving location",self)
        self.file_save.setGeometry(30, 200, 180, 30)
        self.save_data=' '
        self.file_save.clicked.connect(lambda: self.saveFileDialog())
        #Freq de début
        labelfreqbeg = QLabel("Begining frequency (Hz)", self)
        labelfreqbeg.setAlignment(Qt.AlignCenter)
        labelfreqbeg.setGeometry(220, 30, 180, 30)
        begin_freq = QLineEdit(self)
        begin_freq.setMaxLength(20)
        begin_freq.setPlaceholderText("Enter value")
        begin_freq.setGeometry(220, 60, 180, 30)
        #Freq de fin
        labelfreqend = QLabel("Ending frequency (Hz)", self)
        labelfreqend.setAlignment(Qt.AlignCenter)
        labelfreqend.setGeometry(410, 30, 180, 30)
        end_freq = QLineEdit(self)
        end_freq.setMaxLength(20)
        end_freq.setPlaceholderText("Enter value")
        end_freq.setGeometry(410, 60, 180, 30)
        #Sampling rate                             To get info in line edit : fe.text()
        labelsr = QLabel("Sampling rate (Hz)", self)
        labelsr.setAlignment(Qt.AlignCenter)
        labelsr.setGeometry(600, 30, 180, 30)
        self.sr = QLineEdit(self)
        self.sr.setMaxLength(20)
        self.sr.setPlaceholderText("Enter value")
        self.sr.setGeometry(600, 60, 180, 30)
        #Génération du sweep 
        gen_sweep = QPushButton("Sweep generator",self)
        gen_sweep.setGeometry(800, 60, 180, 70)
        gen_sweep.clicked.connect(lambda : self.sweep())
        #Launch IR
        ir = QPushButton("Create IR",self)
        ir.setGeometry(220,100,560,30)
        ir.clicked.connect(lambda : self.programme(self.sweep_data,self.response_data,self.save_data,begin_freq.text(),end_freq.text(),self.sr.text()))
        # spectre
        self.labelgraph = QLabel("Your IR",self)
        self.labelgraph.setAlignment(Qt.AlignCenter)
        self.labelgraph.setGeometry(30,250,940,30)
        self.ir_graph=pg.PlotWidget(self)
        self.ir_graph.setGeometry(30, 280, 940, 400)
        self.ir_graph.setLabel('bottom', 'Time (s)')
        self.ir_graph.setBackground('w')

        #plot left
        # self.labelgraphL = QLabel("Left Channel",self)
        # self.labelgraphL.setAlignment(Qt.AlignCenter)
        # self.labelgraphL.setGeometry(30,250,900,20)
        # self.irL_graph=pg.PlotWidget(self)
        # self.irL_graph.setGeometry(30, 270, 900, 180)
        # self.irL_graph.setLabel('left', 'Amplitude')
        # self.irL_graph.setLabel('bottom', 'Time (s)')
        # self.irL_graph.setBackground('w')
        # #plot right
        # self.labelgraphR = QLabel("Right Channel",self)
        # self.labelgraphR.setAlignment(Qt.AlignCenter)
        # self.labelgraphR.setGeometry(30,450,900,20)
        # self.irR_graph=pg.PlotWidget(self)
        # self.irR_graph.setGeometry(30, 470, 900, 180)
        # self.irR_graph.setLabel('left', 'Amplitude')
        # self.irR_graph.setLabel('bottom', 'Time (s)')
        # self.irR_graph.setBackground('w')
        # self.irL_graph.setVisible(False)
        # self.irR_graph.setVisible(False)
        # self.labelgraphL.setVisible(False)
        # self.labelgraphR.setVisible(False)
  
    def graph(self, data):
        if data.ndim == 1:
            #)MainWindow.resize(self,1000,700)
            self.labelgraph.setVisible(True)
            self.ir_graph.setVisible(True)
            i=-1
            while abs(data[i]) <= 0.00002:
                i -= 1
            data_mono = data[0:i]
            normdata_mono = data_mono/np.max(abs(data))
            t = [x/int(self.sr.text()) for x in range(len(data_mono))]
            self.ir_graph.clear()
            pen = pg.mkPen(color = 'b')
            self.ir_graph.plot(t[0:len(data_mono)],normdata_mono,pen=pen)
            # y ticks
            ay = self.ir_graph.getAxis('left')
            yticks = [(0,'Mono')]
            ay.setTicks([yticks])

            self.ir_graph.showGrid(x=True, y=True)
            return
        
        if data.ndim==2:
            self.labelgraph.setVisible(True)
            self.ir_graph.setVisible(True)
            i = -1
            while abs(data[i,0]) <= 0.00005 and abs(data[i,1]) <= 0.00005:
                i -= 1
            dataL = data[0:i,0]
            dataR = data[0:i,1]
            normdataL = dataL/np.max(abs(data))
            normdataR = dataR/np.max(abs(data))
            t = [x/int(self.sr.text()) for x in range(len(data[:,0]))]
            # plot
            self.ir_graph.clear()
            # left
            pen = pg.mkPen(color = 'r')
            self.ir_graph.plot(t[0:len(dataL)],normdataL + 2.2,pen=pen)
            # right
            pen = pg.mkPen(color = 'b')
            self.ir_graph.plot(t[0:len(dataR)],normdataR,pen=pen)
            
            # y axis ticks
            ay = self.ir_graph.getAxis('left')
            yticks = [(0,'Right'),(2.2,'Left')]
            ay.setTicks([yticks])

            self.ir_graph.showGrid(x=True, y=True)
            # self.irL_graph.setYRange(-1.1*max(abs(dataL)),1.1*max(abs(dataL)))
            # plot right
            # self.irR_graph.clear()
            # pen = pg.mkPen(color = 'b')
            # self.irR_graph.plot(t[0:len(dataR)],dataR/max(abs(dataL)),pen=pen)
            # self.irR_graph.setYRange(-1.1*max(abs(dataR)),1.1*max(abs(dataR)))
            return 
        
    
    def openFileSweep(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"Select sweep file", "","*wav", options=options)
        test = read(fileName)[1] 
        if test.ndim == 2 :
            _,c = np.shape(test)
            while c > 2 :
                warn=QMessageBox()
                warn.setText('File contains too much channels. Please choose a mono or stereo file only.')
                warn.setTitle('Error')
                fileName, _ = QFileDialog.getOpenFileName(self,"Please choose a mono or stereo file only", "",'*.wav', options=options)
                test = read(fileName)[1] 
                if test.ndim == 1:
                    c=1
                else :
                    _,c = np.shape(test)
        self.sweep_data=fileName
        
        self.file_sweep.setText(Path(fileName).name)

        return
    
    def openFileResponse(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"Select response file", "","*.wav", options=options)
        test = read(fileName)[1] 
        if test.ndim == 2 :
            _,c = np.shape(test)
            while c > 2 :
                warn=QMessageBox()
                warn.setText('File contains too much channels. Please choose a mono or stereo file only.')
                warn.setTitle('Error')
                fileName, _ = QFileDialog.getOpenFileName(self,"Please choose a mono or stereo file only", "","*.wav", options=options)
                test = read(fileName)[1] 
                if test.ndim == 1:
                    c=1
                else :
                    _,c = np.shape(test)
        self.response_data=fileName
        self.file_response.setText(Path(fileName).name)
        return
        
    def saveFileDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Select saving location","","*.wav", options=options)
        self.save_data=fileName
        self.file_save.setText(Path(fileName).name)
        return

    def normalize(self,values):
        return values / np.max(values)
    
    def deconvolver(self,sweeppath, recordpath, targetname,f1, f2, sr):
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
        rateout, outfile = read(recordpath, mmap=False)
    
        # sweep read
        ress, ess1 = read(sweeppath)
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
            while abs(irL[i]) < 10:
                i += 1
            irL = irL[i:len(irL)]
                # right channel
            irR = irR[i:len(irR)]
            # merge in one stereo wav file
            ir = np.column_stack((irL, irR))
            # normalization of the IR for int16 encoding
            normir = self.normalize(ir)*32767
            normirL = normir[:,0]
            normirR = normir[:,1]
            i = -1
            while abs(normirL[i]) <= 13:
                i -= 1
            normirL = normirL[0:i]
            normirR = normirR[0:i]
            normir = np.column_stack((normirL, normirR))
            normir = normir.astype(np.int16)
            # save wav file
            write(targetname + '.wav', sr, normir)
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
            # normalization of the IR for int16 encoding
            normir = self.normalize(ir)*32767
            normir = normir.astype(np.int16)
            i = -1
            while abs(normir[i]) <= 13:
                i -= 1
            normir = normir[0:i]
            # save wav file
            write(targetname +'.wav', rateout, normir)
            return (normir)
            
    def programme(self,sweep_data,response_data,save_data,begin_freq,end_freq,sr):
        
        # on choppe les paramètres nécessaires au calcul
        sweeppath = sweep_data
        recordpath = response_data
        targetname = save_data
        f1 = float(begin_freq)
        f2 = float(end_freq)
        sr=int(sr)
        ir = self.deconvolver(sweeppath, recordpath, targetname, f1, f2, sr)
        # output depending of the format in mode
            # if stereo
        if ir.ndim == 2:
            # stereo IR
            self.graph(ir)
            # if mono
        elif ir.ndim == 1:
            # mono IR
            self.graph(ir)
    
    def sweep(self):
        dialog=Sweep_Window(self)
        dialog.show()

#Creating Sweep generator Window
class Sweep_Window(QMainWindow):
    def __init__(self, parent=None):
        super(Sweep_Window, self).__init__(parent)
        Sweep_Window.resize(self, 1500, 800)
        self.setWindowTitle('Sweep generator')
        #Saving file
        labelsave = QLabel("Select saving location", self)
        labelsave.setAlignment(Qt.AlignCenter)
        labelsave.setGeometry(280, 30, 180, 30)
        self.file_save = QPushButton("Saving location",self)
        self.file_save.setGeometry(280, 60, 180, 30)
        self.save_data=' '
        self.file_save.clicked.connect(lambda: self.saveSweepDialog())
        #Freq de début
        labelfreqbeg = QLabel("Begining Frequency (Hz)", self)
        labelfreqbeg.setAlignment(Qt.AlignCenter)
        labelfreqbeg.setGeometry(470, 30, 180, 30)
        begin_freq = QLineEdit(self)
        begin_freq.setMaxLength(20)
        begin_freq.setPlaceholderText("Enter value")
        begin_freq.setGeometry(470, 60, 180, 30)
        #Freq de fin
        labelfreqend = QLabel("Ending frequency (Hz)", self)
        labelfreqend.setAlignment(Qt.AlignCenter)
        labelfreqend.setGeometry(660, 30, 180, 30)
        end_freq = QLineEdit(self)
        end_freq.setMaxLength(20)
        end_freq.setPlaceholderText("Enter value")
        end_freq.setGeometry(660, 60, 180, 30)
        #Duration of sweep
        sweeplabel=QLabel("Duration (s)",self)
        sweeplabel.setAlignment(Qt.AlignCenter)
        sweeplabel.setGeometry(850, 30, 180, 30)
        self.T = QLineEdit(self)
        self.T.setPlaceholderText("Enter value")
        self.T.setGeometry(850, 60, 180, 30)
        #Sampling rate                             To get info in line edit : fe.text()
        labelsr = QLabel("Sampling frequency (Hz)", self)
        labelsr.setAlignment(Qt.AlignCenter)
        labelsr.setGeometry(1040, 30, 180, 30)
        self.sr = QLineEdit(self)
        self.sr.setMaxLength(20)
        self.sr.setPlaceholderText("Enter value")
        self.sr.setGeometry(1040, 60, 180, 30)
        
        #Generatin sweep
        gen_sweep = QPushButton("Generate ESS",self)
        gen_sweep.setGeometry(280,100,940,60)       
        gen_sweep.clicked.connect(lambda : self.sweep(begin_freq.text(), end_freq.text(),self.sr.text(),self.T.text(),self.save_data))
        #spectre
        self.labelgraph = QLabel("Your ESS",self)
        self.labelgraph.setAlignment(Qt.AlignCenter)
        self.labelgraph.setGeometry(30,220,1440,20)
        self.graph=pg.PlotWidget(self)
        self.graph.setGeometry(30, 250, 1440, 520)
        self.graph.setLabel('left', 'Amplitude')
        self.graph.setLabel('bottom', 'Time (s)')
        self.graph.setBackground('w')
        
        #spectrogramme
        self.Spectro_button=QPushButton("Spectrogram of your ESS",self)
        self.Spectro_button.setGeometry(280,160,940,60)
        self.Spectro_button.setVisible(False)
        self.x=0
        #self.Spectro_button.clicked.connect(lambda : self.Spectrogram(self.x,int(self.sr.text())))
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
        f, t, Sxx = signal.spectrogram(x, fs)
        plt.pcolormesh(t, f, Sxx)
        plt.ylabel('Frequency [Hz]')
        plt.xlabel('Time [sec]')
        plt.title("Spectrogram of your ESS")
        plt.colorbar()
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

# main
# sweeppath = "C:\IRs\Test\Sweeps\Sweep20to20k-44,1k-10sec.wav"
# recordpath = "C:\IRs\Test\Réponses réelles\Test Chainsaw 12 EQed.wav"
# targetname = "C:\IRs\Cabs\Homemade\Test Chainsaw 12 EQed.wav"
def main():
    app = QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()