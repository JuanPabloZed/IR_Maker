# IR Maker®\nDevelopped by Nathan Zwahlen, Benjamin Quiédeville & Hugo Perrier.
# V1.0.0
#Deconvolver maison
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io.wavfile import read, write
from scipy.signal import convolve, spectrogram
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import pyqtgraph as pg
from Sweep_Window import Sweep_Window
import sys



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
        self.begin_freq = QLineEdit(self)
        self.begin_freq.setPlaceholderText("Enter value")
        self.begin_freq.setMaxLength(20)
        self.begin_freq.setGeometry(220, 60, 180, 30)
        #Freq de fin
        labelfreqend = QLabel("Ending frequency (Hz)", self)
        labelfreqend.setAlignment(Qt.AlignCenter)
        labelfreqend.setGeometry(410, 30, 180, 30)
        self.end_freq = QLineEdit(self)
        self.end_freq.setMaxLength(20)
        self.end_freq.setPlaceholderText("Enter value")
        self.end_freq.setGeometry(410, 60, 180, 30)
        #Sampling rate                             
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
        # Create IR
        ir = QPushButton("Create IR",self)
        ir.setGeometry(220,100,560,60)
        ir.clicked.connect(lambda : self.programme(self.sweep_data,self.response_data,self.save_data,self.begin_freq.text(),self.end_freq.text(),self.sr.text()))
        # plot
        self.labelgraph = QLabel("Your IR",self)
        self.labelgraph.setAlignment(Qt.AlignCenter)
        self.labelgraph.setGeometry(30,250,940,30)
        # plot mono
        self.ir_graphmono=pg.PlotWidget(self)
        self.ir_graphmono.setGeometry(30, 280, 940, 400)
        self.ir_graphmono.setLabel('bottom', 'Time (s)')
        self.ir_graphmono.setLabel('left', 'Amplitude')
        self.ir_graphmono.setBackground('w')
        # plot stereo
        self.ir_graphstereo=pg.PlotWidget(self)
        self.ir_graphstereo.setGeometry(30, 280, 940, 400)
        self.ir_graphstereo.setLabel('bottom', 'Time (s)')
        self.ir_graphstereo.setLabel('left', 'Amplitude')
        self.ir_graphstereo.setBackground('w')
        self.ir_graphstereo.setVisible(False)
        # Bouton spectro/fft
        self.sp_button = QPushButton('Spectro ou FFT',self)
        self.sp_button.setGeometry(800,180,180,70)
        self.sp_button.setVisible(False)
        # spectro
        self.ir_spectro=pg.PlotWidget(self)
        self.ir_spectro.setGeometry(30, 280, 940, 400)
        self.ir_spectro.setVisible(False)
        self.ir_spectro.setLabel('bottom','Time (s)')
        self.ir_spectro.setLabel('left','Frequency (Hz)')
        self.ir_spectro.setBackground('w')
        # FFT
        self.ir_fft=pg.PlotWidget(self)
        self.ir_fft.setGeometry(30, 280, 940, 400)
        self.ir_fft.setLabel('bottom','Frequency (Hz)')
        self.ir_fft.setLabel('left','Amplitude (dB)')
        self.ir_fft.setBackground('w')
        self.ir_fft.setLogMode(x=True)
        self.ir_fft.showGrid(x=True,y=True)
        self.ir_fft.setVisible(False)

    def graphIR(self, data):
        if data.ndim == 1:
            self.ir_spectro.setVisible(False)
            self.ir_fft.setVisible(False)
            self.ir_graphstereo.setVisible(False)
            self.ir_graphmono.setVisible(True)
            i=-1
            while abs(data[i]) <= 0.00002:
                i -= 1
            data_mono = data[0:i]
            normdata_mono = data_mono/np.max(abs(data))
            t = [x/int(self.sr.text()) for x in range(len(data_mono))]
            self.ir_graphmono.clear()
            pen = pg.mkPen(color = 'b')
            self.ir_graphmono.plot(t[0:len(data_mono)],normdata_mono,pen=pen)
            # y ticks
            aymono = self.ir_graphmono.getAxis('left')
            yticksmono = [(0,'Mono')]
            aymono.setTicks([yticksmono])
            self.ir_graphmono.setXRange(0,t[-1])
            self.ir_graphmono.showGrid(x=True, y=True)
            return
        
        if data.ndim==2:
            self.ir_spectro.setVisible(False)
            self.ir_fft.setVisible(False)
            self.ir_graphmono.setVisible(False)
            self.ir_graphstereo.setVisible(True)
            i = -1
            while abs(data[i,0]) <= 0.00005 and abs(data[i,1]) <= 0.00005:
                i -= 1
            dataL = data[0:i,0]
            dataR = data[0:i,1]
            normdataL = dataL/np.max(abs(data))
            normdataR = dataR/np.max(abs(data))
            t = [x/int(self.sr.text()) for x in range(len(data[:,0]))]
            # plot
            self.ir_graphstereo.clear()
            # left
            pen = pg.mkPen(color = 'r')
            self.ir_graphstereo.plot(t[0:len(dataL)],normdataL + 2.2,pen=pen)
            # right
            pen = pg.mkPen(color = 'b')
            self.ir_graphstereo.plot(t[0:len(dataR)],normdataR,pen=pen)
            
            # y axis ticks
            aystereo = self.ir_graphstereo.getAxis('left')
            yticksstereo = [(0,'Right'),(2.2,'Left')]
            aystereo.setTicks([yticksstereo])
            self.ir_graphstereo.setXRange(0,t[-1])
            self.ir_graphstereo.showGrid(x=True, y=True)
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
        sr = int(sr)
        ir = self.deconvolver(sweeppath, recordpath, targetname, f1, f2, sr)
        # output depending of the format in mode
            # if stereo
        if ir.ndim == 2:
            # stereo IR
            self.graphIR(ir)
            self.sp_button.setText('Spectrogram')
            self.sp_button.setVisible(True)
            self.sp_button.clicked.connect(lambda : self.spectroIR(ir,sr))
            # if mono
        elif ir.ndim == 1:
            # mono IR
            self.graphIR(ir)
            self.sp_button.setText('Spectrum')
            self.sp_button.setVisible(True)
            self.sp_button.clicked.connect(lambda : self.fftIR(ir,sr))
    
    def fftIR(self,irfile,srate):
        self.ir_graphstereo.setVisible(False)
        self.ir_graphmono.setVisible(False)
        self.ir_spectro.setVisible(False)
        self.ir_fft.setVisible(True)
        # get signal & compute FFT
        npoutfile = np.asarray(irfile)
        pad_length = next_power_of_2(next_power_of_2(len(npoutfile)))
        padded_npoutfile = np.pad(npoutfile,(0,pad_length-len(npoutfile)),'constant',constant_values=(0,0))
        h_panned_npoutfile = padded_npoutfile*np.blackman(pad_length)
        # fft_outfile = np.fft.rfft(h_panned_npoutfile)
        fft_outfile = np.fft.rfft(h_panned_npoutfile)
        f = np.fft.rfftfreq(pad_length,1/srate)
        
        pen = pg.mkPen(color = 'r')
        self.ir_fft.plot(f,smooth(20*np.log10(abs(fft_outfile)),45)-np.max(20*np.log10(abs(fft_outfile))),pen=pen)
        self.ir_fft.setXRange(np.log10(int(self.begin_freq.text())),np.log10(int(self.end_freq.text())))
        self.sp_button.setText('Temporal signal')
        self.sp_button.clicked.connect(lambda : self.replotIRmono())

    def replotIRmono(self):
        self.ir_fft.setVisible(False)
        self.ir_graphstereo.setVisible(False)
        self.ir_spectro.setVisible(False)
        self.ir_graphmono.setVisible(True)

        self.sp_button.setText('Spectrum')
        self.sp_button.clicked.connect(lambda : self.replotIRfft())
        return
    
    def replotIRfft(self):
        self.ir_graphstereo.setVisible(False)
        self.ir_graphmono.setVisible(False)
        self.ir_spectro.setVisible(False)
        self.ir_fft.setVisible(True)

        self.sp_button.setText('Temporal signal')
        self.sp_button.clicked.connect(lambda : self.replotIRmono())
    
    def spectroIR(self,file,srate):
        self.ir_fft.setVisible(False)
        self.ir_graphstereo.setVisible(False)
        self.ir_graphmono.setVisible(False)
        self.ir_spectro.setVisible(True)
        # spectro de la moyenne des deux canaux
        IR = (np.asarray(file[:,0]) + np.asarray(file[:,1]))/2
        f,t,Sxx = spectrogram(IR,fs=srate,nfft=len(IR)//50,nperseg=len(IR)//400,scaling='spectrum')
        Sxx = 20*np.log10(np.matrix.transpose(Sxx))
        img = pg.ImageItem()
        img.setImage(Sxx)
        tr = pg.Qt.QtGui.QTransform()
        tr.scale(t[-1] / np.size(Sxx, axis=0), f[-1] / np.size(Sxx, axis=1))  
        img.setTransform(tr)
        self.ir_spectro.setLimits(xMin=0, xMax=t[-1], yMin=f[0], yMax=f[-1])
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
        self.ir_spectro.addItem(img)
        self.ir_spectro.showGrid(x=True,y=True)
        self.sp_button.setText('Temporal signal')
        self.sp_button.clicked.connect(lambda : self.replotIRstereo())
        return

    def replotIRstereo(self):
        self.ir_fft.setVisible(False)
        self.ir_graphmono.setVisible(False)
        self.ir_spectro.setVisible(False)
        self.ir_graphstereo.setVisible(True)
        self.sp_button.setText('Spectrogram')
        self.sp_button.clicked.connect(lambda : self.replotSpectro())
        return

    def replotSpectro(self):
        self.ir_fft.setVisible(False)
        self.ir_graphstereo.setVisible(False)
        self.ir_graphmono.setVisible(False)
        self.ir_spectro.setVisible(True)
        self.sp_button.setText('Temporal signal')
        self.sp_button.clicked.connect(lambda : self.replotIRstereo())
        return

    def sweep(self):
        dialog=Sweep_Window(self)
        dialog.show()

def next_power_of_2(n):
    return 1 << (int(np.log2(n - 1)) + 1)

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
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