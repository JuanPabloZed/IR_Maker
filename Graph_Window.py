from scipy import signal
from scipy.io.wavfile import read
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
import numpy as np
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import pyqtgraph
import sys
from funcs import *

# def next_power_of_2(n):
#     return 1 << (int(np.log2(n - 1)) + 1)
mode = 'spectro'

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)

class GraphWindow(QMainWindow):
    def __init__(self, parent=None):
        super(GraphWindow,self).__init__(parent)
        self.setWindowTitle('Test')
        # GraphWindow.resize(self,1000,460)
        sc = MplCanvas(self, width=15, height=6, dpi=100)

        # self.graph = pyqtgraph.PlotWidget(self)
        # self.graph.setGeometry(30,30,940,400)
        self.plotgraph(sc,mode)
    
    def plotgraph(self,canvas,mode):
        # data pour plot spectro
        
            # x ticks for fft
        majfticks_fft=[10**power for power in range(0,5)]
        minfticks_fft=[]
        for power in range(0,5):
            minfticks_fft += [10**power * x for x in range(2,10)]
        if mode == 'fft':
        # fft
            # data pour plot normal
            sr_fft,irfile = read(R'c:\IRs\Test\cab.wav.wav')
            npoutfile = np.asarray(irfile)/np.max(irfile)
            pad_length = next_power_of_2(next_power_of_2(len(npoutfile)))
            padded_npoutfile = np.pad(npoutfile,(0,pad_length-len(npoutfile)),'constant',constant_values=(0,0))
            h_panned_npoutfile = padded_npoutfile*np.blackman(pad_length)
            fft_outfile = 20*np.log10(abs(np.fft.rfft(h_panned_npoutfile)))
            fft_outfile = smooth(fft_outfile,80)# - np.max(fft_outfile)
            f_fft = np.fft.rfftfreq(pad_length,1/sr_fft)
            # indexs pour valeurs de f_fft dans les bornes 20-22000
            i=0
            j=-1
            while f_fft[i] < 20:
                i+=1
            while f_fft[j] > 22000:
                j-=1
            limmin = 1.02*np.min(fft_outfile[i:j])
            limmax = 0.9*np.max(fft_outfile[i:j])
        # display data            
            # plot
            canvas.axes.clear()
            canvas.axes.semilogx(f_fft,fft_outfile,lw=0.7)
            canvas.axes.fill_between(f_fft,-1000,y2=fft_outfile,alpha=0.25)
            canvas.axes.grid(which='both',alpha=0.25)
            canvas.axes.set_xticks(majfticks_fft,minor=False)
            canvas.axes.set_xticks(minfticks_fft,minor=True)
            canvas.axes.set_xticklabels(majfticks_fft,minor=False,rotation=50)
            canvas.axes.set_xticklabels(minfticks_fft,minor=True,rotation=50)
            canvas.axes.set_title('Spectrum')
            canvas.axes.set_xlabel('Frequency (Hz)')
            canvas.axes.set_ylabel('Amplitude (dB)')
            canvas.axes.set_xlim(20,22000)
            canvas.axes.set_ylim(limmin,limmax)
            canvas.fig.tight_layout()
            
            
        
        # spectro
        elif mode == 'spectro':
            # data pour spectro
            canvas.axes.clear()
            sr_spectro,data = read(R'c:\IRs\_Verbs\Kohle Reverberant Room\KM183 Tight.wav')
            onedata = (np.asarray(data[:,0]) + np.asarray(data[:,1]))/2
            canvas.axes.specgram(onedata,Fs=sr_spectro,
                                 NFFT=len(onedata)//100,
                                 noverlap=len(onedata)//170,
                                 scale='dB') 
            canvas.axes.grid(which='both',alpha=1)
            canvas.axes.set_yscale('symlog')
            canvas.axes.set_yticks(majfticks_fft,minor=False)
            canvas.axes.set_yticks(minfticks_fft,minor=True)
            canvas.axes.set_yticklabels(majfticks_fft,minor=False,fontsize=7)
            canvas.axes.set_yticklabels(minfticks_fft,minor=True,fontsize=7)
            canvas.axes.set_ylim(20,22000)
            canvas.axes.set_title('Spectrogram')
            canvas.axes.set_xlabel('Time (s)')
            canvas.axes.set_ylabel('Frequency (Hz)')
            canvas.fig.tight_layout()

        # toolbar + layout
        toolbar = NavigationToolbar2QT(canvas, self)
        layout = QVBoxLayout()
        layout.addWidget(toolbar)
        layout.addWidget(canvas)
        # widget
        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)
        self.show()      



def main():
    app = QApplication(sys.argv)
    main = GraphWindow()
    main.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()