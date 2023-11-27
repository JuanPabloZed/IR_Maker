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

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow,self).__init__(parent)
        self.setWindowTitle('Test')
        # MainWindow.resize(self,1000,460)
        sc = MplCanvas(self, width=15, height=4, dpi=100)

        # self.graph = pyqtgraph.PlotWidget(self)
        # self.graph.setGeometry(30,30,940,400)
        self.plotgraph(sc)
    
    def plotgraph(self,canvas):
        # data pour plot spectro
        sr,data = read(R'c:\IRs\_Verbs\Kohle Reverberant Room\KM183 Tight.wav')
        onedata = np.asarray(data[:,0]/max(data[:,0]))
        f_spectro, t_spectro, Sxx = signal.spectrogram(onedata,fs=sr,nfft=len(onedata)//50,nperseg=len(onedata)//400,scaling='spectrum')
        Sxx = 20*np.log10(np.matrix.transpose(Sxx))

        # data pour plot normal
        srate,irfile = read(R'c:\IRs\Test\cab.wav.wav')
        npoutfile = np.asarray(irfile)
        pad_length = next_power_of_2(next_power_of_2(len(npoutfile)))
        padded_npoutfile = np.pad(npoutfile,(0,pad_length-len(npoutfile)),'constant',constant_values=(0,0))
        h_panned_npoutfile = padded_npoutfile*np.blackman(pad_length)
        fft_outfile = 20*np.log10(abs(np.fft.rfft(h_panned_npoutfile)))
        fft_outfile = smooth(fft_outfile,80) - np.max(fft_outfile)
        f_fft = np.fft.rfftfreq(pad_length,1/srate)

        # display data
        # x ticks for fft
        majxticks_fft=[10**power for power in range(0,5)]
        minxticks_fft=[]
        for power in range(0,5):
            minxticks_fft += [10**power * x for x in range(2,10)]
        # plot
        canvas.axes.semilogx(f_fft,fft_outfile,lw=0.7)
        canvas.axes.grid(which='both')
        canvas.axes.set_xticklabels(majxticks_fft,minor=False,rotation=50)
        canvas.axes.set_xticklabels(minxticks_fft,minor=True,rotation=50)
        canvas.axes.set_xlim(20,22000)
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
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()