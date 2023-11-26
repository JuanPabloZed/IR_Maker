from scipy import signal
from scipy.io.wavfile import read
import matplotlib.pyplot as plt
import numpy as np
import pyqtgraph

def next_power_of_2(n):
    return 1 << (int(np.log2(n - 1)) + 1)

sr,data = read('c:\IRs\Test\jsp48.wav')
onedata = np.asarray(data[:,0]/max(data[:,0]))

f, t, Sxx = signal.spectrogram(onedata,fs=sr,nfft=len(onedata)//50,nperseg=len(onedata)//400,scaling='spectrum')
Sxx = 20*np.log10(np.matrix.transpose(Sxx))
# Interpret image data as row-major instead of col-major
pyqtgraph.setConfigOptions(imageAxisOrder='row-major')

pyqtgraph.mkQApp()
win = pyqtgraph.GraphicsLayoutWidget()
# A plot area (ViewBox + axes) for displaying the image
p1 = win.addPlot()

# Item for displaying image data
img = pyqtgraph.ImageItem()
img.setImage(Sxx)
tr = pyqtgraph.Qt.QtGui.QTransform()
tr.scale(t[-1] / np.size(Sxx, axis=1), f[-1] / np.size(Sxx, axis=0))  
img.setTransform(tr)
p1.setLimits(xMin=0, xMax=t[-1], yMin=np.log10(20), yMax=np.log10(20000))

hist = pyqtgraph.HistogramLUTItem()
hist.setImageItem(img)
win.addItem(hist)
win.show()
hist.setLevels(np.min(Sxx), np.max(Sxx))
hist.gradient.restoreState(
        {'mode': 'rgb',
         'ticks': [(1.0, (253, 231, 36, 255)),
                   (0.85, (94, 201, 97, 255)),
                   (0.65, (32, 144, 140, 255)),
                   (0.47, (58, 82, 139, 255)),
                   (0.0, (68, 1, 84, 255))]})
ay = p1.getAxis('left')
ay.setLogMode(True)
p1.plot(img)
p1.setLogMode(y=True)
# Add labels to the axis
p1.setLabel('bottom', "Time", units='s')
# If you include the units, Pyqtgraph automatically scales the axis and adjusts the SI prefix (in this case kHz)
p1.setLabel('left', "Frequency", units='Hz')

# Plotting with Matplotlib in comparison
plt.pcolormesh(t, f, Sxx)
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.colorbar()
plt.show()

pyqtgraph.Qt.QtGui.QApplication.instance().exec_()