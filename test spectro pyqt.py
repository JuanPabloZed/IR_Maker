from scipy.signal import spectrogram
from scipy.io.wavfile import read
import numpy as np
import pyqtgraph

def next_power_of_2(n):
    return 1 << (int(np.log2(n - 1)) + 1)

sr,data = read('c:\IRs\Test\jsp48.wav')
onedata = np.asarray(data[:,0]/max(data[:,0]))
f,t,Sxx = spectrogram(onedata,fs=sr,nfft=len(onedata)//50,nperseg=len(onedata)//400)
Sxx = 20*np.log10(Sxx)
# Interpret image data as row-major instead of col-major
pyqtgraph.setConfigOptions(imageAxisOrder='row-major')

pyqtgraph.mkQApp()
win = pyqtgraph.GraphicsLayoutWidget()
# A plot area (ViewBox + axes) for displaying the image
p1 = win.addPlot()

# Item for displaying image data
img = pyqtgraph.ImageItem()
p1.addItem(img)
# Add a histogram with which to control the gradient of the image
hist = pyqtgraph.HistogramLUTItem()
# Link the histogram to the image
hist.setImageItem(img)
# If you don't add the histogram to the window, it stays invisible, but I find it useful.
# win.addItem(hist)
# Show the window
win.show()
# Fit the min and max levels of the histogram to the data available
hist.setLevels(np.min(Sxx), np.max(Sxx))
# This gradient is roughly comparable to the gradient used by Matplotlib
# You can adjust it and then save it using hist.gradient.saveState()
hist.gradient.restoreState(
        {'mode': 'rgb',
         'ticks': [(1.0, (253, 231, 36, 255)),
                   (0.85, (94, 201, 97, 255)),
                   (0.65, (32, 144, 140, 255)),
                   (0.47, (58, 82, 139, 255)),
                   (0.0, (68, 1, 84, 255))]})
# Sxx contains the amplitude for each pixel
img.setImage(Sxx)
tr = pyqtgraph.Qt.QtGui.QTransform()
tr.scale(t[-1] / np.size(Sxx, axis=1), f[-1] / np.size(Sxx, axis=0))  
img.setTransform(tr)
# Scale the X and Y Axis to time and frequency (standard is pixels)
# img.scale(t[-1]/Sxx.shape[0], f[-1]/Sxx.shape[1])
# Limit panning/zooming to the spectrogram
p1.setLimits(xMin=0, xMax=t[-1], yMin=0, yMax=f[-1])
# Add labels to the axis
p1.setLabel('bottom', "Time")
# If you include the units, Pyqtgraph automatically scales the axis and adjusts the SI prefix (in this case kHz)
p1.setLabel('left', "Frequency")

pyqtgraph.Qt.QtGui.QApplication.instance().exec_()