import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import read
from scipy.signal import spectrogram
def next_power_of_2(n):
    return 1 << (int(np.log2(n - 1)) + 1)


sr,data = read('c:\IRs\Test\e.wav.wav')
onedata = data[:,0]/max(data[:,0])
npoutfile = np.asarray(onedata)
pad_length = next_power_of_2(next_power_of_2(len(npoutfile)))
padded_npoutfile = np.pad(npoutfile,(0,pad_length-len(npoutfile)),'constant',constant_values=(0,0))
f,t,Sxx = spectrogram(padded_npoutfile,fs=sr,nfft=256,noverlap=128,scaling='density')
plt.pcolormesh(t,f,20*np.log10(Sxx))
plt.show()