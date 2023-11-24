import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import read
from scipy.signal import spectrogram

def next_power_of_2(n):
    return 1 << (int(np.log2(n - 1)) + 1)


sr,data = read('c:\IRs\Test\jsp48.wav')
onedata = data[:,0]/max(data[:,0])
npoutfile = np.asarray(onedata)
pad_length = next_power_of_2(next_power_of_2(len(npoutfile)))
padded_npoutfile = np.pad(npoutfile,(0,pad_length-len(npoutfile)),'constant',constant_values=(0,0))
f,t,Sxx = spectrogram(padded_npoutfile,fs=sr,nfft=len(padded_npoutfile)//50,nperseg=len(padded_npoutfile)//400)
plt.yscale('symlog')
plt.pcolormesh(t,f,20*np.log10(Sxx))
cb = plt.colorbar(extend='both')
cb.remove()
plt.clim(-300,-80)
plt.ylim(6,22000)
plt.xlim(0,len(onedata)/sr)
plt.show()