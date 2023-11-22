import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch
from scipy.fft import rfft, rfftfreq
from scipy.io.wavfile import read

def next_power_of_2(n):
    return 1 << (int(np.log2(n - 1)) + 1)

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

path = "c:\IRs\Test\est.wav"

srate, file = read(path)
npoutfile = np.asarray(file)
pad_length = next_power_of_2(next_power_of_2(len(npoutfile)))
padded_npoutfile = np.pad(npoutfile,(0,pad_length-len(npoutfile)),'constant',constant_values=(0,0))
h_panned_npoutfile = padded_npoutfile*np.blackman(pad_length)
# fft_outfile = np.fft.rfft(h_panned_npoutfile)
fft_outfile = np.fft.rfft(h_panned_npoutfile)
f = np.fft.rfftfreq(pad_length,1/srate)
# graph
plt.clf()
plt.semilogx(f,smooth(20*np.log10(fft_outfile),45))
plt.xlim(20,max(f))
plt.title('Spectrum of the IR')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude (dB)')
plt.grid()
plt.show()