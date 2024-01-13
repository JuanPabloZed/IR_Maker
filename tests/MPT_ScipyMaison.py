import numpy as np
from funcs import next_power_of_2 as npow2,normalize,smooth
from scipy.io.wavfile import read
from scipy.signal import minimum_phase
from soundfile import write
import matplotlib.pyplot as plt

sr,input = read(R'c:\IRs\Developping\Cabs\The Extreme\NZ Recorded\IR\M160-HP3-Bal - IR.wav',mmap=False) 
input = normalize(np.asarray(input))
n = len(input)
n_fft = npow2(npow2(n))

# scipy engine
# zero-pad; calculate the DFT
h_temp = np.abs(np.fft.rfft(input, n_fft))
# take 0.25*log(|H|**2) = 0.5*log(|H|)
h_temp += 1e-7 * h_temp[h_temp > 0].min()  # don't let log blow up
np.log(h_temp, out=h_temp)
h_temp *= 0.5
# IDFT
h_temp = np.fft.irfft(h_temp).real
# multiply pointwise by the homomorphic filter
# lmin[n] = 2u[n] - d[n]
win = np.zeros(n_fft)
win[0] = 1
stop = (len(input) + 1) // 2
win[1:stop] = 2
if len(input) % 2:
    win[stop] = 1
h_temp *= win
h_temp = np.fft.irfft(np.exp(np.fft.rfft(h_temp)))
h_minimum = h_temp.real
n_out = n//2 + len(input) % 2
mpt = normalize(h_minimum[:n_out])

write(r'c:\IRs\Developping\Cabs\The Extreme\NZ Recorded\IR\M160-HP3-Bal_MPT.wav',mpt,sr,'PCM_24')


# FFT plot data preparation
padded_input = np.pad(input,(0,n_fft-n),'constant',constant_values=(0,0))
padded_mpt = np.pad(mpt,(0,npow2(npow2(len(mpt)))-len(mpt)),'constant',constant_values=(0,0))

h_input = padded_input*np.blackman(n_fft)
h_mpt = padded_mpt*np.blackman(npow2(npow2(len(mpt))))

inputfft = np.fft.rfft(h_input)
inputfft = smooth(20*np.log10(abs(inputfft)),45)-np.max(20*np.log10(abs(inputfft)))
mptfft = np.fft.rfft(h_mpt)
mptfft = smooth(20*np.log10(abs(mptfft)),45)-np.max(20*np.log10(abs(mptfft)))

f_input = np.fft.rfftfreq(n_fft,1/sr)
f_mpt = np.fft.rfftfreq(npow2(npow2(len(mpt))),1/sr)

# plot
plt.figure()
plt.semilogx(f_input,inputfft,'b-')
plt.semilogx(f_mpt,mptfft,'r-')
plt.grid(which='both')
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(input,'b-')
plt.plot(mpt,'r-')
plt.grid(which='both')
plt.tight_layout()
plt.xlim(0,100)
plt.show()