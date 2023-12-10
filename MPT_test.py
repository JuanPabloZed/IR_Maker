import numpy as np
from funcs import next_power_of_2 as npow2,normalize,smooth
from scipy.io.wavfile import read
from scipy.signal import minimum_phase
from soundfile import write
import matplotlib.pyplot as plt

sr,input = read(R'c:\IRs\Developping\Cabs\The Extreme\NZ Recorded\IR\M160-HP1-Ctr - IR.wav',mmap=False)

input = np.asarray(input)
print(input.shape)
n = len(input)
nfft = npow2(npow2(n))
halfnfft = int(nfft/2)
kernel = np.zeros((nfft,))
sc_kernel = np.zeros((2*nfft,))
sc_kernel[0:len(input)] = input
kernel[0:len(input)] = input
kcplx = np.zeros((halfnfft+1,))

kcplx = np.fft.rfft(kernel)
arr_magn = np.abs(kcplx)
kcplx = np.log(arr_magn)

kernel = np.fft.irfft(kcplx)
kernel = -1*kernel

kcplx = np.fft.rfft(kernel)

for i in range(0,halfnfft+1):
    real = np.cos(np.imag(kcplx[i])) * arr_magn[i]
    imag = np.sin(np.imag(kcplx[i])) * arr_magn[i]*1j
    kcplx[i] = real + imag

kernel = np.fft.irfft(kcplx)
mpt = kernel[0:len(kernel)//2]

input = normalize(input)
mpt = normalize(mpt)
mpt[3:] = mpt[0:-3]
mpt[0:3] = mpt[[6,5,4]]

atten=2
sc_mpt = minimum_phase((sc_kernel),n_fft=4*nfft)
sc_mpt = smooth(sc_mpt,atten)
sc_mpt = normalize(sc_mpt)

padded_input = np.pad(input,(0,nfft-n),'constant',constant_values=(0,0))
padded_mpt = np.pad(mpt,(0,npow2(npow2(len(mpt)))-len(mpt)),'constant',constant_values=(0,0))
padded_sc_mpt = np.pad(sc_mpt,(0,npow2(npow2(len(sc_mpt)))-len(sc_mpt)),'constant',constant_values=(0,0))
h_sc_mpt = padded_sc_mpt*np.blackman(npow2(npow2(len(sc_mpt))))
h_input = padded_input*np.blackman(nfft)
h_mpt = padded_mpt*np.blackman(npow2(npow2(len(mpt))))


inputfft = np.fft.rfft(h_input)
inputfft = smooth(20*np.log10(abs(inputfft)),45)-np.max(20*np.log10(abs(inputfft)))
mptfft = np.fft.rfft(h_mpt)
mptfft = smooth(20*np.log10(abs(mptfft)),45)-np.max(20*np.log10(abs(mptfft)))
sc_mptfft = np.fft.rfft(h_sc_mpt)
sc_mptfft = smooth(20*np.log10(abs(sc_mptfft)),45)-np.max(20*np.log10(abs(sc_mptfft)))

f_input = np.fft.rfftfreq(nfft,1/sr)
f_mpt = np.fft.rfftfreq(npow2(npow2(len(mpt))),1/sr)
f_sc_mpt = np.fft.rfftfreq(npow2(npow2(len(sc_mpt))),1/sr)



write(r'c:\IRs\Developping\Cabs\The Extreme\NZ Recorded\IR\MPT.wav',mpt,sr,'PCM_24')
write(r'c:\IRs\Developping\Cabs\The Extreme\NZ Recorded\IR\SC_MPT.wav',sc_mpt,sr,'PCM_24')

plt.Figure()
plt.semilogx(f_input,inputfft,'b-')
plt.semilogx(f_mpt,mptfft,'r-')
plt.semilogx(f_sc_mpt,sc_mptfft,'g-')
plt.xlim(20,22000)
plt.grid()
plt.show()

plt.Figure()
plt.plot(input,'b-')
plt.plot(mpt,'r-')
plt.plot(sc_mpt,'g-')
plt.grid(which='both')
plt.xlim(0,15)
plt.show()