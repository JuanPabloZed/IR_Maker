import numpy as np
from funcs import next_power_of_2 as npow2,normalize,smooth
from scipy.io.wavfile import read,write
from scipy.signal import minimum_phase
import matplotlib.pyplot as plt

sr,input = read(R'c:\IRs\_Cabs\KOHLE\Rainbows & Chainsaws\48\Rainbows\05 Bulldozer Death.wav',mmap=False)

input = np.asarray(input)

n = len(input)
nfft = npow2(npow2(n))

mptransformed = np.asarray(minimum_phase(input,method='homomorphic',n_fft=nfft))
mptransformed[0:9] = [mptransformed[i]*1.5-(0.05*i) for i in range(0,9)]
mptransformed = normalize(mptransformed)
mptransformed[1:-1]=mptransformed[0:-2]
mptransformed[0]=0
Smptransformed = smooth(mptransformed,4)

padded_input = np.pad(input,(0,nfft-n),'constant',constant_values=(0,0))
padded_mpt = np.pad(mptransformed,(0,npow2(npow2(len(mptransformed)))-len(mptransformed)),'constant',constant_values=(0,0))
h_input = padded_input*np.blackman(nfft)
h_mpt = padded_mpt*np.blackman(npow2(npow2(len(mptransformed))))

inputfft = np.fft.rfft(h_input)
inputfft = smooth(20*np.log10(abs(inputfft)),45)-np.max(20*np.log10(abs(inputfft)))
mptfft = np.fft.rfft(h_mpt)
mptfft = smooth(20*np.log10(abs(mptfft)),45)-np.max(20*np.log10(abs(mptfft)))

f_input = np.fft.rfftfreq(nfft,1/sr)
f_mpt = np.fft.rfftfreq(npow2(npow2(len(mptransformed))),1/sr)

input = normalize(input)

write(r'c:\IRs\Test\MPT.wav',sr,normalize(mptransformed))
write(r'c:\IRs\Test\MPT2.wav',sr,normalize(Smptransformed))
plt.Figure()
plt.semilogx(f_input,inputfft,'b-')
plt.semilogx(f_mpt,mptfft,'r-')
plt.grid()
plt.show()

plt.Figure()
plt.plot(input,'b-')
plt.plot(mptransformed,'r-')
plt.grid()
plt.xlim(0,140)
plt.show()