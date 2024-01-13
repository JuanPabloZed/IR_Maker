import numpy as np
from funcs import next_power_of_2 as npow2,normalize,smooth
from scipy.io.wavfile import read
from scipy.signal import minimum_phase, firwin2
from soundfile import write
import matplotlib.pyplot as plt

def padNhamm(input):
    nfft = 8*npow2(len(input))
    padded_input = np.pad(input,(0,nfft-len(input)),'constant',constant_values=(0,0))
    h_input = padded_input*np.blackman(nfft)
    return h_input


sr,input = read(R'C:/IRs/_Cabs/Homemade/The Extreme/04 MeanFace.wav',mmap=False)
input = np.asarray(input)
newinputlen = npow2(npow2(npow2(len(input))))
input = np.pad(input,(0,newinputlen-len(input)),'constant',constant_values=(0,0))
print(input.shape)
n = len(input)
nfft = 8*npow2(n)

#fft of input for linear phase transform
h_input = padNhamm(input)
inputfft = np.fft.rfft(h_input)
f_input = np.fft.rfftfreq(nfft,1/sr)

#linear phase transform
input_lin = firwin2(len(input) + 1,f_input,inputfft,fs=sr)


# scipy MPT
sc_mpt = minimum_phase(input_lin,n_fft=nfft)
sc_mpt = normalize(sc_mpt)
sc_mpt[1:-1]=sc_mpt[0:-2]
sc_mpt[0]=-0.1

h_sc_mpt = padNhamm(sc_mpt)

# fft curves generation
inputfft = smooth(20*np.log10(abs(inputfft)),100)-np.max(20*np.log10(abs(inputfft)))
sc_mptfft = np.fft.rfft(h_sc_mpt)
sc_mptfft = smooth(20*np.log10(abs(sc_mptfft)),100)-np.max(20*np.log10(abs(sc_mptfft)))
# f axis
f_sc_mpt = np.fft.rfftfreq(8*npow2(len(sc_mpt)),1/sr)

input = input/max(abs(input))

write(r'C:/IRs/_Cabs/Homemade/The Extreme\SC_MPT.wav',sc_mpt,sr,'PCM_24')

plt.Figure()
plt.semilogx(f_input,inputfft,'b-')
plt.semilogx(f_sc_mpt,sc_mptfft,'g-')
plt.xlim(20,22000)
plt.grid()
plt.show()

plt.Figure()
plt.plot(input,'b-')
plt.plot(sc_mpt,'g-')
plt.grid(which='both')
plt.xlim(0,50)
plt.show()