import matplotlib.pyplot as plt
import numpy as np
from scipy.io.wavfile import read
from funcs import next_power_of_2,smooth

# data 
srate,irfile = read(R'c:\IRs\Test\cab.wav.wav')
npoutfile = np.asarray(irfile)
pad_length = next_power_of_2(next_power_of_2(len(npoutfile)))
padded_npoutfile = np.pad(npoutfile,(0,pad_length-len(npoutfile)),'constant',constant_values=(0,0))
h_panned_npoutfile = padded_npoutfile*np.blackman(pad_length)
fft_outfile = np.fft.rfft(h_panned_npoutfile)
fft_outfile = 20*np.log10(fft_outfile)
fft_outfile_smoothed = smooth(fft_outfile,80)
f = np.fft.rfftfreq(pad_length,1/srate)

i=0
j=-1
while f[i] < 20:
    i+=1
while f[j] > 20000:
    j-=1
limmin = 0.9*min(fft_outfile_smoothed[i:j])
limmax = 1.02*max(fft_outfile_smoothed[i:j])


# figure creation
# x ticks
r = range(1,10)
tix=[]
for power in range(0,5):
    tix += [10**power * x for x in r]

fig = plt.figure()
plt.semilogx(f,fft_outfile_smoothed,'r',lw=0.75) 
plt.xticks(tix,tix,rotation=50,ha='right',fontsize=7)
plt.xlim(20,22000)
plt.ylim(limmin,limmax)
plt.grid(which='both')   
fig.tight_layout()
plt.show()