
# initialization display
print('IR Maker®\nDevelopped by Nathan Zwahlen, Benjamin Quiédeville & Hugo Perrier.')
print('IR Maker is an application made to compute IRs using ESS deconvolution method.')
print('Distributed under XXXX License')
print('Initializing...')

#Homemade Deconvolver
import numpy as np
from scipy.io.wavfile import read, write
from scipy.signal import convolve
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
from warnings import filterwarnings
import os

clear = lambda: os.system('cls')

def normalize(values):
    return values / np.max(values)

def sweepgenerator(f1, f2, T, sr, savepath):
    """
    Generates an Exponential Sine Sweep (ESS) and saves it in a .wav file.
    -----
    INPUTS :
        - f1 (float) : initial frequency for the sweep (Hz)
        - f2 (float) : final frequency for the sweep (Hz)
        - T (int) : duration of the sweep (sec)
        - sr (int) : sample rate of the wav file
        - savepath (str) : path where the ESS will be saved

    OUTPUTS : 
        None, the function only saves the wav file of the ESS.
    """
    R = np.log(f2/f1)   # sweep rate
    time = np.array([i/sr for i in range(T*sr)]) # time array for the graph
    x = np.array([np.sin(2*np.pi*f1*T/R*(np.exp(t/sr*R/T)-1)) for t in range(T*sr)]) # genertion of the sweep
    plt.plot(time, x, linewidth = 0.5) # plot the sweep to show the user the result
    # graph enhancement
    plt.title('Your sine sweep')
    plt.xlabel('Time(s)')
    plt.ylabel('Amplitude')
    print('Done. To continue, please close the graph window.')
    plt.show()
    write(savepath, sr, x) # save the file

    return

def sweep():
    # sweep parameters
    f1 = float(input('Inital frequency for the sweep (Hz) :\n>>> '))
    f2 = float(input('Final frequency for the sweep (Hz) :\n>>> '))
    sr = int(input('Sample rate for the sweep (Hz) :\n>>> '))
    T = int(input('Duration of the sweep (seconds) :\n>>> '))
    savepath = filedialog.asksaveasfilename(title = "Select a location for the ESS file (don't forget the .wav)")
    print('Generating ESS...')
    # sweep generation
    sweepgenerator(f1, f2, T, sr, savepath)
    return

def deconvolver(sweeppath, recordpath, targetname,f1, f2, sr):
    """
    Calulates the IR of a system in which the sine sweep contained in sweeppath has been emitted.
    Uses theoretical inverse convolution of an Exponential Sine Sweep (ESS) to compute the IR.
    Creates a mono or stereo file depending on the format of the recording of the system response contained in recordpath.
    Saves the IR in a wav file named by targetname.
    -----
    INPUTS :
        - sweeppath (str) : path to the sine sweep wav file (MUST BE ESS)
        - recordpath (str) : path to the system's response to the sine sweep
        - targetname (str) : name of the IR wav file (don't forget the '.wav' at the end)
        - f1 (float) : initial frequency for the sweep
        - f2 (float) : final frequency for the sweep
        - sr (int) : IR sample rate (44100Hz, 48000Hz, 88200Hz, 96000Hz, 192000Hz)
    OUTPUTS : 
        - ir (or irL & irR if stereo output), the ndarray containing the IR
    """
    # load files
    rateout, outfile = read(recordpath, mmap=False)

    # sweep read
    ress, ess1 = read(sweeppath)
    # if sweep in stereo format
    if ess1.ndim == 2 :
        ess = ess1[:,0] # converting to mono
    else : 
        ess = ess1  # don't touch anything
 
    # Sweep parameters for theoretical inverse convolution
    R = np.log(f2/f1)   # sweep rate
    T = len(ess)/ress   # duration (s)

    # Theoretical inverse convolution of the sweep
    t = np.array([x/ress for x in range(len(ess))])
    k = np.exp(t*R / T)
    invess = np.flip(ess)/k
    # sweep normalization
    norminvess = normalize(invess)
    # if recording in stereo
    if outfile.ndim == 2:
        # removing the last zeros of the recording
        i = -1
        while outfile[i , 1] == 0:
            i -= 1
        outfileL = outfile[0:len(outfile)+i , 0]
        outfileR = outfile[0:len(outfile)+i , 1]

        # normalization of the recording
        normoutfileL = normalize(outfileL)
        normoutfileR = normalize(outfileR)
        
        # "deconvolving"
        irL = convolve(normoutfileL, norminvess)
        irR = convolve(normoutfileR, norminvess)
        ir = np.column_stack((irL, irR))
        # normalize ir channels
        normir = normalize(ir)
        normirL = normir[:,0]
        normirR = normir[:,1]
        # remove the last zeros on the IR 
        i = -1
        while abs(normirL[i]) <= 0.00005 and abs(normirR[i]) <= 0.00005:
            i -= 1
        normirL = normirL[0:i]
        normirR = normirR[0:i]
        # remove first zeros
        i = 0
        while abs(normirL[i]) <= 0.05 and abs(normirR[i]) <=  0.05 : 
            i += 1
        normirL = normirL[i:len(normirL)]
        normirR = normirR[i:len(normirR)]

        # merge in one stereo wav file
        normir = np.column_stack((normirL, normirR))
        # save wav file
        write(targetname, sr, normir)
        return (normir)
        
    #if recording in mono
    elif outfile.ndim == 1:
        # removing the last zeros of the recording
        i = -1
        while outfile[i] == 0:
            i -= 1
        outfile = outfile[0:len(outfile)+i]

        # normalization of the recording
        normoutfile = normalize(outfile)

        # "deconvolving"
        ir = convolve(normoutfile, norminvess)
        # normalization of the IR
        normir = normalize(ir)

        # removing last zeros of the IR
        i = -1
        while abs(normir[i]) <= 0.0001:
            i -= 1
        normir = normir[0:i]
        # removing first zeros
        i = 0
        while abs(normir[i]) <= 0.1:
            i += 1
        normir = normir[i:len(normir)]
        # save wav file
        write(targetname, rateout, normir)
        return (normir)

def ir():
    # get necessary parameters and file paths
    sweeppath = filedialog.askopenfilename(title = "Select the sine sweep file")
    recordpath = filedialog.askopenfilename(title = "Select the system's response file (cannot contain more than 2 channels)")
    # test if recording is in good format
    test = read(recordpath)[1] 
    if test.ndim == 2 :
        _,c = np.shape(test)
        while c > 2 :
            print('File contains too much channels. Please choose a mono or stereo file only.')
            recordpath = filedialog.askopenfilename(title = "Please choose a mono or stereo file only")
            test = read(recordpath)[1] 
            if test.ndim == 1:
                c=1
            else :
                _,c = np.shape(test)

    targetname = filedialog.asksaveasfilename(title = "Select a location for the IR file") + '.wav'
    f1 = float(input('Initial frequency of the sine sweep (Hz) :\n>>> '))
    f2 = float(input('Final frequency of the sine sweep (Hz) :\n>>> '))
    sr = int(input("IR sample rate (Hz) :\n>>> "))
    # IR format
    #mode = input("'mono' --> Mono IR \n'stereo' --> Stereo IR\n>>> ")
    print('Calculation in progress...') # display for user
    ir = deconvolver(sweeppath, recordpath, targetname, f1, f2, sr)
    # output depending of the format in mode
        # if stereo
    if ir.ndim == 2:
        # stereo IR
        irL = ir[:,0]
        irR = ir[:,1]
        t = [x/sr for x in range(len(irL))] 
        # graph to show the result to the user
        plt.subplot(2,1,1)
        plt.title('Your stereo IR (top : left, down : right) :')
        plt.ylabel('Amplitude')
        plt.plot(t, irL)
        plt.subplot(2,1,2)
        plt.plot(t, irR)
        # if mono
    elif ir.ndim == 1:
        t = [x/sr for x in range(len(ir))] 
        # graph to show the result to the user
        plt.plot(t, ir)
        plt.title('Your mono IR')

    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    print('Done. To continue, please close the graph window.') #display for user
    plt.show()

def main():
    filterwarnings("ignore") #ignore warnings in user display
    clear()
    # Sweep generator mode or IR creating mode
    purpose = input("Sweep generator or IR creator ? \n'sweep' --> Create an ESS\n'ir' --> Create an IR\n'done' --> Close app\n>>> ")
    tk.Tk() # open filedialog window manager (for IDE)
    clear()
    while purpose != 'done':
        # if sweep generator mode
        if purpose == 'sweep':
            sweep()
        #if IR creator mode
        elif purpose == 'ir':    
            ir()
        clear()
        purpose = input("Do you want to do anything else ? \n'sweep' --> create another ESS\n'ir' --> generate an IR \n'done' --> close\n>>> ")
        clear()
    return

if __name__ == '__main__':
    main()