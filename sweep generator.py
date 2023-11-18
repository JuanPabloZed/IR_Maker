from scipy.io.wavfile import write
import numpy as np
import matplotlib.pyplot as plt



def sweepgenerator(f1, f2, T, sr, savepath):
    """
    Generates an Exponential Sine Sweep (ESS) and saves it in a .wav file.
    -----
    INPUTS :
        - f1 (float) : initial frequency for the sweep (Hz)
        - f2 (float) : final frequency for the sweep (Hz)
        - T (float) : duration of the sweep (sec)
        - sr (int) : sample rate of the wav file
        - savepath (str) : path where the ESS will be saved

    OUTPUTS : 
        None, the function only saves the wav file of the ESS.
    """
    R = np.log(f2/f1)   # sweep rate
    time = np.array([i/sr for i in range(T*sr)]) # time array for the graph
    x = np.array([np.sin(2*np.pi*f1*T/R*(np.exp(t/sr*R/T)-1)) for t in range(T*sr)]) # genertion of the sweep
    plt.plot(time, x) # plot the sweep to show the user the result
    plt.show()
    write(savepath, sr, x) # save the file

    return