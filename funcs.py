from numpy import log2,ones,convolve

def normalize(data):
    if data.ndim == 1:
        normdata = data/max(abs(data))
    elif data.ndim == 2:
        maxL = max(abs(data[:,0]))
        maxR = max(abs(data[:,1]))
        normdata = data/max([maxL,maxR])
    return normdata

def next_power_of_2(n):
    return 1 << (int(log2(n - 1)) + 1)

def smooth(y, box_pts):
    box = ones(box_pts)/box_pts
    y_smooth = convolve(y, box, mode='same')
    return y_smooth
