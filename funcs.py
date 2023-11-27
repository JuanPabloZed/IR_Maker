from numpy import log2,ones,convolve

def next_power_of_2(n):
    return 1 << (int(log2(n - 1)) + 1)

def smooth(y, box_pts):
    box = ones(box_pts)/box_pts
    y_smooth = convolve(y, box, mode='same')
    return y_smooth