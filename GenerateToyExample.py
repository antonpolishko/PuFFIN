import sys
from Puffin import *
import numpy as np


def weighted_values( probabilities, size):
    probabilities = np.array(probabilities)
    bins = np.add.accumulate( (probabilities+0.)/probabilities.sum())
    return np.digitize(random_sample(size), bins)


def main():
    if len(sys.arg) > 1:
        fileNameOut = sys.argv[1]
    else:
        fileNameOut = 'toyExample.bed'
    mu = [] # list to store peak centers
    w = [] # list to store peak weights
    sig = [] # list to store peak variation
    mu.append(300) # location of the main peak
    w.append(10) # weight of the main peak
    sig.append =
    mu.append(225) # location of the first secondary peak
    mu.append(375) # location of the second secondary peak
    w.append(2)
    w.append(1)
    choosePeak = weighted_values(w, 20)









if __name__ == '__main__':
    main()
