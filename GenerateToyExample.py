import sys
from Puffin import *
import numpy as np


def weighted_values(probabilities, size):
    probabilities = np.array(probabilities)
    bins = np.add.accumulate((probabilities + 0.) / probabilities.sum())
    return np.digitize(random_sample(size), bins)


def main():
    if len(sys.argv) > 1:
        fileNameOut = sys.argv[1]
    else:
        fileNameOut = 'toyExample.bed'
    mu = []  # list to store peak centers
    w = []  # list to store peak weights
    sig = []  # list to store peak variation
    mu.append(300.)  # location of the main peak
    w.append(80.)  # weight of the main peak
    sig.append(40.)
    mu.append(250.)  # location of the first secondary peak
    w.append(10.)
    sig.append(20.)
    mu.append(350.0)  # location of the second secondary peak
    w.append(10.)
    sig.append(20.)
    choosePeaks = weighted_values(w, 10)
    listOfMidpoints = []  # location of fragment midpoints
    print w, mu, sig
    # for peak in choosePeaks:
    #     x = np.random.normal(mu[peak], sig[peak])
    #     listOfMidpoints.append(x)
    for i in xrange(20):
        x = np.random.normal(mu[0], sig[0])
        listOfMidpoints.append(x)
    # for i in xrange(10):
    #     x = np.random.normal(mu[1], sig[1])
    #     listOfMidpoints.append(x)
    #     x = np.random.normal(mu[2], sig[2])
    #     listOfMidpoints.append(x)
    try:
        f = open(fileNameOut, 'w')
        for loc in listOfMidpoints:
            print >>f, '{0}   {1:d}   {2:d}'.format('chr1', int(loc - 73), int(loc + 73))
        f.close()
    except Exception, e:
        raise e

if __name__ == '__main__':
    main()
