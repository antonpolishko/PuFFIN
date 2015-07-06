import numpy as np
import math
import pylab as pl
import pickle
import sys


def NucPos(listNucsIn):
    nucRes = np.array([0, 0, 0, 0, 0])
    nucChild = []
    numLevels = len(listNucsIn)
    listNucs = list(listNucsIn)
    flag = True
    justAdded = 0
    lastAddedNuc = 0
    while flag:
        flag = False
        justAdded = 0
        for curLevel in range(numLevels):
            nucs0 = listNucs[curLevel]
            newCol = np.ones((len(nucs0), 1), dtype='int') * curLevel
            nucs0 = np.hstack((nucs0, newCol))
            if len(nucs0) > 0:
                # add first and last nucleosomes
                nucs = np.vstack(
                    (np.array([0, 0, 0, 0, 0]), nucs0, [nucs0[-1, 0] + 200, 0, 0, 0, 0]))
                dist = np.roll(nucs[:, 0], -1, 0) - nucs[:, 0]
                ind = (dist > 147) * (np.roll(dist, 1, 0) > 147)
                ind1 = np.delete(ind, 0)
                ind = np.delete(ind1, len(ind1) - 1)
                nucs = nucs[1:-1, :]
                justAdded = ind.sum()
                if justAdded > 0:
                    flag = True
                    nucRes = np.vstack((nucRes, nucs[ind, :]))
                    break
        if justAdded > 0:
                # Remove nucleosomes that overlap with solution from the
                # list
            for l in range(justAdded):
                lastAddedNuc += 1
                nucToProcess = nucRes[lastAddedNuc][0]
                nucOver = [0, 0, 0, 0]
                for i in range(curLevel):
                    nucs = listNucs[i]
                    if len(nucs) > 0:
                        ind = abs(nucs[:, 0] - nucToProcess) < 147
                        nucOver = np.vstack((nucOver, nucs[ind, :]))
                        listNucs[i] = nucs[~ind, :]
                for i in range(curLevel, numLevels):
                    nucs = listNucs[i]
                    if len(nucs) > 0:
                        ind = abs(nucs[:, 0] - nucToProcess) < 147
                        listNucs[i] = nucs[~ind, :]
                nucChild.append(nucOver[1:])
    return nucRes[1:, :], nucChild


def ReadBED(fileName):
    """
    Reading BED file <fileName> that contains pair-end data-points in the format
    (chromosome name, left-most location of the mate, right-most location of the mate, mapping score)
    Returns list of points (left-most position of the alignment, size, direction) and list of sizes
    The coordinates are the same 0/1-based as the input in <fileName>
    """
    try:
        # read the input bam file that contains the whole experiment
        with open(fileName, 'rU') as inputFile:
            print "Starting process  file", fileName
            sizes = []
            point = []
            for line in inputFile:
                read = line.split()
                direction = 1
                posLeft = int(read[1])
                readSize = abs(int(read[2]) - int(read[1]))
                flag = False
                if (len(read) > 3):
                    readScore = int(read[3])
                    flag = True
                # filter out reads with bad mapping score
                if (flag) and (readScore > 10):
                    point.append([posLeft, readSize, direction])
                    sizes.append(readSize)
                else:
                    if (~flag):
                        point.append([posLeft, readSize, direction])
                        sizes.append(readSize)
        print "Number of points ", len(point)
        return np.array(point), sizes
    except:
        print "Can't read input file", fileName


def saveVar(var, fileName):
    """
    just a one-liner style wrapper for pickle.dump module
    uses binary files and -1 protocol
    """
    from pickle import dump
    try:
        with open(fileName, 'wb') as fout:
            dump(var, fout, -1)
            print "saving successful"
    except Exception:
        print "Smth went wrong during dumping..."


def loadVar(fileName):
    """
    A one-liner style wrapper for pickle.load
    reads file <fileName> as binary
    """
    try:
        with open(fileName, 'rb') as fin:
            A = pickle.load(fin)
            print "reading done...", fileName
        return A
    except Exception:
        print "Loading failed ...", fileName
        return float('nan')


def Precompute(alpha, sizeLen):
    # Create matrix that stores precomputed templates
    B = np.zeros((1000, sizeLen))
    alpha = float(alpha)

    def _Gauss(x, mu, sigma):
    # Gauss function at point x with parameters (mu,sigma^2)
        x = float(x)
        mu = float(mu)
        return math.exp(-(x - mu) ** 2 / sigma ** 2)

    for i in range(10, 1000):
        vec = np.zeros(sizeLen)
        sigma = alpha * i
        mu = sizeLen / 2
        coef = 1. / \
            (sigma * math.sqrt(math.pi) * math.erf(mu / sigma)
             - 2 * mu * _Gauss(0, mu, sigma))
        shift = _Gauss(0, mu, sigma)
        for j in range(sizeLen):
            vec[j] = vec[j] + coef * (_Gauss(float(j), mu, sigma) - shift)
        B[i, :] = vec
    return B


def BuildSignals(dataChr, curves=None):
    '''
    takes points for given chromosome (leftPos, size, direction, [probability])
    and creates raw signal, signal with curves Y, signal with curved normalized by probability
    '''
    chrSize = int(max(dataChr[:, 0])) + 2000  # estimate chromosome size
    if not curves is None:
        numCurves = len(curves)
    else:
        numCurves = 0
    signal = np.zeros([numCurves, chrSize])
    count = 0  # number of processed reads
    for x in dataChr:
        mu = int(x[0] + 0.5 * x[1])
        lx = int(max(x[0], 0))
        ly = int(min(x[0] + x[1], chrSize))
        lxS = max(mu - 500, 0)
        lxY = - (mu - 500) + lxS
        lyS = min(mu + 499, chrSize)
        lyY = int(mu + 499 - lyS)
        delta = ly - lx
        count += 1
        if delta < 1000:
            for i in range(numCurves):
                signal[i][lxS:lyS + 1] += curves[i][x[1]][lxY:1000 - lyY]
    print 'Signal populated, now normalazing'
    for i in range(numCurves):
        signal[i] = signal[i] * 1000000. / count
    return signal


def nucdet(curve, delta, curveOrig):
    """
    Takes function and places nucleosomes at peak with size calculated
    as distance to the closest deep
    //the code is based on peakdetction function https://gist.github.com/250860
    """
    from numpy import NaN, Inf, arange, isscalar, array, asarray
    maxtab = []
    mintab = []
    nucs = []

    x = arange(len(curve))

    v = asarray(curve)

    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    center = NaN
    left = 0
    right = NaN

    lookformax = True
    lookforbound = False

    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx - delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
                center = mxpos
                lookforbound = True
        else:
            if this > mn + delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
                if lookforbound:
                    sizeL = center - left
                    sizeR = mnpos - center
                    size = min(sizeL, sizeR)
                    left = mnpos
                    lookforbound = False
                    if (curveOrig[center] > delta) and (curve[center] > -0.5):
                        nucs.append((center, size, sizeL, sizeR))

    return array(nucs)



def NucsScores(nucs, inputPoints, numPointsTreshold=0, adjustScore=1):
    """
    Takes an array of (center, 1/2*size) nucleosomes and set of points
    and calculates the scores for nucleosomes as a number of points beloning to each nucleosome
    """
    points = inputPoints[:, 0] + 0.5 * inputPoints[:, 1]
    print len(points)
    res = []
    count = 0
    for nuc in nucs:
        # filter all points that are withing nucleosome boundaries
        tempInd = points[points >= nuc[0] - nuc[1]]
        setOfPoints = tempInd[tempInd <= nuc[0] + nuc[1]]
        if len(setOfPoints) > numPointsTreshold:
            try:
                score = len(setOfPoints) * adjustScore
                temp = np.concatenate((nuc, [score, np.std(setOfPoints)]))
                count += score
                res.append(temp)
            except Exception:
                print "bad"
                pass
    return np.array(res)


def NucsAdjust(nucs, inputPoints):
    """
    Takes an array of (center, 1/2*size) nucleosomes and set of points
    and recalculates the centers for nucleosomes as a centroid of corresponding points
    """
    points = inputPoints[:, 0] + 0.5 * inputPoints[:, 1]
    res = []
    res2 = []
    for nuc in nucs:
        # filter all points that are withing nucleosome boundaries
        tempInd = points[points >= nuc[0] - nuc[1]]
        setOfPoints = tempInd[tempInd <= nuc[0] + nuc[1]]
        tempInd2 = points[points >= nuc[0] - 73]
        setOfPoints2 = tempInd2[tempInd2 <= nuc[0] + 73]
        if len(setOfPoints) > 0:
            try:
                temp = nuc.copy();
                temp[0] = np.mean(setOfPoints)
                res.append(temp)
            except Exception:
                print "bad"
                pass
        if len(setOfPoints2) > 0:
            try:
                temp = nuc.copy();
                temp[0] = np.mean(setOfPoints2)
                res2.append(temp)
            except Exception:
                print "bad"
                pass
    return np.array(res), np.array(res2)


def NucPlace(listOfNucsIn, listOfSizes=None):
    """
    Picking set of non overlapping nucleosomes. Each set of nuclesomes in the <listOfNucs> should contain
    sorted list of (center position, size, etc)

    The list itself should contain sorted list of sets with "thinner"-to-"fatter" nucleosomes
    """
    listOfNucs = list(listOfNucsIn)
    if (listOfSizes == None):
        listOfSizes = NucSizeCurves(listOfNucs)
    chrSize = len(listOfSizes[0])
    Nucleos = []
    curLevel = 0  # current level of which curve produces properly spaced peaks
    # current position on the genome, to the left everything is already
    # processed
    curPosition = 0
    flagRun = True  # flag to run the main cycle
    while flagRun:
        # remove nuclesomes to the left of the processed boundary
        for i in range(len(listOfNucs)):
            nucs = listOfNucs[i]
            while (len(listOfNucs[i]) > 0 and listOfNucs[i][0][0] <= curPosition):
                listOfNucs[i] = np.delete(listOfNucs[i], (0), axis=0)
        # pick level that satisfies
        level = 0
        for nucs in listOfNucs:
            # check whether potential nucleosome is within 146bp of it's
            # neighbors
            if (len(nucs) > 0 and min(listOfSizes[level][nucs[0][0] - 2:nucs[0][0] + 2] > 146)):
                break
            else:
                level += 1
        # level picked, so we already now the nucleosome
        if (len(nucs) == 0):
            flagRun = False
            break
        mu = nucs[0][0]
        size = nucs[0][1]
        if len(Nucleos) > 0:
            # check whether previous nucleosome is overlapping with the
            # candidate one
            if mu - Nucleos[-1][0] < 140:
                if curLevel < level:
                    # delete prev nucleosome since is "covered" by the
                    # candidate, add candidate to the list
                    Nucleos.pop()
                    Nucleos.append([mu, size, level])
                else:
                    # nothing to do here,  prev nucleosome already covering the
                    # candidate
                    pass
            else:
                Nucleos.append([mu, size, level])

        else:
            # no nucleosomes -> just add candidate to the list
            Nucleos.append([mu, size, level])
        # update current position
        curPosition = mu + 140  # we add some offset
        # print curPosition, len(Nucleos), listOfNucs[0][0], listOfNucs[15][0]
        curLevel = level
        if curPosition < chrSize or len(listOfNucs[0]) > 0:
            flagRun = True
        else:
            flagRun = False

    return Nucleos