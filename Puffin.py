import numpy as np
import math
import pylab as pl
import pickle
import sys


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


def ReadBAM(fileName):
    try:
        # read the input bam file that contains the whole experiment
        print "Starting process bam file", fileName
        inputBam = Samfile(fileName, 'rb')
        print "file opened..."
        alignments = []
        points = []
        sizes = []
        chrDict = {}
        chrInd = []
        count = 0
        for reference in inputBam.header['SQ']:
            # populate the dictionary of chr_index:"chr reference" with new
            # sequence
            chrDict[reference['SN']] = count
            chrInd.append(reference['SN'])
            count += 1
            point = []
            print reference['SN'], inputBam.count(reference['SN'])
            # add alignments of particular sequence to the list of alignments
            alignmentChr = inputBam.fetch(reference['SN'])
            # populate points with 0-based position of the left-most read
            for read in alignmentChr:
                try:
                    posLeft = 0
                    readSize = 0
                    direction = 0
                    # check whether the alignment is properly paired and
                    # processed only once
                    if read.is_proper_pair and read.is_read1:
                        if read.is_reverse:
                            posLeft = read.mpos
                            readSize = abs(read.isize)
                            direction = -1
                        else:
                            posLeft = read.pos
                            readSize = read.isize
                            direction = 1
                        # filter out the alignments with not proper size
                        # filterSize
                        if readSize > 0 and readSize < 1500 and read.tags[0] == ('XT', 'U'):
                            point.append([posLeft, readSize, direction])
                            sizes.append(readSize)
                except:
                    pass
            # appending a set of points of a processed chromosome to the list
            # of all points
            print "Appeinding chr points to the list...", len(point)
            points.append(np.array(point))
            alignments.append(alignmentChr)
        inputBam.close()
        return points, chrDict, chrInd, sizes
    except:
        print "Can't read input BAM file"


def MakeUnique(data):
    """
    Takes input points (location, size, direction) and leaves only unique points
    """
    try:
        data_new = []
        for chrom in data:
            chrom_new = np.unique(chrom[:, 0] + 1j * chrom[:, 2] * chrom[:, 1])
            chrom_uniq = []
            for row in chrom_new:
                chrom_uniq.append([row.real, abs(row.imag), np.sign(row.imag)])
            data_new.append(np.array(chrom_uniq))
        return data_new
    except Exception, e:
        print "Something went wrong with making points unique"


def ReadInput(fileName):
    """
    Reading file <fileName> that contains pair-end data-points in the format
    (left-most position of the read, left-most location of the mate, size, direction)
    Returns list of points (left-most position of the alignment, size, direction) and list of sizes
    The coordinates are the same 0/1-based as the input in <fileName>
    """
    try:
        # read the input bam file that contains the whole experiment
        with open(fileName, 'r') as inputFile:
            print "Starting process  file", fileName
            sizes = []
            point = []
            for line in inputFile:
                try:
                    posLeft = 0
                    readSize = 0
                    direction = 0
                    read = line.split()
                    if read[3] == "1":
                        posLeft = int(read[0])
                        readSize = abs(int(read[2]))
                        direction = 1
                    else:
                        posLeft = int(read[1])
                        readSize = abs(int(read[2]))
                        direction = -1
                    # filter out the alignments with not proper size
                    # filterSize
                    if readSize > 0 and readSize < 1500:
                        point.append([posLeft, readSize, direction])
                        sizes.append(readSize)
                except:
                    print "Failed to read ", fileName
                    pass
            # appending a set of points of a processed chromosome to the list
            # of all points
        print "Number of points ", len(point)
        return np.array(point), sizes
    except:
        print "Can't read input file", fileName


def importdata(fileName, skipRows=0):
    """
    Reading matrix from a file <fileName>
    It's just a wrapper for mlab.load function
    """
    try:
        A = pl.mlab.load(fileName, skiprows=skipRows)
        return A
    except Exception:
        print "Can't read ", fileName
        return float('nan')


def exportdata(variable, fileName):
    """
    Output matrix in np.array into a file
    """
    if variable.ndim <= 2:
        try:
            with open(fileName, 'w') as fout:
                for line in variable:
                    if line.size() > 1:
                        for element in line:
                            print>>fout, element,
                    else:
                        print>>fout, line,
                    print>> fout
        except Exception:
            print "Failed to output variable to ", fileName
    else:
        print "Input variable is not a matrix (1 or 2 dim)"


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


def ProcessAlignments(alignments, chrDict=None):
    points = []
    for chr in alignments:
        point = []
        for read in chr:
            try:
                d = {}
                for x, y in read.tags():
                    d.setdefault(x, []).append(y)
                if d('XT')[0] == ['U']:
                    point.append([read.pos, read.isize])
                else:
                    print read.tags
            except:
                pass
        points.append(point)
    return points, chrDict


def Clustering(data, epsilon):
    # data = maine[:][0]
    # numberOfPoints = 1000
    # samplePoints = data[random_integers(0, len(data), size=numberOfPoints)]
    X = []
    for i in data:
        X.append([i, 0.15])
    X = np.array(X)
    # samplePoints = \
    # maine[:,random_integers(0,maine.shape[1],size=numberOfPoints)].transpose()
    # X = samplePoints
    import np as np
    from scipy.spatial import distance
    from sklearn.cluster import DBSCAN
    from sklearn import metrics

    D = distance.squareform(distance.pdist(X))
    S = 1 - (D / np.max(D))

    #
    # Compute
    db = DBSCAN(eps=epsilon, min_samples=3).fit(S)
    core_samples = db.core_sample_indices_
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    print 'Estimated number of clusters: %d' % n_clusters_

    import pylab as pl
    from itertools import cycle

    # pl.close('all')
    # pl.figure(1)
    # pl.clf()

    # Black removed and is used for noise instead.
    colors = cycle('bgrcmybgrcmybgrcmybgrcmy')
    for k, col in zip(set(labels), colors):
        if k == -1:
            # Black used for noise.
            col = 'k'
            markersize = 6
        class_members = [index[0] for index in np.argwhere(labels == k)]
        cluster_core_samples = [index for index in core_samples
                                if labels[index] == k]
        for index in class_members:
            x = X[index]
            if index in core_samples and k != -1:
                markersize = 14
            else:
                markersize = 6
            pl.plot(x[0], x[1], 'o', markerfacecolor=col,
                    markeredgecolor='k', markersize=markersize)

    # pl.title('Estimated number of clusters: %d' % n_clusters_)
    # pl.show()
    return db


def printToFile(fileName, data, chr):
    '''
    Create a wig file for the function (data) and file <fileName>
    and specified chromosome (chr)
    '''
    fout = open(fileName, 'w')
    print>>fout, 'fixedStep  chrom=' + chr + ' start=1  step=1\n'
    i = 0
    for item in data:
        print>>fout, float(data[i])
        i = i + 1
    fout.close()


def peakdet(v, delta, x=None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html

    Returns two arrays

    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.

    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.

    """
    maxtab = []
    mintab = []
    if x is None:
        x = np.arange(len(v))
    v = np.asarray(v)
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    lookformax = True
    for i in np.arange(len(v)):
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
        else:
            if this > mn + delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab), np.array(mintab)


def ReadInputPos(inputFileName):
    M = []
    with open(inputFileName, 'r') as file:
        for line in file:
            try:
                words = [int(x) for x in line.strip().split()]
                if words[2] > 10 and words[2] < 1000:
                    M.append(words)
            except Exception, e:
                raise e
            else:
                pass
    m = np.array(M)
    return m


def ReadInputNeg(inputFileName):
    M = []
    with open(inputFileName, 'r') as file:
        for line in file:
            try:
                words = [int(x) for x in line.strip().split()]
                if words[2] > -1000 and words[2] < -10:
                    M.append(words)
            except Exception, e:
                raise e
            else:
                pass
    m = np.array(M)
    temp = np.copy(m[:, 0])
    m[:, 0] = m[:, 1]
    m[:, 1] = temp
    m[:, 2] = abs(m[:, 2])
    return m


def BuildSignal(points, chrLength):
    signal = np.zeros(chrLength)
    for i in range(len(points)):
        lx = min(points[i, 0], points[i, 0] + points[i, 2])
        ly = max(points[i, 0], points[i, 0] + points[i, 2])
        delta = math.abs(lx - ly)
        if (delta > 40 and delta < 1500):
            signal[lx:ly + 1] = signal[lx:ly + 1] + 1 / (points[i, 2] + 0.)
    return signal / len(points) * 1000


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
  #  rawSignal = np.zeros(chrSize)
    if not curves is None:
        numCurves = len(curves)
    else:
        numCurves = 0
    signal = np.zeros([numCurves, chrSize])
#    probSignal = np.zeros([numCurves, chrSize])
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
 #       rawSignal[lx:ly + 1] += 1
        count += 1
        if delta < 1000:
            for i in range(numCurves):
                signal[i][lxS:lyS + 1] += curves[i][x[1]][lxY:1000 - lyY]
 #               if len(x) > 3:
  #                  probSignal[i][lxS:lyS + 1] += 1.0 * curves[i][x[2]][lxY:1000 - lyY] * x[3]
    # Normalization of the curves by number of processed reads
#    rawSignal = rawSignal * 1000000 / count
    print 'Signal populated, now normalazing'
    for i in range(numCurves):
        signal[i] = signal[i] * 1000000. / count
 #       tempSum = probSignal[i].sum()
#        if tempSum > 0:
#            probSignal[i] = probSignal[i] / tempSum
#        probSignal[i] = probSignal[i] * 1000000 / count
    return signal


def CurvesToWig(fileName, data, chrInd):
    '''
    Create a wig file for the curves stored in <data> (data with corresponding chrInd (index:chromosome seq)) and filename (fileName)
    '''
    with open(fileName, 'w') as fout:
        i = 0
        for chrom in data:
            print>>fout, 'fixedStep  chrom=' + \
                chrInd[i][10:] + ' start=1  step=1\n'
            i += 1
            for item in chrom:
                print>>fout, float(item)


 # def BuildMaineSignal(points, chrLength, curves):
 # Build Coverage function for pair-end MAINE coverage
 # by adding Gaussian with area under the curve of size 1
 #    signal = np.zeros(chrLength)
 #    for i in range(len(points)):
 #        mu = math.floor(points[i, 0] + 0.5 * points[i, 2])
 #        w = int(abs(points[i, 2]))
 #        lx = int(max(0, 500 - mu))
 #        ly = int(min(999, chrLength - (mu + 499)))
 #        signal[mu - 500 + lx:mu - 500 + ly + 1] = \
 #            signal[mu - 500 + lx:mu - 500 + ly + 1] + curves[w][lx:ly + 1]
 #    return signal / len(points) * 1000


 # def BuildFaireSignal(points, chrLength):
 #    signal = np.zeros(chrLength)
 #    for i in range(len(points)):
 #        lx = min(points[i, 0], points[i, 0] + points[i, 2])
 #        ly = max(points[i, 0], points[i, 0] + points[i, 2])
 #        signal[lx:ly + 1] = signal[lx:ly + 1] + 1 / (ly - lx + 0.0)
 #    return signal / len(points) * 1000


# def CreateCoverage(inputFileName, curves):
#     a = ReadInputPos(inputFileName)
#     b = ReadInputNeg(inputFileName)
#     chrLength = 1000 + \
#         max([max(a[:, 0]), max(a[:, 1]), max(b[:, 0]), max(b[:, 1])])
#     A = BuildMaineSignal(a, chrLength, curves)
#     B = BuildMaineSignal(b, chrLength, curves)
#     C = (A + B) / 2

#     res = []
#     res.append(C)
#     res.append(A)
#     res.append(B)
#     res.append(a)
#     res.append(b)
#     return res


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
                        nucs.append((center, size))

    return array(nucs)


def PlotNucs(nucs, curve, color='red'):
    """
    Takes array of (center, 1/2*size) nucleosomes and plots a set of segments
    """
    pl.plot(curve)
    for nuc in nucs:
        lx = nuc[0] - nuc[1]
        ly = nuc[0] + nuc[1]
        pl.plot(range(lx, ly), curve[lx:ly], color, linewidth=8, alpha=0.5)
        pl.plot(range(lx, ly), np.arange(lx, ly) * 0, color='black')
        pl.plot([lx, ly], [0, 0], 'o', color='black')


def CurvePurify(threshold, curve, curve2):
    """
    Takes coverage profile curve and curve2 and keeps only locations that exceed threshold in curve
    Everything else is padded with zeros
    """
    nlen = min(len(curve), len(curve2))
    res = curve[:nlen].copy()
    res2 = curve2[:nlen].copy()
    res[res < threshold] = 0
    res2[res < threshold] = 0
    return res, res2


def CurveShrink(threshold, curve, curve2):
    """
    Takes coverage profile curve and curve2 and keeps only locations that exceed threshold in curve
    Everything else is padded with zeros
    """
    res = []
    res2 = []
    nlen = min(len(curve), len(curve2))
    for i in range(nlen):
        if curve[i] > threshold and curve2[i] > threshold:
            res.append(curve[i])
            res2.append(curve2[i])
    return np.array(res), np.array(res2)


def CompareCoverage(maine, faire, threshold=5):
    nlen = max(len(maine), len(faire))
    A = np.zeros(nlen)
    B = np.zeros(nlen)
    A[:len(maine)] = maine
    B[:len(faire)] = faire
    return sum(A > threshold), sum(B > threshold), sum((A > threshold) * 1 * (B > threshold) * 1)


def BuildRawSignal(points):
    """
    Takes points (center, size, direction, ...) and creates coverage function
    that represents simple read pile-up (pair ends treated separately)
    """
    chrLength = max(points[:, 0])
    signal = np.zeros(chrLength)
    for i in range(len(points)):
        lx = min(points[i, 0], points[i, 0] + points[i, 2])
        ly = max(points[i, 0], points[i, 0] + points[i, 2])
        if (ly - lx < 40):
            signal[lx:ly + 1] = signal[lx:ly + 1] + 1
        else:
            signal[lx:lx + 40] = signal[lx:lx + 40] + 1
            signal[ly - 40:ly] = signal[ly - 40:ly] + 1
    return signal


def DepthOfCoverage(points):
    """
    Computes the depth of coverage for a set of points
    how many nucleodite bases are in the reads
    """
    covRaw = 0
    countRaw = 0
    covPair = 0
    countPair = 0
    for i in range(len(points)):
        countPair += 1
        lx = min(points[i, 0], points[i, 0] + points[i, 2])
        ly = max(points[i, 0], points[i, 0] + points[i, 2])
        covPair += int(ly - lx)
        if (ly - lx > 40):
            countRaw += 2
            covRaw += 80
        else:
            countRaw += 1
    return covRaw, covPair, countRaw, countPair


def NucSizeCurves(listOfNucs, rawSignal=None):
    """
    Building curves that represent the distance between adjacent peaks(nucsleosome centers). Each set of nuclesomes in the <listOfNucs> should contain
    sorted list of (center position, size, etc)
    """
    numSets = len(listOfNucs)
    if rawSignal == None:
        chrSize = 0
        for nucs in listOfNucs:
            if chrSize < nucs[-1:, 0]:
                chrSize = nucs[-1:, 0] + nucs[-1:, 1]
    else:
        chrSize = len(rawSignal)
    # baseline for all diff curves
    curves = np.ones([numSets, chrSize], dtype=int) * 200
    i = 0
    for nucs in listOfNucs:
        curve = curves[i]
        lx = 0  # position of the boundary of computed nucSize landscape
        for j in range(len(nucs)):
            proc = nucs[j, 0]  # position of the current processing nucleosome
            if (j < len(nucs) - 1):
                next = nucs[j + 1, 0]  # position of the next nucleosome
            else:
                next = len(curve)

            value = min(next - proc, proc - lx)
            if value < curve[proc - 1]:
                curve[lx:next + 1] = value
            else:
                curve[proc:next + 1] = value
            lx = proc  # move the position of the processed nucSize landscape

        i += 1
    return curves


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


def CurveNucsPurify(nucs, curve):
    """
    Takes coverage profile curve and keeps only locations of nucs (center, 1/2*size)
    Everything else is padded with zeros
    """
    res = np.zeros(len(curve))
    for nuc in nucs:
        res[int(nuc[0] - nuc[1] * 0.8):int(nuc[0] + nuc[1] * 0.8) + 1] = curve[
            int(nuc[0] - nuc[1] * 0.8):int(nuc[0] + nuc[1] * 0.8) + 1]
    return res


def NucsPurify(nucs, inputPoints, numPointsTreshold=2):
    """
    Takes an array of (center, 1/2*size) nucleosomes and set of points and recalculates locations, variation, size and scores
    """
    points = inputPoints[:, 0] + 0.5 * inputPoints[:, 1]
    print len(points)
    res = []
    for nuc in nucs:
        # filter all points that are withing nucleosome boundaries
        tempInd = points[points >= nuc[0] - nuc[1]]
        setOfPoints = tempInd[tempInd <= nuc[0] + nuc[1]]
        if len(setOfPoints) > numPointsTreshold:
            mu = int(np.mean(setOfPoints))
            size = int((max(setOfPoints) - min(setOfPoints)) / 2.)
            if size > 0:
                var = int(np.var(setOfPoints) / size * 2.)
            else:
                var = 0
            score = len(setOfPoints)
            res.append([mu, size, var, score])
    return np.array(res)


def NucsGlobalScore(nucs, inputPoints, numPointsTreshold=2):
    """
    Takes an array of (center, 1/2*size) nucleosomes and set of points and calculates the scores for nucleosomes as a number of points beloning to each nucleosome
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
            score = len(setOfPoints)
            temp = np.concatenate((nuc, [int(score)]))
            count += score
            res.append(temp)
    return np.array(res), count
