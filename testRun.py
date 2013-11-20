from Puffin import *
import numpy as np
import sys


def CreateCurves(number):
    Q = []
    for i in range(number):
        Q.append(Precompute(0.01 + i / number * (0.2 - 0.01)))
    return np.array(Q)


def PrintToBed(nucs, fileName):
    try:
        with open(fileName, 'w') as fout:
            pass  # print to file
    except Exception, e:
        print "Smth went wront during outputing to", fileName


def NucsScores(nucs, inputPoints, numPointsTreshold=2, adjustScore=1):
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
                temp = np.concatenate((nuc, [score, np.var(setOfPoints)]))
                count += score
                res.append(temp)
            except Exception:
                print "bad"
                pass
    return np.array(res)


def Run(fileName, Y):
    B = ReadBED(fileName)
    points = B[0].copy()
    shift = np.min(B[0][:, 0])
    points[:, 0] -= shift
    A = BuildSignals(points, Y)
    print 'signals done'
    listNucs = []
    A.shape
    for i in range(0, len(Y) - 1):
        listNucs.append(
            nucdet((A[i] + 1.0) / (A[len(Y) - 1] + 1.0) - 1, 0.0001, A[i]))
        print 'curve ', i, ' is done...'
    #listSize = NucSizeCurves(listNucs, A[0])
    for nucs in listNucs:
        for nuc in nucs:
            nuc[0] += shift
    return listNucs, B[0]


def main():
    fileName = "DATA/test1.bed"
    Y = loadVar('Q_100.var')
    A, inputPoints = Run(fileName, Y)
    print 'Done reading...'
#    saveVar([A, B], fileName + '.var')
#    saveVar([C, inputPoints], 'OUT/curves/' + fileName + '.var')
    nucs = NucPlace(A)
    print 'Placement done'
    D = NucsScores(nucs, inputPoints, 2, 1)
    del A
    del nucs
    D = np.array(D)
    with open(fileName + '.nucs', 'w') as fout:
        for line in D:
            for elem in line:
                print>>fout, elem,
            print>>fout, ""
    print 'Saving done'
#    saveVar(nucs, 'OUT/nucs/' + fileName + '_nucs.var')
#    print 'Saving done'


if __name__ == '__main__':
    main()
