from Puffin import *
import numpy as np
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
            newCol = np.ones((len(nucs0), 1), dtype='int')
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

fileName = "DATA/test.bed"
Y = loadVar('Qnew.var')

B = ReadBED(fileName)
points = B[0].copy()
shift = 0  # np.min(B[0][:, 0]) - 1000
points[:, 0] -= shift
print 'Done reading...'
A = BuildSignals(points, Y)
print 'signals done'
listNucs = []
listNucs2 = []
for i in range(0, len(Y) - 1):
    listNucs.append(
        nucdet((A[i] + 1.0) / (A[len(Y) - 1] + 1.0) - 1, 0.001, A[0]))
    listNucs2.append(
        nucdet(A[i], 0.0001, A[i]))
    print 'curve ', i, ' is done...'
    #listSize = NucSizeCurves(listNucs, A[0])

#    saveVar([A, B], fileName + '.var')
#    saveVar([C, inputPoints], 'OUT/curves/' + fileName + '.var')
nucs = NucPlace(listNucs)
nucs2 = NucPlace(listNucs2)

a, b = NucPos(listNucs)

for i in range(99):
    vlines(listNucs[i][:, 0], -i - 1, -i, 'black')
    vlines(listNucs[i][:, 0] - listNucs[i][:, 2], -i - 1, -i, 'yellow')
    vlines(listNucs[i][:, 0] + listNucs[i][:, 3], -i - 1, -i, 'green')

for i in range(99):
    # plot(A[i])
    plot((A[i] + 1.) / (A[-1] + 1.) - 1)
