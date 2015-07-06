from Puffin import *
import numpy as np
import sys
import os.path


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
            nucdet(np.log((A[i] + 1.0) / (A[len(Y) - 1] + 1.0)), 0.0001, A[i]))
        print 'curve', i, str(i/0.4)+'% is done...'
    #listSize = NucSizeCurves(listNucs, A[0])
    for nucs in listNucs:
        for nuc in nucs:
            nuc[0] += shift
    return listNucs, B[0]


def main():
    fileName = sys.argv[1]
    if (os.path.isfile('Q.var')):
        Y = loadVar('Q.var')
    else:
        import pregenerateCurves
        print "No curves file found, trying to create one (this step is suppose to take place only once)"
        pregenerateCurves.Run()
        if (os.path.isfile('Q.var')):
            Y = loadVar('Q.var')
        else:
            print "Something went wrong with loading curve file. Try running >python pregenerateCurve.py"        
    print "Looking for input file"
    if (os.path.isfile(fileName)):
        A, inputPoints = Run(fileName, Y)
    else:
        print "input file doens't exist. Quiting..."
        return 1
    print 'Done reading...'
    C = []
    for line in A:
        b, c = NucsAdjust(line, inputPoints)
        C.append(c)
    nucsRes, child = NucPos(C)
    print 'Placement done'
    D = NucsScores(nucsRes,inputPoints, -1, 1)
    del A
    D = np.array(D)
    D = D[:, [0,1,5,6, 4]]
    # with open(fileName + '.nucs', 'w') as fout:
    #     print>>fout, "Location\tPeak_width\tScore\tSTD(fuzziness)\tCurve_level"
    #     for line in D:
    #         for elem in line:
    #             print>>fout, elem,
    #         print>>fout, ""
    np.savetxt(fileName+'.nucs',D, fmt='%.2f', delimiter='\t', newline='\n', header='Location\tPeak_width\tScore\tSTD(fuzziness)\tCurve_level')
    print 'Saving done'


if __name__ == '__main__':
    main()
