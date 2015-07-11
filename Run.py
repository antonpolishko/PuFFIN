from Puffin import *
import numpy as np
import sys
import os.path
import time


def Run(fileName, Y):
    B = ReadBED(fileName)
    points = B[0].copy()
    shift = np.min(B[0][:, 0])
    points[:, 0] -= shift
    A = BuildSignals(points, Y)
    print 'signals done'
    listNucs = []
    A.shape
    print "start placing nucleosome candidates on computed curves..."
    for i in range(0, len(Y) - 1):
        listNucs.append(
            nucdet(np.log((A[i] + 1.0) / (A[len(Y) - 1] + 1.0)), 0.0001, A[i]))
        print 'curve', i, str(int((i+0.0)/len(Y)*100.0))+'% is done...'
    #listSize = NucSizeCurves(listNucs, A[0])
    for nucs in listNucs:
        for nuc in nucs:
            nuc[0] += shift
    print "nucleosome candidates are done, dumping into tmpListNucs file"
    saveVar(listNucs, "tmpListNucs.var")
    return listNucs, B[0]


def main():
    # fileName = 'DATA/puffin_test_3R_10K.bed'
    if (len(sys.argv) > 1):
        fileName = sys.argv[1]
    else:
        print "please provide the input fileName"
        return 0
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
    # Starting extra step to make more precise locations for nucleosome candidates
    # If this step takes to long, comment this section out
    print "Start adjusting nucleosome candidates, if it takes too long, comment out the code..."
    C = []
    i = 0
    for line in A:
        i += 1
        start = time.time()
        c = NucsAdjust(line, inputPoints)
        C.append(c)
        fin = time.time()
        print 'curve', i, str(int((i+0.0)/len(A)*100.0))+'% is done...'
        print "it took " + str(int(fin - start)) + "sec"
    A = C
    saveVar(C, "tmpAdjustedNucs.var")
    # End of optional section
    print "Starting picking final nucleosomes from candidates"
    start = time.time()
    nucsRes, child = NucPos(A)
    print 'Placement done'
    fin = time.time()
    print "it took " + str(int(fin - start)) + "sec"
    saveVar(nucsRes, 'tmpNucsRes.var')
    print "Computing nucleosome scores"
    D = NucsScores(nucsRes,inputPoints, -1)
    D = np.array(D)
    D = D[D[:,0].argsort(axis=0)] #sort the nucleosomes according to location
    D = D[:, [0,1,5,6, 4]]
    # with open(fileName + '.nucs', 'w') as fout:
    #     print>>fout, "Location\tPeak_width\tScore\tSTD(fuzziness)\tCurve_level"
    #     for line in D:
    #         for elem in line:
    #             print>>fout, elem,
    #         print>>fout, ""
    np.savetxt(fileName+'.nucs',D, fmt='%.2f', delimiter='\t', newline='\n', header='Location\tPeak_width\tScore\tSTD(fuzziness)\tCurve_level')
    print 'Saving done'
    return D, inputPoints, A, C


if __name__ == '__main__':
    main()
