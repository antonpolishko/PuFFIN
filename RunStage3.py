from Puffin import *
import numpy as np
import sys
import os.path
import time


def Run(fileName, Y):
    B = ReadBED(fileName)
    listNucs = loadVar("tmpListNucs.var")
    return listNucs, B[0]


def main():
    # fileName = 'DATA/puffin_test_3R_10K.bed'
    if (len(sys.argv) > 1):
        fileName = sys.argv[1]
    else:
        print "please provide the input fileName"
        return 0
    print "Looking for input file"
    if (os.path.isfile(fileName)):
        B = ReadBED(fileName)
        inputPoints = B[0]
    else:
        print "input file doens't exist. Quiting..."
        return 1
    nucsRes = loadVar('tmpNucsRes.var')
    print "Computing nucleosome scores"
    D = NucsScores(nucsRes,inputPoints, -1)
    D = np.array(D)
    D = D[D[:,0].argsort(axis=0)] #sort the nucleosomes according to location
    D = D[:, [0,1,5,6, 4]]
    with open(fileName + '.nucs', 'w') as fout:
        print>>fout, "Location\tPeak_width\tScore\tSTD(fuzziness)\tCurve_level"
        for line in D:
            for elem in line:
                print>>fout, elem,
            print>>fout, ""
    np.savetxt(fileName+'_notAdjusted.nucs',D, fmt='%.2f', delimiter='\t', newline='\n', header='Location\tPeak_width\tScore\tSTD(fuzziness)\tCurve_level')
    print 'Saving done'


if __name__ == '__main__':
    main()
