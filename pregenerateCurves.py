import Puffin as pf
import sys


def Run(fileName="Q.var"):
    Q = []
    print "Start pregenerating curves"
    for i in range(40):
        print str(i/0.4)+"%"
        Q.append(pf.Precompute(0.05 + (i + 0.) / 100. * 1.45, 1000))
    Q.append(pf.Precompute(1.5, 1000))
    pf.saveVar(Q, fileName)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        Run()
    else:
        Run(sys.argv[1])
