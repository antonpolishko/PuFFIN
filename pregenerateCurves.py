import Puffin as pf
import sys


def Run(fileName="Q.var", num_curves=40):
    Q = []
    print("Start pregenerating curves")
    for i in range(num_curves):
        print( str(i/num_curves*100)+"%" )
        Q.append(pf.Precompute(0.03 + (i + 0.) / 100. * 1.45, 1000))
    Q.append(pf.Precompute(1.5, 1000))
    pf.saveVar(Q, fileName)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        Run()
    else:
        Run(sys.argv[1])
