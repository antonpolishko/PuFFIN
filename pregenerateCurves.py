import Puffin as pf
import sys


def main():
    if sys.argv[1] is None:
        fileName = "Q.var"
    else:
        fileName = sys.argv[1]
    Q = []
    for i in range(100):
        print i
        Q.append(pf.Precompute(0.05 + (i + 0.) / 100. * 1.45, 1000))
    pf.saveVar(Q, fileName)


if __name__ == '__main__':
    main()
