import Puffin as pf
import sys


def main():
    Q = []
    for i in range(100):
        Q.append(pf.Precompute(0.05 + (i + 0.) / 100. * 1.95))
    pf.saveVar(Q, 'Q_100.var')


if __name__ == '__main__':
    main()
