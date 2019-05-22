from __future__ import print_function

from matplotlib import pylab

if __name__ == "__main__":

    import sys

    with open(sys.argv[1]) as infile:
        X = []
        Y = []
        for line in infile.readlines():
            t = line.split()
            X.append(int(t[0]))
            Y.append(float(t[1]))

    print(min(Y), min(X))

    ax = pylab.scatter(X, Y)
    pylab.xlim((1, 1500))
    pylab.xscale("log")
    pylab.ylim((1e-6, 1e-1))
    pylab.yscale("log")
    pylab.xlabel("Number of observations")
    pylab.ylabel("KL Divergence")

    pylab.show()
