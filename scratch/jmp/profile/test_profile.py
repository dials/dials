from __future__ import division
from dials.algorithms.integration.profile import ProfileFitting2
from scitbx.array_family import flex

def evaluate_gaussian(x, a, x0, sx):

    from math import exp

    assert(len(x) == len(x0))
    assert(len(x) == len(sx))

    g = 0.0
    for xi, x0i, sxi in zip(x, x0, sx):
        g += (xi - x0i)**2 / (2.0 * sxi**2)

    return a * exp(-g)

def gaussian(size, a, x0, sx):

    from scitbx.array_family import flex

    result = flex.double(flex.grid(size))

    index = [0] * len(size)
    while True:
        result[index[::-1]] = evaluate_gaussian(index, a, x0, sx)
        for j in range(len(size)):
            index[j] += 1
            if index[j] < size[j]:
                break
            index[j] = 0
            if j == len(size) - 1:
                return result


#class ProfileFitting2:
#
#    def __init__(self, p, c, b, maxiter=50, eps=1e-3):
#
#        v = b.deep_copy() + 0.0001
#
#        I0 = 0
#        self.niter_ = maxiter
#        for niter in range(maxiter):
#            I = flex.sum((c - b) * p / v) / flex.sum(p*p / v)
#            v = b + I * p
#            if abs(I - I0) < eps:
#                self.niter_ = niter
#                self.eps = abs(I - I0)
#                break
#            I0 = I
#
#        self.intensity_ = I
#        self.variance_ = flex.sum(v)
#
#    def intensity(self):
#        return self.intensity_
#
#    def variance(self):
#        return self.variance_
#
#    def niter(self):
#        return self.niter_
#
#    def eps(self):
#        return self.eps_;


# Create profile
p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
s = flex.sum(p)
p = p / s

# Copy profile
c = gaussian((9, 9, 9), 10, (4, 4, 4), (2, 2, 2))
b = flex.double(flex.grid(9, 9, 9), 0)

# Fit
fit = ProfileFitting2(p, c, b)
I = fit.intensity()
V = fit.variance()
N = fit.niter()

print flex.sum(c), I, V, N
