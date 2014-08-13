


# def test(C, B, S):

#   Fsr = sum([s * (c - mu) / mu for s, c, mu in zip(S, C, MU)])
#   Fs2 = sum([c * (s / mu)**2 for s, c, mu in zip(S, C, MU)])
#   Fbr = sum([b * (c - mu) / mu for b, c, mu in zip(B, C, MU)])
#   Fb2 = sum([c * (b / mu)**2 for b, c, mu in zip(B, C, MU)])
#   sigma_Sa = (sqrt(Fsr**2 + Fs2) + Fsr) / Fs2
#   sigma_Ba = (sqrt(Fbr**2 + Fb2) + Fbr) / Fb2

# def test(C, Bs, Ss):
#   R = [s / b for s, b in zip(Ss, Bs)]
#   S0 = sum([c*(r**0) for c, r in zip(C, R)])
#   S1 = sum([c*(r**1) for c, r in zip(C, R)])
#   S2 = sum([c*(r**2) for c, r in zip(C, R)])
#   S3 = sum([c*(r**3) for c, r in zip(C, R)])
#   Fb = sum(Bs)
#   Fs = sum(Ss)

#   a = Fb*S3 - Fs*S2
#   b = Fb*S2 - Fs*S1
#   c = Fb*S1 - Fs*S0
#   Ra = (c / b) * (1 + a*c/(b**2))
#   Ba = (S0 - Ra*S1+Ra*Ra*S2) / Fb
#   Sa = Ra * Ba

#   print "S0: ", S0
#   print "S1: ", S1
#   print "S2: ", S2
#   print "S3: ", S3
#   print "Fb: ", Fb
#   print "Fs: ", Fs
#   print "Ra: ", Ra
#   print "Ba: ", Ba
#   print "Sa: ", Sa

def iteration(c, b, s, B, S):
  from scitbx import matrix
  H = - matrix.sqr((
    sum([ss*ss*cc/(B*bb+S*ss)**2 for cc, bb, ss in zip(c, b, s)]),
    sum([ss*bb*cc/(B*bb+S*ss)**2 for cc, bb, ss in zip(c, b, s)]),
    sum([ss*bb*cc/(B*bb+S*ss)**2 for cc, bb, ss in zip(c, b, s)]),
    sum([bb*bb*cc/(B*bb+S*ss)**2 for cc, bb, ss in zip(c, b, s)])))

  Hinv = H.inverse()

  delta = matrix.col((
    sum([ss*cc/(B*bb+S*ss) - ss for cc, bb, ss in zip(c, b, s)]),
    sum([bb*cc/(B*bb+S*ss) - bb for cc, bb, ss in zip(c, b, s)])))

  x = matrix.col((S, B))
  return x - Hinv * delta

def test(counts, background_shape, signal_shape):
  S = 1
  B = 1
  for i in range(20):
    S, B = iteration(counts, background_shape, signal_shape, B, S)
    print S, B
  # from dials.array_family import flex
  # Fs = sum(signal_shape)
  # Fb = sum(background_shape)
  # Rs = flex.double(flex.grid(100, 100))
  # for B in range(1, 101):
  #   for S in range(1, 101):
  #     Fb2 = sum([b*c/(B*b+S*s) for c, b, s in zip(counts, background_shape, signal_shape)])
  #     Fs2 = sum([s*c/(B*b+S*s) for c, b, s in zip(counts, background_shape, signal_shape)])
  #     R = (Fb2 - Fb)**2 + (Fs2 - Fs)**2
  #     Rs[B-1,S-1] = R
  # from matplotlib import pylab
  # pylab.imshow(flex.log(Rs).as_numpy_array())
  # pylab.show()




if __name__ == '__main__':

  from numpy.random import poisson
  from math import exp

  mean = [5*exp(-(x - 10.0)**2 / (2*3.0**2)) for x in range(20)]
  signal = mean
  # signal = [poisson(m, 1)[0] for m in mean]

  background = [4 for x in range(20)]
  # background = [poisson(50, 1)[0] for x in range(20)]

  total = [b + s for b, s in zip(background, signal)]

  print "Signal: ", sum(signal)
  print "Background: ", sum(background)
  print "R: ", sum(signal) / sum(background)

  C = total
  B = [1.0 / 20 for i in range(20)]
  S = [m / sum(mean) for m in mean]
  print test(C, B, S)

  # from matplotlib import pylab
  # pylab.plot(signal)
  # # pylab.plot(background)
  # pylab.plot(total)
  # pylab.show()
