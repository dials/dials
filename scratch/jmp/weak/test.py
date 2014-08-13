


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
# def iteration(c, b, s, B, S):
#   from scitbx import matrix

#   sum_ss = 0
#   sum_sb = 0
#   sum_bb = 0
#   sum_s = 0
#   sum_b = 0
#   TINY = 1e-5
#   for i in range(len(c)):
#     cc = c[i]
#     ss = s[i]
#     bb = b[i]
#     den = max(TINY, B*bb + S*ss)
#     sum_s += ss*cc/den - ss
#     sum_b += bb*cc/den - bb
#     sum_ss += ss*ss*cc/den**2
#     sum_bb += bb*bb*cc/den**2
#     sum_sb += bb*ss*cc/den**2

#   H = - matrix.sqr((sum_ss, sum_sb,
#                     sum_sb, sum_bb))
#   try:
#     Hinv = H.inverse()
#   except Exception, e:
#     import traceback
#     traceback.print_exc()
#     print H
#     print B, S
#     raise e

#   delta = matrix.col((sum_s, sum_b))

#   x = matrix.col((S, B))
#   return x - Hinv * delta
# def iteration(c, b, s, B, S):
#   from scitbx import matrix

#   TINY = 1e-7
#   sum_ss = 0
#   sum_sb = 0
#   sum_bb = 0
#   sum_s = 0
#   sum_b = 0
#   for i in range(len(c)):
#     cc = c[i]
#     ss = s[i]
#     bb = b[i]
#     den = B*bb + S*ss
#     assert(den > TINY)
#     sum_s += ss*cc/den - ss
#     sum_b += bb*cc/den - bb
#     sum_ss += ss*ss*cc/den**2
#     sum_bb += bb*bb*cc/den**2
#     sum_sb += bb*ss*cc/den**2

#   H = - matrix.sqr((sum_ss, sum_sb,
#                     sum_sb, sum_bb))
#   print c
#   print b
#   print s
#   print H
#   Hinv = H.inverse()

#   delta = matrix.col((sum_s, sum_b))
#   dd = delta.transpose() * Hinv * delta
#   den = 1.0 - dd[0]

#   d = Hinv * delta / den

#   x = matrix.col((S, B))
#   return x - d



def ml_poisson2_step(c, b, s, B, S):
  step = 1.0
  sum_b = sum(b)
  sum_s = sum(s)
  sum_sc = sum([cc*ss/(B*bb+S*ss) for cc, bb, ss in zip(c, b, s)])
  sum_bc = sum([cc*bb/(B*bb+S*ss) for cc, bb, ss in zip(c, b, s)])
  dldb = sum_bc - sum_b
  dlds = sum_sc - sum_s
  DB = step * dldb
  DS = step * dlds
  B1 = B + DB
  S1 = S + DS
  print B1, S1
  if B1 < 0:
    B1 = B / 2
  if S1 < 0:
    S1 = S / 2
  print B1, S1
  return B1, S1





def iteration(c, b, s, B, S):
  from scitbx import matrix


  K1 = 0
  K2 = 0.05 / 0.1111111111

  eps = 1.0
  sum_b = sum(b)
  sum_s = sum(s)
  sum_sc = sum([cc*ss/(B*bb+S*ss) for cc, bb, ss in zip(c, b, s)])
  sum_bc = sum([cc*bb/(B*bb+S*ss) for cc, bb, ss in zip(c, b, s)])
  lb = 1.0 / (B + K1*S) + K2 / (S + K2*B)
  ls = K1 / (B + K1*S) + 1.0 / (S + K2*B)
  dldb = sum_bc - sum_b + lb
  dlds = sum_sc - sum_s + ls
  DB = eps * dldb
  DS = eps * dlds
  B1 = B + DB
  S1 = S + DS
  print B1, S1
  # print ":", dldb, dlds
  # exit(0)
  # if B1 < 0:
  #   B1 = B / 2
  # if S1 < 0:
  #   S1 = S / 2
  # print B1, S1
  return B1, S1




def evaluate(c, b, s, B, S):
  DB = sum([bb*cc/(B*bb+S*ss)-bb for cc, bb, ss in zip(c, b, s)])
  DS = sum([ss*cc/(B*bb+S*ss)-ss for cc, bb, ss in zip(c, b, s)])
  return DB, DS

def value(counts, background_shape, signal_shape):
  from math import log, exp
  background_shape = [float(b) / sum(background_shape) for b in background_shape]
  signal_shape = [float(s) / sum(signal_shape) for s in signal_shape]
  S = 1
  B = 1
  EPS = 1e-10
  print "Start"
  SSS = [S]
  BBB = [B]
  n = 0
  while(True):
    B1, S1 = iteration(counts, background_shape, signal_shape, B, S)
    # print S1, B1, S, B
    # print " ", S1, B1
    if n >= 100:#(S1-S)**2 + (B1-B)**2 < EPS:
      break
    S = S1
    B = B1
    SSS.append(S)
    BBB.append(B)
    n += 1

  from matplotlib import pylab
  pylab.plot(SSS, BBB)
  pylab.show()
  print counts
  print background_shape
  print signal_shape
  DB, DS = evaluate(counts, background_shape, signal_shape, 0.5236486486486487, 2.920795795795796)
  print "Derivatives: ", DB, DS
  DB, DS = evaluate(counts, background_shape, signal_shape, B1, S1)
  print "Derivatives: ", DB, DS

  return S1, B1

def sigmas(c, b, s, B, S):
  from scitbx import matrix
  invsb = sum([cc*bb*bb / (B*bb + S*ss)**2 for cc, bb, ss in zip(c, b, s)])
  invss = sum([cc*ss*ss / (B*bb + S*ss)**2 for cc, bb, ss in zip(c, b, s)])
  return 1.0 / invsb, 1.0 / invss

def test(counts, background_shape, signal_shape):
  S = 1
  B = 1
  for i in range(20):
    S, B = iteration(counts, background_shape, signal_shape, B, S)
    print S, B

  sigma_b, sigma_s = sigmas(counts, background_shape, signal_shape, B, S)
  print sigma_b, sigma_s
  from dials.array_family import flex
  Fs = sum(signal_shape)
  Fb = sum(background_shape)
  Rs = flex.double(flex.grid(100, 100))
  for B in range(1, 101):
    for S in range(1, 101):
      Fb2 = sum([b*c/(B*b+S*s) for c, b, s in zip(counts, background_shape, signal_shape)])
      Fs2 = sum([s*c/(B*b+S*s) for c, b, s in zip(counts, background_shape, signal_shape)])
      R = (Fb2 - Fb)**2 + (Fs2 - Fs)**2
      Rs[B-1,S-1] = Fs2 - Fs
  from matplotlib import pylab
  pylab.imshow(Rs.as_numpy_array())
  pylab.show()

def plot_prob_for_zero(c, b, s):
  from math import log, exp, factorial
  from dials.array_family import flex
  L = flex.double(flex.grid(100, 100))
  MASK = flex.bool(flex.grid(100, 100))
  c = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
  b = [bb / sum(b) for bb in b]
  s = [ss / sum(s) for ss in s]
  for BB in range(0, 100):
    for SS in range(0, 100):
      B = 0 + BB / 10000.0
      S = 0 + SS / 40.0
      LL = 0
      for i in range(len(b)):
        if B*b[i] + S*s[i] <= 0:
          MASK[BB, SS] = True
          LL = -999999
          break
        else:
          LL += c[i]*log(B*b[i]+S*s[i]) - log(factorial(c[i])) - B*b[i] - S*s[i]

      L[BB, SS] = LL
  index = flex.max_index(L)
  i = index % 100
  j = index // 100
  B = 0 + j / 10000.0
  S = 0 + i / 40.0
  print flex.max(L), B, S
  from matplotlib import pylab
  import numpy
  im = numpy.ma.masked_array(flex.exp(L).as_numpy_array(), mask=MASK.as_numpy_array())
  pylab.imshow(im)
  pylab.show()
  exit(0)


def plot_valid(b, s):
  from dials.array_family import flex

  b = [0.1, 0.2, 0.3, 0.4, 0.5]
  s = [0.1, 0.3, 0.5, 0.3, 0.1]


  v1 = flex.bool(flex.grid(100, 100))
  v2 = flex.bool(flex.grid(100, 100))
  v3 = flex.bool(flex.grid(100, 100))
  r = [float(ss) / float(bb) for ss, bb in zip(s, b)]

  for BB in range(0, 100):
    for SS in range(0, 100):
      B = -5.0 + BB / 10.0
      S = -5.0 + SS / 10.0
      V1 = True
      V2 = True
      V3 = True
      for i in range(len(b)):
        if B*b[i]+S*s[i] <= 0:
          V1 = False
          break
      for i in range(len(b)):
        if B*b[i] <= -S*s[i]:
          V2 = False
          break

      v1[BB,SS] = V1
      v2[BB,SS] = V2

  from matplotlib import pylab
  pylab.imshow(v1.as_numpy_array())
  pylab.show()
  exit(0)


def paper_test(B, S):

  from numpy.random import poisson
  from math import exp

  background_shape = [1 for i in range(20)]
  signal_shape = [1 if i >= 6 and i < 15 else 0 for i in range(20)]

  background = [poisson(bb * B,1)[0] for bb in background_shape]
  signal = [poisson(ss * S, 1)[0] for ss in signal_shape]
  # background = [bb * B for bb in background_shape]
  # signal = [ss * S for ss in signal_shape]
  total = [b + s for b, s in zip(background, signal)]

  # from matplotlib import pylab
  # pylab.plot(total)
  # pylab.plot(signal)
  # pylab.plot(background)
  # pylab.show()

  total = [0, 1, 0, 0, 0, 0, 3, 1, 3, 3, 6, 6, 4, 1, 4, 0, 2, 0, 1, 1]
  total = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]

  # plot_prob_for_zero(total, background_shape, signal_shape)
  # signal_shape = [exp(-(x - 10.0)**2 / (2*3.0**2)) for x in range(20)]
  # signal_shape = [ss / sum(signal_shape) for ss in signal_shape]
  # print signal_shape
  #plot_valid(background_shape, signal_shape)

  B = 155.0 / 296.0
  from dials.array_family import flex
  from math import log, factorial
  V = flex.double(flex.grid(100, 100))
  L = flex.double(flex.grid(100, 100))
  DB = flex.double(flex.grid(100, 100))
  DS = flex.double(flex.grid(100, 100))
  P = flex.double(flex.grid(100, 100))
  Fb = sum(background_shape)
  Fs = sum(signal_shape)
  SV = []
  MASK = flex.bool(flex.grid(100, 100), False)
  for BB in range(100):
    for SS in range(100):
      B = -5.0 + (BB) / 10.0
      S = -5.0 + (SS) / 10.0
      # SV.append(S)
      VV = 0
      LL = 0
      DDB = 0
      DDS = 0
      for i in range(20):
        s = signal_shape[i]
        b = background_shape[i]
        c = total[i]
        if B*b + S*s <= 0:
          # MASK[BB, SS] = True
          LL = 0
          if b == 0:
            DDB += 0
          else:
            DDB += 1e7
          if s == 0:
            DDS += 0
          else:
            DDS += 1e7
          # break
        else:
        # VV += (b + s)*c / (B*b + S*s)
        # LL += c*log(B*b+S*s) - log(factorial(c)) - B*b - S*s
          DDB += c*b/(B*b+S*s) - b
          DDS += c*s/(B*b+S*s) - s
      VV -= (Fb + Fs)
      # print B, S, VV
      # V[BB, SS] = abs(VV)
      L[BB, SS] = LL
      DB[BB,SS] = DDB
      DS[BB,SS] = DDS

  max_ind = flex.max_index(L)
  j = max_ind // 100
  i = max_ind % 100
  print "Approx: ", (j+1) / 20.0, (i+1) / 20.0
  print "Min/Max DB: ", flex.min(DB), flex.max(DB)
  print "Min/Max DS: ", flex.min(DS), flex.max(DS)
  from matplotlib import pylab
  # pylab.imshow(flex.log(V).as_numpy_array(), extent=[0.05, 5.05, 5.05, 0.05])
  # pylab.plot(SV, V)
  # pylab.plot(SV, [0] * 100)
  # pylab.show()
  im = flex.exp(L).as_numpy_array()
  import numpy
  # im = numpy.ma.masked_array(im, mask=MASK.as_numpy_array())
  # pylab.imshow(im)#, extent=[-5.0, 5.0, 5.0, -5.0], origin='lower')
  pylab.imshow(DB.as_numpy_array(), vmin=-100, vmax=100)#, extent=[-5.0, 5.0, 5.0, -5.0], origin='lower')
  pylab.contour(DB.as_numpy_array(), levels=[0], colors=['red'])
  pylab.contour(DS.as_numpy_array(), levels=[0], colors=['black'])
  pylab.show()
  # im = numpy.ma.masked_array(DB.as_numpy_array(), mask=MASK.as_numpy_array())
  # pylab.imshow(im, extent=[-5.0, 5.0, 5.0, -5.0], vmin=-20, vmax=100)
  # pylab.show()
  # im = numpy.ma.masked_array(DS.as_numpy_array(), mask=MASK.as_numpy_array())
  # pylab.imshow(im, extent=[-5.0, 5.0, 5.0, -5.0], vmin=-20, vmax=100)
  # pylab.show()
  # exit(0)

  S1, B1 = value(total, background_shape, signal_shape)
  exit(0)
  try:
    S1, B1 = value(total, background_shape, signal_shape)
    print "Result:"
    print S1, B1
    exit(0)
  except Exception, e:
    raise e
    import sys
    import traceback
    traceback.print_exc()
    from dials.array_family import flex
    Fs = sum(signal_shape)
    Fb = sum(background_shape)
    Rs = flex.double(flex.grid(100, 100))
    print "-----"
    print B, S
    # from matplotlib import pylab
    # pylab.plot(total)
    # pylab.show()
    print background_shape
    print signal_shape
    print total
    from math import exp, factorial, log
    minx = -1
    miny = -1
    minr = 9999
    for BB in range(0, 100):
      for SS in range(0, 100):
        B = -10 + (BB) / 5.0
        S = -10 + (SS) / 5.0
        L = 0
        Fb2 = 0
        Fs2 = 0
        for i in range(len(total)):
          c = total[i]
          b = background_shape[i]
          s = signal_shape[i]
          # P = exp(-(B*b + S*s)) * (B*b+S*s)**c / factorial(c)
          # print P
          # if P > 0:
          #   L += log(P)
          den = B*b + S*s
          num1 = b*c
          num2 = s*c
          if den != 0:
            Fb2 += num1 / den
            Fs2 += num2 / den
        R = (Fb2 - Fb)**2 + (Fs2 - Fs)**2
        if R > 1000:
          R = 0
        # Rs[BB,SS] = L#R#Fs2 - Fs
        Rs[BB,SS] = R#Fs2 - Fs
    from matplotlib import pylab
    pylab.imshow(flex.log(Rs).as_numpy_array(), extent=[-5,5,5,-5])
    pylab.show()
    exit(0)
  exit(0)

  # print S, B, sum(signal), sum(background), S1, B1
  return S, B, sum(signal), sum(background), S1, B1


def single(B, S):
  from numpy.random import poisson
  from math import exp

  background_shape = [1.0 / 20.0 for i in range(20)]
  signal_shape = [exp(-(x - 10.0)**2 / (2*3.0**2)) for x in range(20)]
  signal_shape = [ss / sum(signal_shape) for ss in signal_shape]

  background = [poisson(bb * B,1)[0] for bb in background_shape]
  signal = [poisson(ss * S, 1)[0] for ss in signal_shape]
  # background = [bb * B for bb in background_shape]
  # signal = [ss * S for ss in signal_shape]
  total = [b + s for b, s in zip(background, signal)]

  # from matplotlib import pylab
  # pylab.plot(total)
  # pylab.plot(signal)
  # pylab.plot(background)
  # pylab.show()

  S1, B1 = value(total, background_shape, signal_shape)
  # print S, B, sum(signal), sum(background), S1, B1
  return S, B, sum(signal), sum(background), S1, B1


if __name__ == '__main__':
  from random import uniform

  Sexp = []
  Bexp = []
  Scnt = []
  Bcnt = []
  Sest = []
  Best = []
  for i in range(100):
    B = 1#*20# * 200# * uniform(1, 100)
    S = 2#uniform(1000, 2000)
    # S, B, Stot, Btot, Scal, Bcal = single(B, S)
    S, B, Stot, Btot, Scal, Bcal = paper_test(B, S)
    Sexp.append(S)
    Bexp.append(B)
    Scnt.append(Stot)
    Bcnt.append(Btot)
    Sest.append(Scal)
    Best.append(Bcal)
    print S, B, Stot, Btot, Scal, Bcal
  from matplotlib import pylab
  pylab.scatter(Best, Sest)

#  pylab.scatter(Sexp, [s1 - s2 for s1, s2 in zip(Sexp, Sest)])
  pylab.show()

  # single(20 * 1, 50)
  # single(20 * 5, 50)
  # single(20 * 10, 50)
  # single(20 * 15, 50)
  # single(20 * 20, 50)
  # single(20 * 25, 50)
  # single(20 * 50, 50)
  # single(20 * 100, 50)
  # single(20, 1000)
  # from numpy.random import poisson
  # from math import exp

  # mean = [5*exp(-(x - 10.0)**2 / (2*3.0**2)) for x in range(20)]
  # signal = mean
  # # signal = [poisson(m, 1)[0] for m in mean]

  # background = [4 for x in range(20)]
  # # background = [poisson(50, 1)[0] for x in range(20)]

  # total = [b + s for b, s in zip(background, signal)]

  # print "Signal: ", sum(signal)
  # print "Background: ", sum(background)
  # print "R: ", sum(signal) / sum(background)

  # C = total
  # B = [1.0 / 20 for i in range(20)]
  # S = [m / sum(mean) for m in mean]
  # print test(C, B, S)

  # # from matplotlib import pylab
  # # pylab.plot(signal)
  # # # pylab.plot(background)
  # # pylab.plot(total)
  # # pylab.show()
