from scitbx import matrix






# def f(x):
#   from math import log
#   return -(log(x) - x)

# def df(x):
#   return -(1.0 / x - 1)


# x0 = 10
# B0 = 1

# for i in range(1000):
#   p = -df(x0) * B0
#   alpha = 1.0
#   x1 = x0 + alpha * p
#   s = alpha*p
#   y = df(x1) - df(x0)
#   B1 = B0 + (s*y + y*B0*y)*(s*s) / (s*y)**2 - (B0*y*s + s*y*B0) / (s*y)
#   B0 = B1
#   x0 = x1
#   print B1, x0
def func(c, b, s, B, S):
  from math import log, factorial
  sumc = sum([log(factorial(cc)) for cc in c])
  sumb = sum(b)
  sums = sum(s)
  suml = 0
  for i in range(len(c)):
    if B*b[i] + S*s[i] <= 0:
      return 0, True
    suml += c[i] * log(B*b[i] + S*s[i])
  return -(suml - sumc - B*sumb - S*sums), False

def db(c, b, s, B, S):
  sumb = sum(b)
  sumc = 0
  for i in range(len(c)):
    assert(B*b[i] + S*s[i] > 0)
    sumc += c[i]*b[i] / (B*b[i] + S*s[i])
  return sumc - sumb

def ds(c, b, s, B, S):
  sums = sum(s)
  sumc = 0
  for i in range(len(c)):
    assert(B*b[i] + S*s[i] > 0)
    sumc += c[i]*s[i] / (B*b[i] + S*s[i])
  return sumc - sums

def df(c, b, s, x0):
  return matrix.col((
    -db(c, b, s, x0[0], x0[1]),
    -ds(c, b, s, x0[0], x0[1])))

def show_image(c,b,s, BB=None, SS=None):

  import numpy.ma
  from dials.array_family import flex
  N = 100
  im = flex.double(flex.grid(N, N))
  mask = flex.bool(flex.grid(N, N))
  for j in range(N):
    for i in range(N):
      B = -1.0 + j * 10.0 / N
      S = -1.0 + i * 10.0 / N
      im[j,i], mask[j,i] = func(c,b,s,B,S)
      im[j,i] = -im[j,i]

  masked_im = numpy.ma.array(
    # im.as_numpy_array(),
    flex.exp(im).as_numpy_array(),
    mask = mask.as_numpy_array())
  mask2 = flex.bool(flex.grid(N, N))
  indices = []
  for i in range(len(mask)):
    if mask[i] == False:
      indices.append(i)
  indices = flex.size_t(indices)
  ind = flex.max_index(im.select(indices))
  ind = indices[ind]
  maxy = -1.0 + (ind % N) * 10.0 / N
  maxx = -1.0 + (ind // N) * 10.0 / N
  from matplotlib import pylab
  pylab.imshow(masked_im, origin='bottom', extent=[-1.0, 9.0, -1.0, 9.0])
  if YY is not None and XX is not None:
    pylab.plot(YY, XX)
  pylab.scatter([maxy], [maxx])
  pylab.show()

  # pylab.plot(masked_im[11,:])
  # pylab.show()
  # exit(0)


def linesearch(c, b, s, x0, p):
  c1 = 1e-4
  c2 = 0.5

  def psi(alpha):
    return func(c, b, s, x0[0]+alpha*p[0], x0[1]+alpha*p[1])[0]

  def dpsi(alpha):
    sums = sum(s)
    sumb = sum(b)
    sumc = 0
    for i in range(len(c)):
      den = (x0[0]*b[i]+x0[1]*s[i])+alpha*(p[0]*b[i]+p[1]*s[i])
      assert(den > 0)
      sumc += c[i]*(p[0]*b[i]+p[1]*s[i]) / den
    return -(sumc - p[0]*sumb - p[1]*sums)

  def zoom(a0, a1, d0, f0, fa0):
    j = 0
    while (j < 10):
      aj = (a1 + a0) / 2.0
      faj = psi(aj)
      if faj > f0 + c1*aj*d0 or faj >= fa0:
        a1 = aj
      else:
        daj = dpsi(aj)
        if abs(daj) <= -c2 * d0:
          return aj
        if daj*(a1 - a0) >= 0:
          a1 = a0
        a0 = aj
      j += 1

  a0 = 0
  am = 1
  a1 = am / 2.0
  f0 = psi(a0)
  d0 = dpsi(a0)
  fa0 = f0
  da0 = d0
  i = 1
  alpha = None
  while (True):
    fa1 = psi(a1)
    da1 = dpsi(a1)
    if (fa1 > f0 + c1*a1*d0) or (fa1 >= fa0 and i > 1):
      alpha = zoom(a0, a1, d0, f0, fa0)
      break
    if abs(da1) <= -c2*d0:
      alpha = a1
      break
    if da1 >= 0:
      alpha = zoom(a1, a0, d0, f0, fa1)
      break
    a0 = a1
    fa0 = fa1
    da0 = da1
    a1 = (a1 + am) / 2.0
    i += 1
    # print a1, f0, fa0, fa1
  # print alpha
  return alpha
  # exit(0)


def linesearch2(c, b, s, x0, p):
  def psi(alpha):
    return func(c, b, s, x0[0]+alpha*p[0], x0[1]+alpha*p[1])[0]

  def dpsi(alpha):
    sums = sum(s)
    sumb = sum(b)
    sumc = 0
    for i in range(len(c)):
      den = (x0[0]*b[i]+x0[1]*s[i])+alpha*(p[0]*b[i]+p[1]*s[i])
      assert(den > 0)
      sumc += c[i]*(p[0]*b[i]+p[1]*s[i]) / den
    return -(sumc - p[0]*sumb - p[1]*sums)

  def d2psi(alpha):
    sumc = 0
    for i in range(len(c)):
      den = (x0[0]*b[i]+x0[1]*s[i])+alpha*(p[0]*b[i]+p[1]*s[i])
      assert(den > 0)
      sumc += c[i]*(p[0]*b[i]+p[1]*s[i])**2 / den**2
    return (sumc)


  f = psi(0)
  fp = dpsi(0)
  fpp = d2psi(0)
  K = sum([b[i]*p[0]+s[i]*p[1] for i in range(len(c))])
  A = (K - fp) / fpp - 0
  C = (K - fp)**2 / fpp
  print A, C, K
  print C / K - A

  F1 = []
  F2 = []
  from math import log
  for a in range(100):
    aa = 0.1 + -A + a / 10.0
    F2.append(-C*log(A+aa) + K*aa)
    F1.append(psi(aa))
  from matplotlib import pylab
  pylab.plot(F1)
  pylab.plot(F2)
  pylab.show()

  return C / K - A

def sign(a):
  if a < 0:
    return -1
  else:
    return 1


def linesearch3(c, b, s, x0, p):

  def psi(alpha):
    return func(c, b, s, x0[0]+alpha*p[0], x0[1]+alpha*p[1])[0]

  def dpsi(alpha):
    sums = sum(s)
    sumb = sum(b)
    sumc = 0
    for i in range(len(c)):
      den = (x0[0]*b[i]+x0[1]*s[i])+alpha*(p[0]*b[i]+p[1]*s[i])
      assert(den > 0)
      sumc += c[i]*(p[0]*b[i]+p[1]*s[i]) / den
    return -(sumc - p[0]*sumb - p[1]*sums)

  alpha0 = 0
  alpha1 = 1
  d0 = dpsi(alpha0)
  f0 = psi(alpha0)
  f00 = f0
  d00 = d0

  c1 = 1e-4
  c2 = 0.9

  while (True):
    d1 = dpsi(alpha1)
    f1 = psi(alpha1)
    if sign(d1) != sign(d0):
      break
    alpha0 = alpha1
    d0 = d1
    f0 = f1
    alpha1 *= 1.6128
    print alpha0, alpha1

  assert(sign(d0) != sign(d1))
  i = 0
  while (i < 10):
    alpha = 0.5 * (alpha0 + alpha1)
    d = dpsi(alpha)
    f = psi(alpha)
    if sign(d) == sign(d0):
      alpha0 = alpha
      d0 = d
      f0 = f
    else:
      alpha1 = alpha
      d1 = d
      f1 = f
    bound = f00 + c1*alpha*d00
    print alpha, alpha0, alpha1, d, bound, d00
    if (f <= f00 + c1*alpha*d00 and abs(d) <= c2*abs(d00)):
      break

# print alpha0, alpha1, f, f00, f00 +c1*alpha*d00

    i += 1
  return alpha


def linesearch4(c, b, s, x0, p):


  c1 = 1e-4
  c2 = 0.9

  def psi(alpha):
    return func(c, b, s, x0[0]+alpha*p[0], x0[1]+alpha*p[1])[0]

  def dpsi(alpha):
    sums = sum(s)
    sumb = sum(b)
    sumc = 0
    for i in range(len(c)):
      den = (x0[0]*b[i]+x0[1]*s[i])+alpha*(p[0]*b[i]+p[1]*s[i])
      assert(den > 0)
      sumc += c[i]*(p[0]*b[i]+p[1]*s[i]) / den
    return -(sumc - p[0]*sumb - p[1]*sums)

  def zoom(alpha0, alpha1, f00, d00, f0):
    while (True):
      alpha = 0.5*(alpha0 + alpha1)
      f1 = psi(alpha)
      d1 = dpsi(alpha)
      if f1 > f00 + c1*alpha*d00 or (f1 >= f0 and i > 0):
        alpha1 = alpha
      else:
        if (abs(d1) <= -c2*d00):
          return alpha
        if d1 * (alpha1 - alpha0) >= 0:
          alpha1 = alpha0
        alpha0 = alpha



  max_iter = 10
  alpha0 = 0
  alpha1 = 1
  d00 = dpsi(alpha0)
  f00 = psi(alpha0)
  f0 = f00
  d0 = d00
  print d00
  assert(d00 < 0)
  for i in range(max_iter):
    f1 = psi(alpha1)
    if f1 > f00 + c1*alpha1*d00 or (f1 >= f0 and i > 0):
      alpha1 = zoom(alpha0, alpha1, f00, d00, f0)
      break
    d1 = dpsi(alpha1)
    if abs(d1) <= -c2*d00:
      break
    if d1 >= 0:
      alpha1 = zoom(alpha1, alpha0, f00, d00, f1)
      break
    f0 = f1
    d0 = d1
    alpha0 = alpha1
    alpha1 *= 1.618
  return alpha1



c = [0, 1, 0, 0, 0, 0, 3, 1, 3, 3, 6, 6, 4, 1, 4, 0, 2, 0, 1, 1]
# c = [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
b = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
s = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
b = [bb / sum(b) for bb in b]
s = [ss / sum(s) for ss in s]

x0 = matrix.col((0.2, 0.2))
B0 = matrix.sqr((1, 0, 0, 1))

K1 = 0
K2 = 0.45

XX = [x0[0]]
YY = [x0[1]]

from dials.algorithms.integration.profile import MLPoisson2Stepper
from dials.array_family import flex

stepper = MLPoisson2Stepper(
  flex.double(c),
  flex.double(b),
  flex.double(s), x0)

# try:
#   for i in range(10):
#     stepper.step()
#     print stepper.X()
#     XX.append(stepper.X()[0])
#     YY.append(stepper.X()[1])
# except Exception, e:
#   print e
#   pass
# show_image(c, b, s)
try:
  for i in range(5):
    DFX0 = df(c,b,s,x0)
    p = -B0*df(c, b, s, x0)
    print "P: ", tuple(p)
    d = x0 + p

    #alpha = 0.5#1.0

    den1 = (p[0]+p[1]*K1)
    den2 = (K2*p[0]+p[1])
    alpha_max1 = -(x0[0] + x0[1]*K1) / den1
    alpha_max2 = -(K2*x0[0] + x0[1]) / den2
    print "Alpha Max: ", tuple(p), alpha_max1, alpha_max2

    alpha = linesearch4(c, b, s, x0, p)
    print alpha
    # if alpha > 0.45:
    #   print i, tuple(B0), tuple(p), den1, den2, alpha1, alpha2, alpha

    # sumc = sum(c)
    # sumb = sum(b)
    # sums = sum(s)
    # sump = d[0]*sumb+d[1]*sums
    # alpha = sumc / sump
    # print alpha

    # alpha = 1.0
    S = alpha * p
    x1 = x0 + S
    Y = df(c, b, s, x1) - df(c, b, s, x0)
    ST = S.transpose()
    YT = Y.transpose()
    YTS = (YT * S)[0]
    I = matrix.sqr((1, 0, 0, 1))
    B11 = (I - S*YT / YTS)
    B12 = (I - Y*ST / YTS)
    B13 = S*ST / YTS
    B1 = B11*B0*B12+B13

    # STY = (ST * Y)[0]
    # YTB0Y = (YT * B0 * Y)[0]
    # SST = S * ST
    # B11 = (STY + YTB0Y)*(SST)/(STY)**2
    # B12 = (B0*Y*ST+S*YT*B0)/(STY)
    # B1 = B0 + B11 - B12
    x0 = x1
    B0 = B1
    print S*YT
    print Y*ST
    print S*ST
    print B0

    # print tuple(x0), alpha, alpha1, alpha2, tuple(p), tuple(DFX0)
    XX.append(x0[0])
    YY.append(x0[1])
    print tuple(x0)
except Exception, e:
  # raise
  pass


show_image(c, b, s, XX, YY)

# from matplotlib import pylab
# pylab.plot(YY, XX)
# pylab.show()
