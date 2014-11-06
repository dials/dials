from scitbx.array_family import flex
import scitbx.lbfgs
import math

def residual(Oij, Gi, Ij):
  '''Oij, Gi, Ij are *dictionaries*'''
  f = 0.0
  for i, j in Oij:
    f += (Oij[(i, j)] - Gi[i] * Ij[j]) ** 2
  return f

def gradients_Gi(Oij, Gi, Ij):
  '''Oij, Gi, Ij are *dictionaries*'''
  gGi = { }
  for i in Gi:
    gGi[i] = 0.0
    for j in Ij:
      if (i, j) in Oij:
        gGi[i] += -2 * Ij[j] * Oij[(i, j)] + 2 * Gi[i] * (Ij[j] ** 2)

  return gGi

def gradients_Ij(Oij, Gi, Ij):
  '''Oij, Gi, Ij are *dictionaries*'''
  gIj = { }
  for j in Ij:
    gIj[j] = 0.0
    for i in Gi:
      if (i, j) in Oij:
        gIj[j] += -2 * Gi[i] * Oij[(i, j)] + 2 * (Gi[i] ** 2) * Ij[j]

  return gIj

class scaler:
  def __init__(self, Oij):
    # first scrape through observation dictionary Oij & work out how
    # many Gi & Ij we have

    self._i = []
    self._j = []

    for (i, j) in Oij:
      if not i in self._i:
        self._i.append(i)
      if not j in self._j:
        self._j.append(j)

    self._Oij = Oij
    self.x = flex.double([1 for i in self._i] + [1 for j in self._j])
    scitbx.lbfgs.run(target_evaluator=self)

    return

  def compute_functional_and_gradients(self):
    Gi, Ij = self.get_Gi_Ij()
    f = residual(self._Oij, Gi, Ij)
    gGi = gradients_Gi(self._Oij, Gi, Ij)
    gIj = gradients_Ij(self._Oij, Gi, Ij)
    g = flex.double()
    for _i in self._i:
      g.append(gGi[_i])
    for _j in self._j:
      g.append(gIj[_j])
    return f, g

  def get_Gi_Ij(self):
    # unpack self.x to Gi and Ij
    Gi = { }
    for i, _i in enumerate(self._i):
      Gi[_i] = self.x[i]
    Ij = { }
    for j, _j in enumerate(self._j):
      Ij[_j] = self.x[j + len(self._i)]
    return Gi, Ij

def test(nG, nI):
  from random import random

  Gi_t = [random() for i in range(nG)]
  Ij_t = [random() for j in range(nI)]

  Oij = { }
  for i in range(nG):
    for j in range(nI):
      Oij[(i, j)] = Gi_t[i] * Ij_t[j]

  s = scaler(Oij)

  Gi, Ij = s.get_Gi_Ij()

  scale_G = sum([Gi[i] for i in Gi]) / sum(Gi_t)
  scale_I = sum([Ij[j] for j in Ij]) / sum(Ij_t)

  rms_G = sum([(Gi[i] / scale_G - Gi_t[i]) ** 2 for i in sorted(Gi)])
  rms_I = sum([(Ij[j] / scale_I - Ij_t[j]) ** 2 for j in sorted(Ij)])
  return math.sqrt(rms_G / len(Gi)), math.sqrt(rms_I / len(Ij)), \
      scale_G * scale_I

def meansd(values):
  mean = sum(values) / len(values)
  var = sum([(v - mean) ** 2 for v in values]) / len(values)
  return mean, math.sqrt(var)

if __name__ == '__main__':
  print ' nG    nI   rmsG  rmsI'
  for pnG in range(3, 9):
    nG = 2 ** pnG
    for pnI in range(3, 9):
      nI = 2 ** pnI
      rmsG, rmsI, s = test(nG, nI)
      print '%5d %5d %.3f %.3f %.3f' % (nG, nI, rmsG, rmsI, s)
