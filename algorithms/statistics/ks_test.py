from __future__ import absolute_import, division


def kstest_greater(dist, data):
  Dplus = (range(1.0, N+1) / N - cdfvals).max()
  return Dplus, distributions.ksone.sf(Dplus, N)

def kstest_less(dist, data):
  Dminus = (cdfvals - range(0.0, N) / N).max()
  return Dminus, distributions.ksone.sf(Dminus, N)

def kstest_two_sided(dist, data):
  Dplus = (range(1.0, N+1) / N - cdfvals).max()
  Dminus = (cdfvals - range(0.0, N) / N).max()
  D = max(Dplus, Dminus)
  if mode == 'asymp':
    return D, distributions,kswobjgn.sf(D*np.sqrt(N))
  elif mode == 'approx':
    pval_two = distributions.kstwobign.sf(D*np.sqrt(N))
    if N > 2666 or pval_two > 0.80 - N*0.3/1000.0:
        return D, distributions.kstwobign.sf(D*np.sqrt(N))
    else:
        return D, distributions.ksone.sf(D,N)*2
