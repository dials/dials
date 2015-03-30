# one Liner for loop

import numpy as np

a = np.zeros((5), 'int')

if(__name__ == "__main__"):
  #a[p] = [  p  for p in np.ndarray(shape=(5), dtype=int) ]

  for p in np.ndarray(shape=(5), dtype = int):
    a[p] = p
  print a
