


if __name__ == '__main__':
  import numpy as np
  from numpy import ma

  data = np.array([
    [1, 2, 3, 4],
    [1, 2, 3, 4],
    [1, 2, 3, 4],
    [1, 2, 3, 4]
  ])

  mask = np.array([
    [False, True, True, False],
    [False, True, True, False],
    [False, True, True, False],
    [False, True, True, False]
  ])

  masked_data = ma.MaskedArray(data=data, mask=mask)

  fft_data =  np.fft.fft2(masked_data)
