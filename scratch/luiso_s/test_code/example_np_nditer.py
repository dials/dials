import numpy as np
import time

def tst_func(n_in):
  n_out = 2 * n_in
  return n_out



if(__name__ == "__main__"):

  width = 4
  height = 5

  a = np.arange(width * height).reshape(width, height)
  b = np.zeros( (width, height, 3), 'uint8')

  #b[:,:] = a[:,:]
  print ("a =")
  print (a)
  print ("b =")
  print (b)

  b[:,:,0] = a[:,:]
  b[:,:,1] = a[:,:]
  b[:,:,2] = a[:,:]

  print
  print ("after applying")
  print ("after b[:,:,0] = a[:,:]")
  print ("after b[:,:,1] = a[:,:]")
  print ("after b[:,:,2] = a[:,:]")

  print ("b =")
  print (b)

  time1 = time.time()
  for x in np.nditer(b, op_flags=['readwrite'], flags=['external_loop']):
    x[...] = tst_func(x)
  time2 = time.time()
  time_dif = time2 - time1

  print
  print ("<<<< after"), (time_dif), ("time >>>>")
  print ("a =")
  print (a)
  print ("b =")
  print (b)
