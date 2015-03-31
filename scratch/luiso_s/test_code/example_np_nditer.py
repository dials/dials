import numpy as np
import time

def red_func(itnts, i_min, rim1, rim2, i_max, blk_siz):

  if( itnts > i_min and itnts < rim1 ):
    # first block (black to red)
    uni = float(itnts - i_min) / blk_siz
    red_out = int(uni * 255.0)

  elif( itnts >= rim1 and itnts <= i_max ):
    # second block (red to yellow)
    # or
    # third block (yellow to white)
    red_out = 255

  else:
    # defaulting to black
    red_out = 0

  return red_out


def green_func(itnts, i_min, rim1, rim2, i_max, blk_siz):

  if( itnts > i_min and itnts < rim1 ):
    # first block (black to red)
    green_out = 0

  elif( itnts >= rim1 and itnts <= rim2 ):
    # second block (red to yellow)
    uni = float(itnts - rim1) / blk_siz
    green_out = int(uni * 255.0)

  elif( itnts > rim2 and itnts <= i_max ):
    # third block (yellow to white)
    green_out = 255
  else:
    # defaulting to black
    green_out = 0

  return green_out


def blue_func(itnts, i_min, rim1, rim2, i_max, blk_siz):

  if( itnts > i_min and itnts < rim1 ):
    # first block (black to red)
    blue_out = 0

  elif( itnts >= rim1 and itnts <= rim2 ):
    # second block (red to yellow)
    blue_out = 0

  elif( itnts > rim2 and itnts <= i_max ):
    # third block (yellow to white)
    uni = float(itnts - rim2) / blk_siz
    blue_out = int(uni * 255.0)

  else:
    # defaulting to black
    blue_out = 0

  return blue_out


if(__name__ == "__main__"):

  width = 4
  height = 5

  a = np.arange(width * height).reshape(width, height)
  b = np.zeros( (width, height, 3), 'uint8')

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

  ints_min = np.min(a)
  ints_max = np.max(a)
  print ("ints_min ="), (ints_min)
  print ("ints_max ="), (ints_max)

  block_size = float(ints_max - ints_min) / 3.0
  rim_n1 = ints_min + block_size
  rim_n2 = ints_max - block_size

  fast_way = '''
  for x in np.nditer(b[:,:,0], op_flags=['readwrite'], flags=['external_loop']):
    x[...] = red_func(x, ints_min, rim_n1, rim_n2, ints_max, block_size)
  for x in np.nditer(b[:,:,1], op_flags=['readwrite'], flags=['external_loop']):
    x[...] = green_func(x, ints_min, rim_n1, rim_n2, ints_max, block_size)
  for x in np.nditer(b[:,:,2], op_flags=['readwrite'], flags=['external_loop']):
    x[...] = blue_func(x, ints_min, rim_n1, rim_n2, ints_max, block_size)
  '''


  for x in np.nditer(b[:,:,0], op_flags=['readwrite']):
    x[...] = red_func(x, ints_min, rim_n1, rim_n2, ints_max, block_size)
  for x in np.nditer(b[:,:,1], op_flags=['readwrite']):
    x[...] = green_func(x, ints_min, rim_n1, rim_n2, ints_max, block_size)
  for x in np.nditer(b[:,:,2], op_flags=['readwrite']):
    x[...] = blue_func(x, ints_min, rim_n1, rim_n2, ints_max, block_size)




  time2 = time.time()
  time_dif = time2 - time1

  print
  print ("<<<< after"), (time_dif), ("time >>>>")
  print ("a =")
  print (a)
  print ("b =")
  print (b)
