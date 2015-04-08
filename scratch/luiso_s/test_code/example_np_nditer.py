import numpy as np
import time

def red_func(itnts, i_min, rim1, rim2, i_max, blk_siz):

  if( itnts > i_min and itnts < rim1 ):
    # first block (black to red)
    uni = (itnts - i_min) / blk_siz
    red_out = (uni * 255.0)

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
  #not_serializable = '''
  if( itnts > i_min and itnts < rim1 ):
    # first block (black to red)
    green_out = 0

  elif( itnts >= rim1 and itnts <= rim2 ):
    # second block (red to yellow)
    uni = (itnts - rim1) / blk_siz
    green_out = (uni * 255.0)

  elif( itnts > rim2 and itnts <= i_max ):
    # third block (yellow to white)
    green_out = 255
  else:
    # defaulting to black
    green_out = 0

  return green_out


def blue_func(itnts, i_min, rim1, rim2, i_max, blk_siz):

  if( itnts > i_min  and itnts <= rim2 ):
    # first block (black to red)
    # or
    # second block (red to yellow)
    blue_out = 0

  elif( itnts > rim2 and itnts <= i_max ):
    # third block (yellow to white)
    uni = (itnts - rim2) / blk_siz
    blue_out = (uni * 255.0)

  else:
    # defaulting to black
    blue_out = 0

  return blue_out


if(__name__ == "__main__"):

  width = 40
  height = 50

  a = np.arange(width * height).reshape(width, height)
  b = np.zeros( (width, height, 3), 'uint8')

  b_red   = np.zeros( (width, height), 'uint8')
  b_green = np.zeros( (width, height), 'uint8')
  b_blue  = np.zeros( (width, height), 'uint8')

  print ("a =")
  print (a)
  print ("b =")
  print (b)

  b_red[:,:]   = a[:,:]
  b_green[:,:] = a[:,:]
  b_blue[:,:]  = a[:,:]



  print
  print ("after applying")
  print ("b[:,:,0] = a[:,:]")
  print ("b[:,:,1] = a[:,:]")
  print ("b[:,:,2] = a[:,:]")

  print ("b =")
  print (b)

  time1 = time.time()

  ints_min = np.min(a)
  ints_max = np.max(a)

  block_size = float(ints_max - ints_min) / 3.0
  rim_n1 = ints_min + block_size
  rim_n2 = ints_max - block_size

  '''

  for x in np.nditer(b_red[:,:], op_flags=['readwrite'], flags=['external_loop']):
    x[...] = red_func(x, ints_min, rim_n1, rim_n2, ints_max, block_size)

  for x in np.nditer(b_green[:,:], op_flags=['readwrite'], flags=['external_loop']):
    x[...] = green_func(x, ints_min, rim_n1, rim_n2, ints_max, block_size)

  for x in np.nditer(b_blue[:,:], op_flags=['readwrite'], flags=['external_loop']):
    x[...] = blue_func(x, ints_min, rim_n1, rim_n2, ints_max, block_size)

  '''

  for x in np.nditer(b_red[:,:], op_flags=['readwrite']):
    x[...] = red_func(x, ints_min, rim_n1, rim_n2, ints_max, block_size)

  for x in np.nditer(b_green[:,:], op_flags=['readwrite']):
    x[...] = green_func(x, ints_min, rim_n1, rim_n2, ints_max, block_size)

  for x in np.nditer(b_blue[:,:], op_flags=['readwrite']):
    x[...] = blue_func(x, ints_min, rim_n1, rim_n2, ints_max, block_size)


  b[:,:,0] = b_red[:,:]
  b[:,:,1] = b_green[:,:]
  b[:,:,2] = b_blue[:,:]


  time2 = time.time()
  time_dif = time2 - time1

  print
  print ("<<<< after"), (time_dif), ("time >>>>")
  print ("a =")
  print (a)
  print ("b =")
  print (b)
  print ("<<<< after"), (time_dif), ("time >>>>")
