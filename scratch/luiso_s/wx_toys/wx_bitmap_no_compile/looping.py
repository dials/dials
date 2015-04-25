import numpy as np

def loops_2d(data2d_scale):

  width = np.size( data2d_scale[0:1, :] )
  height = np.size( data2d_scale[:, 0:1] )
  img_array = np.empty( (height ,width, 3), 'int')

  img_array_r = np.empty( (height, width), 'int')
  img_array_g = np.empty( (height, width), 'int')
  img_array_b = np.empty( (height, width), 'int')

  scaled_i = np.empty( (height, width), 'int')

  red_byte = np.empty( (255 * 3), 'int')
  green_byte = np.empty( (255 * 3), 'int')
  blue_byte = np.empty( (255 * 3), 'int')

  for i in xrange(255):
    red_byte[i] = i
    green_byte[i + 255] = i
    blue_byte[i + 255 * 2 ] = i


  red_byte[255:255 * 3] = 255
  green_byte[0:255] = 0
  green_byte[255 * 2 : 255 * 3] = 255
  blue_byte[0:255 * 2] = 0

  blue_byte[764] = 255
  red_byte[764] = 255
  green_byte[764] = 255

  scaled_i[:,:] = data2d_scale[:,:]

  img_array_r[:,:] = scaled_i[:,:]
  for x in np.nditer(img_array_r[:,:], op_flags=['readwrite'], flags=['external_loop']):
    x[...] = red_byte[x]

  img_array_g[:,:] = scaled_i[:,:]
  for x in np.nditer(img_array_g[:,:], op_flags=['readwrite'], flags=['external_loop']):
    x[...] = green_byte[x]

  img_array_b[:,:] = scaled_i[:,:]
  for x in np.nditer(img_array_b[:,:], op_flags=['readwrite'], flags=['external_loop']):
    x[...] = blue_byte[x]

  img_array[:, :, 0] = img_array_r[:,:]
  img_array[:, :, 1] = img_array_g[:,:]
  img_array[:, :, 2] = img_array_b[:,:]
  print ("From Python call")

  return img_array
