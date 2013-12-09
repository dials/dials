from __future__ import division

def print_refl(refl):
  sbox = refl.shoebox
  mask = refl.shoebox_mask
  ndim = sbox.nd()
  dims = sbox.focus()

  from dials.algorithms.shoebox import MaskCode

  print '-' * 80

  print 'HKL:      %d %d %d' % refl.miller_index
  print 'S vector: %.3f %.3f %.3f' % refl.beam_vector
  print 'Bounding: %d %d %d %d %d %d' % refl.bounding_box
  print 'Frame:    %.3f' % refl.frame_number
  print 'Position: %.3f %.3f' % refl.image_coord_px

  print '-' * 80

  for k in range(dims[0]):
    for j in range(dims[1]):
      for i in range(dims[2]):
        if mask[k, j, i] & MaskCode.BackgroundUsed and \
          mask[k, j, i] & MaskCode.Background:
          print '%4d#' % int(sbox[k, j, i]),
        elif mask[k, j, i] & MaskCode.BackgroundUsed:
          print '%4d*' % int(sbox[k, j, i]),
        elif mask[k, j, i] & MaskCode.Background:
          print '%4d-' % int(sbox[k, j, i]),
        else:
          print '%4d ' % int(sbox[k, j, i]),
      print
    print '-' * 80

  return

if __name__ == '__main__':
  import cPickle as pickle
  import sys
  from dials.model.data import ReflectionList # implicit import

  for j, refl in enumerate(pickle.load(open(sys.argv[1]))):
    print '*' * 80
    print 'Reflection %d' % j
    print '*' * 80
    print_refl(refl)
