# LIBTBX_SET_DISPATCHER_NAME dev.dials.scale_down_image

from __future__ import absolute_import, division

# FIXME set a random seed here to allow results to be exactly reproduced.

if __name__ == '__main__':
  import sys
  if len(sys.argv) != 4:
    raise RuntimeError, '%s in_image.cbf out_image.cbf scale_factor' % \
      sys.argv[0]
  in_image = sys.argv[1]
  out_image = sys.argv[2]
  scale_factor = float(sys.argv[3])
  from dials.util.scale_down_image import scale_down_image
  scale_down_image(in_image, out_image, scale_factor)
