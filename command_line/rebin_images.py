# LIBTBX_SET_DISPATCHER_NAME dev.dials.rebin_images

from __future__ import absolute_import, division, print_function

def rebin_images(in_template, out_template, start, end):
  from dials.util.rebin_images import main
  in_images = [in_template % j for j in range(start, end + 1)]
  out_images = [out_template % j for j in range(start, end + 1)]
  main(in_images, out_images)
  return

if __name__ == '__main__':
  import sys
  if len(sys.argv) != 5:
    raise RuntimeError('%s in_\%04d.cbf out_\%04d.cbf start end' % \
      sys.argv[0])
  in_template = sys.argv[1]
  out_template = sys.argv[2]
  start = int(sys.argv[3])
  end = int(sys.argv[4])
  rebin_images(in_template, out_template, start, end)
