from __future__ import absolute_import, division

def display(shoebox):

  from matplotlib import pylab
  from math import sqrt, ceil
  n = shoebox.shape[0]

  ncol = int(ceil(sqrt(n)))
  nrow = int(ceil(n / ncol))
  vmin = min(shoebox.reshape(-1))
  vmax = max(shoebox.reshape(-1))
  #if vmax < 10:
  #    vmax = 10
  for i in range(n):
    pylab.subplot(nrow, ncol, i+1)
    pylab.imshow(shoebox[i], vmin=0, vmax=vmax, interpolation='none')
  pylab.show()

def display2(xyz1, xyz2, xyz3):

  ymax = max([max(xyz1), max(xyz2), max(xyz3)])

  from matplotlib import pylab
  pylab.subplot(3, 1, 1)
  pylab.title('XDS Obs - XDS Cal')
  pylab.ylim([0, ymax])
  pylab.scatter(range(len(xyz1)), xyz1, color='blue')
  pylab.subplot(3, 1, 2)
  pylab.title('My Obs - XDS Cal')
  pylab.ylim([0, ymax])
  pylab.scatter(range(len(xyz2)), xyz2, color='blue')
  pylab.subplot(3, 1, 3)
  pylab.title('My Obs - XDS Obs')
  pylab.ylim([0, ymax])
  pylab.scatter(range(len(xyz3)), xyz3, color='blue')
  pylab.show()

if __name__ == '__main__':

  from dials.algorithms.image.centroid import CentroidMaskedImage3d
  from scitbx import matrix
  from scitbx.array_family import flex
  import cPickle as pickle
  from dials.algorithms.background import NormalDiscriminator
  from dials.algorithms.background import MeanSubtractor

  filename = "/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data/integrate_hkl.pickle"
  reflections = pickle.load(open(filename, 'r'))

  discriminate = NormalDiscriminator()
  subtract = MeanSubtractor()


  xyz_diff_xds_cal = []
  xyz_diff_obs_cal = []
  xyz_diff_obs_xds = []

  for r in reflections:

    if matrix.col(r.centroid_position).length() > 0:

      max_val = max(r.shoebox)

      r.shoebox_mask = discriminate(r.shoebox)

      #display(r.shoebox.as_numpy_array())
      #display(r.shoebox_mask.as_numpy_array())

      subtract(r)

      shoebox = r.shoebox
      mask = (r.shoebox_mask == 2).as_1d().as_int()
      mask.reshape(flex.grid(shoebox.all()))

      #centroid = CentroidImage3d(shoebox.as_double())
      try:
        centroid = CentroidMaskedImage3d(shoebox.as_double(), mask)
      except Exception:
        continue

      bbox = r.bounding_box
      offset = matrix.col((bbox[0], bbox[2], bbox[4]))
      xyz1 = matrix.col(r.centroid_position)
      xyzcal = matrix.col(r.image_coord_px + (r.frame_number,))
      xyz2 = offset + matrix.col(centroid.mean())
      hkl = r.miller_index

      str1 = '({0:.2f}, {1:.2f}, {2:.2f})'.format(*xyzcal)
      str2 = '({0:.2f}, {1:.2f}, {2:.2f})'.format(*xyz1)
      str3 = '({0:.2f}, {1:.2f}, {2:.2f})'.format(*xyz2)
      print '{0:>16} {1:>24} {2:>24} {3:>24} {4:>6}'.format(hkl, str1, str2, str3, max_val)

      xyz_diff_xds_cal.append((xyz1 - xyzcal).length())
      xyz_diff_obs_cal.append((xyz2 - xyzcal).length())
      xyz_diff_obs_xds.append((xyz2 - xyz1).length())



  #display2(xyz_diff_xds_cal, xyz_diff_obs_cal, xyz_diff_obs_xds)

    #display(shoebox)
