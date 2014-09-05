

if __name__ == '__main__':

  import sys
  from dials.array_family import flex
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory

  rtable = flex.reflection_table.from_pickle(sys.argv[1])
  mask = flex.bool([x == (0, 0, 0) for x in rtable['xyzcal.mm']])
  rtable.del_selected(mask)
  mask = flex.bool([h == (0, 0, 0) for h in rtable['miller_index']])
  rtable.del_selected(mask)
  rtable = rtable[0:5000]
  shoebox = rtable['shoebox']
  xcal, ycal, zcal = rtable['xyzobs.px.value'].parts()


  exlist = ExperimentListFactory.from_json_file(sys.argv[2])
  beam = exlist[0].beam
  detector = exlist[0].detector
  scan = exlist[0].scan

  from scitbx import matrix
  s0 = matrix.col(beam.get_s0())

  counts = []
  diff_angle_d = []
  diff_angle_m = []


  for index in range(len(rtable)):
    xc = xcal[index]
    yc = ycal[index]
    zc = zcal[index]
    sbox = shoebox[index]
    s1c = matrix.col(detector[0].get_pixel_lab_coord((xc, yc)))
    phic = scan.get_angle_from_array_index(zc)
    for j in range(sbox.data.all()[1]):
      for i in range(sbox.data.all()[2]):
        x = sbox.bbox[0] + i
        y = sbox.bbox[2] + j
        s1ij = matrix.col(detector[0].get_pixel_lab_coord((x,y)))
        da = s1c.angle(s1ij)
        for k in range(sbox.data.all()[0]):
          z = sbox.bbox[4] + k
          phi = scan.get_angle_from_array_index(z)
          db = phic - phi
          c = sbox.data[k,j,i]
          diff_angle_d.extend([da] * int(1000.0*c / flex.sum(sbox.data)))
          diff_angle_m.extend([db] * int(1000.0*c/ flex.sum(sbox.data)))
          # counts.append(sbox.data[k,j,i])
    print index

  print min(diff_angle_d), max(diff_angle_d)
  print min(diff_angle_m), max(diff_angle_m)

  m = sum(diff_angle_m) / len(diff_angle_m)
  v = sum([(d - m)**2 for d in diff_angle_m]) / len(diff_angle_m)
  from math import sqrt
  print m, sqrt(v)

  from matplotlib import pylab
  pylab.hist(diff_angle_d, bins=100)
  pylab.show()
  pylab.hist(diff_angle_m, bins=100)
  pylab.show()
