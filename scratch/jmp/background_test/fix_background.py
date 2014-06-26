from __future__ import division
def create_phil_string(n, cts, bg):
  template = """
    nrefl = %d
    shoebox_size {
      x = 10
      y = 10
      z = 10
    }
    spot_size {
      x = 1
      y = 1
      z = 1
    }
    spot_offset {
      x = 0
      y = 0
      z = 0
    }
    mask_nsigma = 3.0
    counts = %d
    background = %d
    pixel_mask = all *static precise
    background_method = *xds mosflm
    integration_methpd = *xds mosflm
    output {
      over = None
      under = None
      all = all_refl.pickle
    }
    rotation {
      axis {
        x = 0
        y = 0
        z = 0
      }
      angle = 0
    }
    """
  return template % (n, cts, bg)


def create_phil(n, cts, bg):
  from dials.algorithms.simulation.generate_test_reflections import \
    master_phil
  from libtbx.phil import command_line
  cmd = command_line.argument_interpreter(master_params = master_phil)
  working_phil = cmd.process_and_fetch(args = [create_phil_string(n, cts, bg)])
  return working_phil.extract()

def background_xds(rlist):
  from dials.algorithms.background import XdsSubtractor
  from math import sqrt, erf
  background = XdsSubtractor()
  background(None, None, rlist)
#  from dials.algorithms.shoebox import MaskCode
#  from scitbx.array_family import flex
#  for r in rlist:
#    mask = r.shoebox_mask
#    foreground = flex.bool([bool(m & MaskCode.Foreground) for m in mask])
#    pixels = r.shoebox.select(foreground)
#    I = flex.sum(pixels) * (1.0 - erf(3.0 / sqrt(2.0)))
#    for i in range(len(r.shoebox_background)):
#      r.shoebox_background[i] -= I
#  from dials.algorithms.shoebox import MaskCode
#  from scitbx.array_family import flex
#  for r in rlist:
#    mask = r.shoebox_mask
#    background = flex.bool([bool(m & MaskCode.Background) for m in mask])
#    pixels = r.shoebox.select(background)
#    bg = flex.median(pixels)
#    print bg
#    for i in range(len(r.shoebox_background)):
#      r.shoebox_background[i] = bg

  return

def background_inclined(rlist):
  from dials.algorithms.background import InclinedSubtractor
  background = InclinedSubtractor()
  background(None, None, rlist)
  return

def integrate_3d_summation(rlist):
  from dials.algorithms.integration import Summation3d
  integration = Summation3d()
  integration(None, None, rlist)
  return

def generate_reflections(num, min_cts, max_cts, bg):
  from dials.algorithms.simulation.generate_test_reflections import simple_gaussian_spots
  from random import randint
  from dials.model.data import ReflectionList
  from scitbx.array_family import flex
  from dials.algorithms.shoebox import MaskCode
  rlist = ReflectionList()
  expected = []
  for i in range(num):
    cts = randint(min_cts, max_cts)
    phil = create_phil(1, cts, bg)
    expected.append(cts)
    rlist.extend(simple_gaussian_spots(phil))

#  background_inclined(rlist)
  background_xds(rlist)
  integrate_3d_summation(rlist)

  import pickle
  pickle.dump(rlist, open("test.pickle", "w"))

#  for r in rlist:
#    from math import sqrt, erf
#    r.intensity *= 1.0 + (1.0 - erf(3.0 / sqrt(2.0)))


  for r, e in zip(rlist, expected):
    mask = r.shoebox_mask
    background = flex.bool([bool(m & MaskCode.Background) for m in mask])
    foreground = flex.bool([bool(m & MaskCode.Foreground) for m in mask])
    pixels = r.shoebox.select(foreground)
    print e, r.intensity, flex.sum(pixels) - len(pixels) * 0, e / r.intensity


  from math import sqrt
  I = [r.intensity for r in rlist]
  S = [sqrt(r.intensity_variance) for r in rlist]
  Z = [(i - e) / s for i, e, s in zip(I, expected, S)]
  meanz =sum(Z) / len(Z)
  sdevz =sqrt(sum((z - meanz)**2 for z in Z) / len(Z))
  print "MeanZ: %f, SDevZ: %f" % (meanz, sdevz)

  from matplotlib import pylab
  #pylab.ylim(0, 2)
  pylab.scatter(expected, [e / r.intensity for r, e in zip(rlist, expected)])
  pylab.show()


if __name__ == '__main__':
  generate_reflections(1000, 0, 10000, 10)
