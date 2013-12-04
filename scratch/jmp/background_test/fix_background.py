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
      x = -0.5
      y = -0.5
      z = -0.5
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
        z = 1
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
#  background = XdsSubtractor()
#  background(None, None, rlist)
  from dials.algorithms.shoebox import MaskCode
  from scitbx.array_family import flex
  for r in rlist:
    mask = r.shoebox_mask
    background = flex.bool([bool(m & MaskCode.Background) for m in mask])
    pixels = r.shoebox.select(background)
    bg = flex.median(pixels)
    print bg
    for i in range(len(r.shoebox_background)):
      r.shoebox_background[i] = bg

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

  background_xds(rlist)
  integrate_3d_summation(rlist)

  for r, e in zip(rlist, expected):
    mask = r.shoebox_mask
    background = flex.bool([bool(m & MaskCode.Background) for m in mask])
    foreground = flex.bool([bool(m & MaskCode.Foreground) for m in mask])
    pixels = r.shoebox.select(foreground)


    print e, r.intensity, flex.sum(pixels) - len(pixels) * 10.0, e / r.intensity

#  for r in rlist:
#    from math import sqrt, erf
#    r.intensity *= 1.0 + (1.0 - erf(3.0 / sqrt(2.0)))


  from matplotlib import pylab
  pylab.scatter(expected, [e / r.intensity for r, e in zip(rlist, expected)])
  pylab.show()


if __name__ == '__main__':
  generate_reflections(100, 0, 10000, 10)
