#!/usr/bin/env python
#
# asic_stats.py
#
#  Copyright (C) 2013 Lawrence Berkeley National Laboratory (LBNL)
#
#  Author: Aaron Brewster
#
#  This code is distributed under the X license, a copy of which is
#  included in the root directory of this package.
#
from __future__ import division
from scitbx.matrix import col
from scitbx.array_family import flex

help_message = '''

This program is used to calculate statisical measurements of consistency
within CSPAD ASICs after ASIC refinement.

See Hattne et. al. (2014), figure 2a

Examples:

  libtbx.python asic_stats.py refasics.json

  libtbx.python asic_stats.py datablock.json

'''

class Script(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env
    # Create the phil parameters
    phil_scope = parse('''

    ''')

    # Create the option parser
    usage = "usage: %s [options] /path/to/refined/json/file" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      sort_options=True,
      phil=phil_scope,
      read_experiments=True,
      read_datablocks=True,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    from dials.util.options import flatten_experiments, flatten_datablocks
    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)
    datablocks = flatten_datablocks(params.input.datablock)

    # Verify inputs
    if len(experiments) == 0 and len(datablocks) == 0:
      print "No experiments found"
      return

    if len(experiments) > 0 and len(datablocks) > 0:
      print "Please analyze only one datablock or experiment list"
      return

    # Find the detector of interest
    if len(experiments) > 0:
      print "Found %d experiments, doing analysis on experiment 0"%len(experiments)
      detector = experiments[0].detector
    else:
      print "Found %d datablocks, doing analysis on datablock 0"%len(datablocks)
      detector = datablocks[0].unique_detectors()[0]

    # Build a list of all 2x1 sensors in the detector
    sensors = []
    def recursive_find_sensors(group):
      if len(group) == 2:
        # Verify the children of this group are Panels, not PanelGroups
        assert [not hasattr(child, "children") for child in group].count(False) == 0
        sensors.append(group)

      else:
        for child in group:
          recursive_find_sensors(child)

    recursive_find_sensors(detector.hierarchy())

    print "Found %d 2x1 sensors"%len(sensors)

    pixel_gaps = flex.double()
    bottom_gaps = flex.double()
    angles = flex.double()
    sensor_angles = flex.double()

    for sensor in sensors:
      a = sensor[0]
      b = sensor[1]

      # Gget initial vectors, dropping the z axis
      center_a = col(a.get_local_origin()[0:2])
      center_b = col(b.get_local_origin()[0:2])

      pixel_size_a = a.get_pixel_size()[0]
      pixel_size_b = b.get_pixel_size()[0]
      assert pixel_size_a == pixel_size_b
      pixel_size = pixel_size_a

      width_a = a.get_image_size()[0]
      width_b = b.get_image_size()[0]
      assert width_a == width_b
      width = width_a * pixel_size

      # Calculate statistics
      pixel_gaps.append((abs((center_a-center_b).dot(col((1,0)))) - width) / pixel_size)
      bottom_gaps.append(((center_a-center_b).dot(col((0,1))) - 0) / pixel_size)

      slow_a = col(a.get_slow_axis())
      slow_b = col(b.get_slow_axis())

      angles.append(slow_a.angle(slow_b, deg=True))

      a_to_b = center_b - center_a
      sensor_angles.append(a_to_b.angle(col((1,0)), deg=True))

    pg_stats = flex.mean_and_variance(pixel_gaps)
    bg_stats = flex.mean_and_variance(bottom_gaps)
    an_stats = flex.mean_and_variance(angles)
    sa_stats = flex.mean_and_variance(sensor_angles)

    print "Sensor stats (means and standard deviations)"
    print "3 pixel gap                        : %f, %f"%(pg_stats.mean(), pg_stats.unweighted_sample_standard_deviation())
    print "Vertical offset                    : %f, %f"%(bg_stats.mean(), bg_stats.unweighted_sample_standard_deviation())
    print "Angular deviations (between asics) : %f, %f"%(an_stats.mean(), an_stats.unweighted_sample_standard_deviation())
    print "Angular deviations (sensor)        : %f, %f"%(sa_stats.mean(), an_stats.unweighted_sample_standard_deviation())


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

