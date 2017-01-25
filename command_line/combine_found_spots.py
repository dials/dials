#!/usr/bin/env python
#
# dials.combine_found_spots.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
# DIALS_ENABLE_COMMAND_LINE_COMPLETION

from __future__ import absolute_import, division

import libtbx.load_env
import logging
logger = logging.getLogger(libtbx.env.dispatcher_name)

help_message = '''


'''

# Set the phil scope
from libtbx.phil import parse
phil_scope = parse('''

  output {
    reflections = 'combined_strong.pickle'
      .type = str
      .help = "The output filename"

    datablock = 'combined_datablock.json'
      .type = str
      .help = "Save the modified datablock."

    log = 'dials.combine_found_spots.log'
      .type = str
      .help = "The log filename"

    debug_log = 'dials.combine_found_spots.debug.log'
      .type = str
      .help = "The debug log filename"
  }

  input {
    tolerance
        .help = "Tolerances used to determine shared models"
      {

      beam {

        wavelength = 1e-6
          .type = float(value_min=0.0)
          .help = "The wavelength tolerance"

        direction = 1e-6
          .type = float(value_min=0.0)
          .help = "The direction tolerance"

        polarization_normal = 1e-6
          .type = float(value_min=0.0)
          .help = "The polarization normal tolerance"

        polarization_fraction = 1e-6
          .type = float(value_min=0.0)
          .help = "The polarization fraction tolerance"

      }

      detector {

        fast_axis = 1e-6
          .type = float(value_min=0.0)
          .help = "The fast axis tolerance"

        slow_axis = 1e-6
          .type = float(value_min=0.0)
          .help = "The slow axis tolerance"

        origin = 1e-3
          .type = float(value_min=0.0)
          .help = "The origin tolerance"

      }

      goniometer {

        rotation_axis = 1e-6
          .type = float(value_min=0.0)
          .help = "The rotation axis tolerance"

        fixed_rotation = 1e-6
          .type = float(value_min=0.0)
          .help = "The fixed rotation tolerance"

        setting_rotation = 1e-6
          .type = float(value_min=0.0)
          .help = "The setting rotation tolerance"

      }

      scan {

        oscillation = 0.01
          .type = float(value_min=0.0)
          .help = "The oscillation tolerance for the scan"

      }
    }
  }

  verbosity = 1
    .type = int(value_min=0)
    .help = "The verbosity level"

''', process_includes=True)


def combine(datablock_list, reflections_list, params):
  '''
  Combine the found spots.

  '''
  from dxtbx.datablock import BeamComparison
  from dxtbx.datablock import DetectorComparison
  from dxtbx.datablock import GoniometerComparison
  from dxtbx.datablock import DataBlock
  from dxtbx.imageset import ImageSetFactory
  from dials.algorithms.spot_finding import StrongSpotCombiner
  from dials.array_family import flex
  assert len(datablock_list) == len(reflections_list)

  # Get a list of imagesets
  imageset_list = []
  for db in datablock_list:
    iset = db.extract_imagesets()
    assert len(iset) == 1
    imageset_list.append(iset[0])

  compare_beam = BeamComparison(
    wavelength_tolerance=params.input.tolerance.beam.wavelength,
    direction_tolerance=params.input.tolerance.beam.direction,
    polarization_normal_tolerance=params.input.tolerance.beam.polarization_normal,
    polarization_fraction_tolerance=params.input.tolerance.beam.polarization_fraction)
  compare_detector = DetectorComparison(
    fast_axis_tolerance=params.input.tolerance.detector.fast_axis,
    slow_axis_tolerance=params.input.tolerance.detector.slow_axis,
    origin_tolerance=params.input.tolerance.detector.origin)
  compare_goniometer = GoniometerComparison(
    rotation_axis_tolerance=params.input.tolerance.goniometer.rotation_axis,
    fixed_rotation_tolerance=params.input.tolerance.goniometer.fixed_rotation,
    setting_rotation_tolerance=params.input.tolerance.goniometer.setting_rotation)
  scan_tolerance = params.input.tolerance.scan.oscillation

  # The initial models
  format_class = imageset_list[0].reader().get_format_class()
  beam = imageset_list[0].get_beam()
  detector = imageset_list[0].get_detector()
  goniometer = imageset_list[0].get_goniometer()
  scan = imageset_list[0].get_scan()
  template = imageset_list[0].get_template()

  # Check all the models
  for imageset in imageset_list[1:]:
    b = imageset.get_beam()
    d = imageset.get_detector()
    g = imageset.get_goniometer()
    s = imageset.get_scan()
    if not imageset.reader().get_format_class() == format_class:
      raise RuntimeError('Format classes do not match')
    if not imageset.get_template() == template:
      raise RuntimeError('Templates do not match')
    if not compare_beam(beam, b):
      raise RuntimeError('Beam models are too dissimilar')
    if not compare_detector(detector, d):
      raise RuntimeError('Detector models are too dissimilar')
    if not compare_goniometer(goniometer, g):
      raise RuntimeError('Goniometer models are too dissimilar')
    try:
      scan.append(s, scan_tolerance=scan_tolerance)
    except Exception:
      raise RuntimeError('Scans do not match')

  # Get the image range
  image_range = scan.get_image_range()
  image_range = (image_range[0], image_range[1]+1)

  # Create the sweep
  imageset = ImageSetFactory.make_sweep(
    template, range(*image_range),
    format_class,
    beam, detector,
    goniometer, scan)

  # Combine spots
  combiner = StrongSpotCombiner()
  for index, rlist in enumerate(reflections_list, start=1):
    assert rlist['id'].all_eq(0)
    logger.info("Combining %d reflections from reflection list %d" % (
      len(rlist),
      index))
    combiner.add(rlist['shoebox'])
  shoeboxes = combiner.shoeboxes()

  # Calculate the spot centroids and intensities
  logger.info('Combined into %d reflections' % len(shoeboxes))
  centroid = shoeboxes.centroid_valid()
  logger.info('Calculated {0} spot centroids'.format(len(shoeboxes)))
  intensity = shoeboxes.summed_intensity()
  logger.info('Calculated {0} spot intensities'.format(len(shoeboxes)))

  # Construct the reflection table
  reflections = flex.reflection_table(
    flex.observation(
      shoeboxes.panels(),
      centroid,
      intensity),
    shoeboxes)
  reflections['id'] = flex.int(len(reflections), 0)
  reflections.set_flags(
    flex.size_t_range(len(reflections)),
    reflections.flags.strong)

  # Return the datablock and reflections
  return DataBlock([imageset]), reflections


class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] [param.phil] "\
            "datablock.json" \
            % libtbx.env.dispatcher_name

    # Initialise the base class
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_datablocks=True,
      read_reflections=True)

  def run(self):
    '''Execute the script.'''
    from dials.array_family import flex
    from dials.util.options import flatten_datablocks
    from dials.util.options import flatten_reflections
    from time import time
    from dials.util import log
    from libtbx.utils import Sorry
    start_time = time()

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=False)

    # Configure the logging
    log.config(
      params.verbosity,
      info=params.output.log,
      debug=params.output.debug_log)

    from dials.util.version import dials_version
    logger.info(dials_version())

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      logger.info('The following parameters have been modified:\n')
      logger.info(diff_phil)

    # Ensure we have a data block
    datablocks = flatten_datablocks(params.input.datablock)
    reflections = flatten_reflections(params.input.reflections)
    if len(datablocks) == 0 and len(reflections) == 0:
      self.parser.print_help()
      return
    elif len(datablocks) != len(reflections):
      raise Sorry("Must have same number of datablocks and reflection tables")

    # Combine the datablocks and reflections
    datablock, reflections = combine(
      datablocks,
      reflections,
      params)

    # Save the reflections to file
    logger.info('\n' + '-' * 80)
    reflections.as_pickle(params.output.reflections)
    logger.info('Saved {0} reflections to {1}'.format(
        len(reflections), params.output.reflections))

    # Save the datablock
    from dxtbx.datablock import DataBlockDumper
    logger.info('Saving datablocks to {0}'.format(
      params.output.datablock))
    dump = DataBlockDumper(datablocks)
    dump.as_file(params.output.datablock)


    # Print the time
    logger.info("Time Taken: %f" % (time() - start_time))


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
