#!/usr/bin/env python
#
# dials.find_spots.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
# DIALS_ENABLE_COMMAND_LINE_COMPLETION

from __future__ import division

help_message = '''

This program tries to find strong spots on a sequence of images. The program can
be called with either a "datablock.json" file or a sequence of image files (see
help for dials.import for more information about how images are imported). Spot
finding will be done against each logically grouped set of images given. Strong
pixels will be found on each image and spots will be formed from connected
components. In the case of rotation images, connected component labelling will
be done in 3D.

Once a set of spots have been found, their centroids and intensities will be
calculated. They will then be filtered according to the particular preferences
of the user. The output will be a file (strong.pickle) containing a list of spot
centroids and intensities which can be used in the dials.index program. To view
a list of parameters for spot finding use the --show-config option.

Examples::

  dials.find_spots image1.cbf

  dials.find_spots imager_00*.cbf

  dials.find_spots datablock.json

  dials.find_spots datablock.json output.reflections=strong.pickle

'''

# Set the phil scope
from libtbx.phil import parse
phil_scope = parse('''

  output {
    reflections = 'strong.pickle'
      .type = str
      .help = "The output filename"

    shoeboxes = True
      .type = bool
      .help = "Save the raw pixel values inside the reflection shoeboxes."

    datablock = None
      .type = str
      .help = "Save the modified datablock."
              "(usually only modified with hot pixel mask)"

    log = 'dials.find_spots.log'
      .type = str
      .help = "The log filename"

    debug_log = 'dials.find_spots.debug.log'
      .type = str
      .help = "The debug log filename"
  }

  per_image_statistics = False
    .type = bool
    .help = "Whether or not to print a table of per-image statistics."

  verbosity = 1
    .type = int(value_min=0)
    .help = "The verbosity level"

  include scope dials.algorithms.peak_finding.spotfinder_factory.phil_scope

''', process_includes=True)


class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] [param.phil] "\
            "{datablock.json | image1.file [image2.file ...]}" \
            % libtbx.env.dispatcher_name

    # Initialise the base class
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_datablocks=True,
      read_datablocks_from_images=True)

  def run(self):
    '''Execute the script.'''
    from dials.array_family import flex
    from dials.util.options import flatten_datablocks
    from time import time
    from dials.util import log
    from logging import info, debug
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
    info(dials_version())

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      info('The following parameters have been modified:\n')
      info(diff_phil)

    # Ensure we have a data block
    datablocks = flatten_datablocks(params.input.datablock)
    if len(datablocks) == 0:
      self.parser.print_help()
      return
    elif len(datablocks) != 1:
      raise Sorry('only 1 datablock can be processed at a time')

    # Loop through all the imagesets and find the strong spots
    reflections = flex.reflection_table.from_observations(
      datablocks[0], params)

    # Delete the shoeboxes
    if not params.output.shoeboxes:
      del reflections['shoebox']

    # ascii spot count per image plot
    from dials.util.ascii_art import spot_counts_per_image_plot

    for i, imageset in enumerate(datablocks[0].extract_imagesets()):
      ascii_plot = spot_counts_per_image_plot(
        reflections.select(reflections['id'] == i))
      if len(ascii_plot):
        info('\nHistogram of per-image spot count for imageset %i:' %i)
        info(ascii_plot)

    # Save the reflections to file
    info('\n' + '-' * 80)
    reflections.as_pickle(params.output.reflections)
    info('Saved {0} reflections to {1}'.format(
        len(reflections), params.output.reflections))

    # Save the datablock
    if params.output.datablock:
      from dxtbx.datablock import DataBlockDumper
      info('Saving datablocks to {0}'.format(
        params.output.datablock))
      dump = DataBlockDumper(datablocks)
      dump.as_file(params.output.datablock)

    # Print some per image statistics
    if params.per_image_statistics:
      from dials.algorithms.peak_finding import per_image_analysis
      from cStringIO import StringIO
      s = StringIO()
      for i, imageset in enumerate(datablocks[0].extract_imagesets()):
        print >> s, "Number of centroids per image for imageset %i:" %i
        stats = per_image_analysis.stats_imageset(
          imageset, reflections.select(reflections['id'] == i),
          resolution_analysis=False)
        per_image_analysis.print_table(stats, out=s)
      info(s.getvalue())

    # Print the time
    info("Time Taken: %f" % (time() - start_time))


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
