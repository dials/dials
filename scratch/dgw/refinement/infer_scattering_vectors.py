#!/usr/bin/env python
#
# infer_scattering_vectors.py
#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""This script infers and sets the scattering vectors for a set of reflections
by projecting back from the detector intersection to the Ewald sphere."""

from __future__ import division

class ScriptRunner(object):
  """Class to run script."""

  def __init__(self, options, sweep_filename, reflections_filename):
    """Setup the script."""

    # Filename data
    self.output_filename = options.output_file
    self.sweep_filename = sweep_filename
    self.reflections_filename = reflections_filename

  def __call__(self):
    """Run the script."""
    from dials.model.serialize import load, dump
    import cPickle as pickle
    from scitbx import matrix

    # Load the reflection list
    print 'Loading reflections from {0}'.format(self.reflections_filename)
    rlist = pickle.load(open(self.reflections_filename, 'r'))

    # Try to load the models
    print 'Loading models from {0}'.format(self.sweep_filename)

    sweep = load.sweep(open(self.sweep_filename, 'r'))
    beam = sweep.get_beam()
    wavelength = beam.get_wavelength()
    detector = sweep.get_detector()

    # Loop through all the reflections
    for r in rlist:

      panel = detector[r.panel_number]
      lab_coord = panel.get_lab_coord(r.image_coord_mm)
      r.beam_vector = matrix.col(lab_coord).normalize() / wavelength

    # Write out reflections
    if self.output_filename is not None:

      print 'Saving reflections to {0}'.format(self.output_filename)
      pickle.dump(rlist, open(self.output_filename, 'wb'),
          pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':

  from optparse import OptionParser

  # Specify the command line options
  usage  = "usage: %prog [options] sweep.json reflections.pickle"

  # Create an option parser
  parser = OptionParser(usage)

  parser.add_option('-o', '--output-file',
                    dest='output_file', type="string", default=None,
                    help='Destination filename for reflections')

  # Parse the arguments
  options, args = parser.parse_args()

  # Print help if not enough arguments specified, otherwise call function
  if len(args) < 2:
    parser.print_help()
  else:
    # Initialise the script runner
    runner = ScriptRunner(options, sweep_filename=args[0],
                          reflections_filename=args[1])

    # Run the script
    runner()
