#!/usr/bin/env python
#
# import.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dxtbx.datablock import DataBlockFactory, DataBlockDumper

help_message = '''

This program is used to import image data files into a format that can be used
within dials. The program looks at the metadata for each image along with the
filenames to determine the relationship between sets of images. Once all the
images have been analysed, a datablock object is written to file which specifies
the relationship between files. For example if two sets of images which belong
to two rotation scans have been given, two image sweeps will be saved. Images to
be processed are specified as command line arguments. Sometimes, there is a
maximum number of arguments that can be given on the command line and the number
of files may exceed this. In this case image filenames can be input on stdin
delimited by a new line using the -i option (see below for examples).
Alternatively a template can be specified using the template= parameter where
the consecutive digits representing the image numbers in the filenames are
replaced with '#' characters.

Examples::

  dials.import image_*.cbf

  dials.import image_1_*.cbf image_2_*.cbf

  dials.import template=image_1_####.cbf

  find . -name "image_*.cbf" | dials.import

  dials.import <<eof
  image_1.cbf
  image_2.cbf
  eof

'''

# Create the phil parameters
from libtbx.phil import parse
phil_scope = parse('''

  output = datablock.json
    .type = str
    .help = "The output JSON or pickle file"

  compact = False
    .type = bool
    .help = "For JSON output use compact representation"

  input {
    template = None
      .type = str
      .help = "The image sweep template"
      .multiple = True
    reference_geometry = None
      .type = path
      .help = "Experimental geometry from this datablock.json or "
              "experiments.json will override the geometry from the "
              "image headers."
  }

  mosflm_beam_centre = None
    .type = floats(size=2)
    .help = "Override the beam centre from the image headers, following "
            "the mosflm convention."

  lookup {
    mask = None
      .type = str
      .help = "Apply a mask to the imported data"

    gain = None
      .type = str
      .help = "Apply a gain to the imported data"

    pedestal = None
      .type = str
      .help = "Apply a pedestal to the imported data"
  }
''')


class Script(object):
  ''' Class to parse the command line options. '''

  def __init__(self):
    ''' Set the expected options. '''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # Create the option parser
    usage = "usage: %s [options] /path/to/image/files" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      sort_options=True,
      phil=phil_scope,
      read_datablocks_from_images=True,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    from dxtbx.datablock import DataBlockTemplateImporter
    from dials.util.options import flatten_datablocks
    import cPickle as pickle

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    datablocks = flatten_datablocks(params.input.datablock)

    # Load reference geometry
    reference_detector = None
    reference_beam = None
    if params.input.reference_geometry is not None:
      from dxtbx.serialize import load
      try:
        experiments = load.experiment_list(
          params.input.reference_geometry, check_format=False)
        assert len(experiments.detectors()) == 1
        assert len(experiments.beams()) == 1
        reference_detector = experiments.detectors()[0]
        reference_beam = experiments.beams()[0]
      except Exception, e:
        datablock = load.datablock(params.input.reference_geometry)
        assert len(datablock) == 1
        imageset = datablock[0].extract_imagesets()[0]
        reference_detector = imageset.get_detector()
        reference_beam = imageset.get_beam()

    # Check we have some filenames
    if len(datablocks) == 0:

      # Check if a template has been set and print help if not, otherwise try to
      # import the images based on the template input
      if len(params.input.template) == 0:
        self.parser.print_help()
        exit(0)
      else:
        importer = DataBlockTemplateImporter(
          params.input.template,
          options.verbose)
        datablocks = importer.datablocks

    # Check the lookup inputs
    mask_filename = None
    gain_filename = None
    dark_filename = None
    mask = None
    gain = None
    dark = None
    lookup_size = None
    if params.lookup.mask is not None:
      mask_filename = params.lookup.mask
      mask = pickle.load(open(mask_filename))
      if not isinstance(mask, tuple):
        mask = (mask,)
      lookup_size = [m.all() for m in mask]
    if params.lookup.gain is not None:
      gain_filename = params.lookup.gain
      gain = pickle.load(open(gain_filename))
      if not isinstance(gain, tuple):
        gain = (gain,)
      if lookup_size is None:
        lookup_size = [g.all() for g in gain]
      else:
        assert len(gain) == len(lookup_size), "Incompatible size"
        for s, g in zip(lookup_size, gain):
          assert s == g.all(), "Incompatible size"
    if params.lookup.pedestal is not None:
      dark_filename = params.lookup.pedestal
      dark = pickle.load(open(dark_filename))
      if not isinstance(dark, tuple):
        dark = (dark,)
      if lookup_size is None:
        lookup_size = [d.all() for d in dark]
      else:
        assert len(dark) == len(lookup_size), "Incompatible size"
        for s, d in zip(lookup_size, dark):
          assert s == d.all(), "Incompatible size"

    # Loop through the data blocks
    for i, datablock in enumerate(datablocks):

      # Extract any sweeps
      sweeps = datablock.extract_sweeps()

      # Extract any stills
      stills = datablock.extract_stills()
      if not stills:
        num_stills = 0
      else:
        num_stills = len(stills)

      # Set the external lookups
      if lookup_size is not None:
        for imageset in sweeps + stills:
          d = imageset.get_detector()
          assert len(lookup_size) == len(d), "Incompatible size"
          for s, p in zip(lookup_size, d):
            assert s == p.get_image_size()[::-1], "Incompatible size"
          if mask_filename is not None:
            imageset.external_lookup.mask.data = mask
            imageset.external_lookup.mask.filename = mask_filename
          if gain_filename is not None:
            imageset.external_lookup.gain.data = gain
            imageset.external_lookup.gain.filename = gain_filename
          if dark_filename is not None:
            imageset.external_lookup.pedestal.data = dark
            imageset.external_lookup.pedestal.filename = dark_filename

      if reference_beam is not None and reference_detector is not None:
        for sweep in sweeps:
          assert reference_detector.is_similar_to(sweep.get_detector(),static_only=True)
          sweep.set_beam(reference_beam)
          sweep.set_detector(reference_detector)

        for still in stills:
          assert reference_detector.is_similar_to(still.get_detector())
          still.set_beam(reference_beam)
          still.set_detector(reference_detector)

      elif params.mosflm_beam_centre is not None:
        from dxtbx.model.detector_helpers import set_mosflm_beam_centre
        for sweep in sweeps:
          set_mosflm_beam_centre(
            sweep.get_detector(), sweep.get_beam(), params.mosflm_beam_centre)

        for still in stills:
          set_mosflm_beam_centre(
            still.get_detector(), still.get_beam(), params.mosflm_beam_centre)

      # Print some data block info
      print "-" * 80
      print "DataBlock %d" % i
      print "  format: %s" % str(datablock.format_class())
      print "  num images: %d" % datablock.num_images()
      print "  num sweeps: %d" % len(sweeps)
      print "  num stills: %d" % num_stills

      # Loop through all the sweeps
      if options.verbose > 1:
        for j, sweep in enumerate(sweeps):
          print ""
          print "Sweep %d" % j
          print "  length %d" % len(sweep)
          print sweep.get_beam()
          print sweep.get_goniometer()
          print sweep.get_detector()
          print sweep.get_scan()

    # Write the datablock to a JSON or pickle file
    if params.output:
      print "-" * 80
      print 'Writing datablocks to %s' % params.output
      dump = DataBlockDumper(datablocks)
      dump.as_file(params.output, compact=params.compact)

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
