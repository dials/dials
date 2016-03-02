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
as shown in the examples below. Alternatively a template can be specified using
the template= parameter where the consecutive digits representing the image
numbers in the filenames are replaced with '#' characters.

The geometry can be set manually by either specifying a datablock.json file
contaning the reference geometry, by setting the mosflm beam centre of by
setting each variable to be overriden.

Examples::

  dials.import /data/directory-containing-images/

  dials.import image_*.cbf

  dials.import image_1_*.cbf image_2_*.cbf

  dials.import directory/with/images

  dials.import template=image_1_####.cbf

  dials.import directory=directory/with/images

  find . -name "image_*.cbf" | dials.import

  dials.import <<eof
  image_1.cbf
  image_2.cbf
  eof

'''

# Create the phil parameters
from libtbx.phil import parse
phil_scope = parse('''

  output {

    datablock = datablock.json
      .type = str
      .help = "The output JSON or pickle file"

    log = 'dials.import.log'
      .type = str
      .help = "The log filename"

    debug_log = 'dials.import.debug.log'
      .type = str
      .help = "The debug log filename"

    compact = False
      .type = bool
      .help = "For JSON output use compact representation"

  }

  verbosity = 1
    .type = int(value_min=0)
    .help = "The verbosity level"

  input {
    template = None
      .type = str
      .help = "The image sweep template"
      .multiple = True

    directory = None
      .type = str
      .help = "A directory with images"
      .multiple = True

    reference_geometry = None
      .type = path
      .help = "Experimental geometry from this datablock.json or "
              "experiments.json will override the geometry from the "
              "image headers."

    allow_multiple_sweeps = False
      .type = bool
      .help = "If False, raise an error if multiple sweeps are found"
  }

  geometry {

    beam {

      wavelength = None
        .type = float
        .help = "Override the beam wavelength"

      direction = None
        .type = floats(size=3)
        .help = "Override the beam direction"

    }

    detector {

      panel
        .multiple = True
      {
        id = 0
          .type = int
          .help = "The panel number"

        name = None
          .type = str
          .help = "Override the panel name"

        type = None
          .type = str
          .help = "Override the panel type"

        pixel_size = None
          .type = floats(size=2)
          .help = "Override the panel pixel size"

        image_size = None
          .type = ints(size=2)
          .help = "Override the panel image size"

        trusted_range = None
          .type = floats(size=2)
          .help = "Override the panel trusted range"

        thickness = None
          .type = float
          .help = "Override the panel thickness"

        material = None
          .type = str
          .help = "Override the panel material"

        fast_axis = None
          .type = floats(size=3)
          .help = "Override the panel fast axis. Requires slow_axis and origin."

        slow_axis = None
          .type = floats(size=3)
          .help = "Override the panel slow axis. Requires fast_axis and origin."

        origin = None
          .type = floats(size=3)
          .help = "Override the panel origin. Requires fast_axis and slow_axis."

      }

    }

    goniometer {

      rotation_axis = None
        .type = floats(size=3)
        .help = "Override the rotation axis"

      fixed_rotation = None
        .type = floats(size=9)
        .help = "Override the fixed rotation matrix"

      setting_rotation = None
        .type = floats(size=9)
        .help = "Override the setting rotation matrix"

    }

    scan {

      image_range = None
        .type = ints(size=2)
        .help = "Override the image range"

      extrapolate_scan = False
        .type = bool
        .help = "When overriding the image range, extrapolate exposure and epoch information from existing images"

      oscillation = None
        .type = floats(size=2)
        .help = "Override the image oscillation"

    }

    mosflm_beam_centre = None
      .type = floats(size=2)
      .help = "Override the beam centre from the image headers, following "
              "the mosflm convention."
  }

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
    from dxtbx.datablock import DataBlockFactory
    from dxtbx.datablock import DataBlockTemplateImporter
    from dials.util.options import flatten_datablocks
    from dials.util import log
    from logging import info, debug
    import cPickle as pickle
    from libtbx.utils import Sorry

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=False)
    datablocks = flatten_datablocks(params.input.datablock)

    # Configure logging
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
      if len(params.input.template) > 0:
        importer = DataBlockTemplateImporter(
          params.input.template,
          options.verbose)
        datablocks = importer.datablocks
      elif len(params.input.directory) > 0:
        datablocks = DataBlockFactory.from_filenames(params.input.directory)
      else:
        self.parser.print_help()
        exit(0)

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

    # Only allow a single datablock
    if len(datablocks) > 1:
      raise Sorry("More than 1 datablock found")

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

      for imageset in sweeps + stills:
        if imageset.get_beam() == None or imageset.get_detector() == None:
          raise Sorry('''
            Imageset contains no beam or detector model. This means you will be
            unable to process your data.

            Possible causes of this error are:
               - A problem reading the images with one of the dxtbx format classes
               - A lack of header information in the file itself.
          ''')

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

      # Override the geometry. Use the reference geometry first, if set.
      # Otherwise use the mosflm beam centre and finally look to see if
      # any items have been otherwise overridden
      if reference_beam is not None and reference_detector is not None:
        for sweep in sweeps:
          assert reference_detector.is_similar_to(
            sweep.get_detector(),
            static_only=True)
          sweep.set_beam(reference_beam)
          sweep.set_detector(reference_detector)
        for still in stills:
          assert reference_detector.is_similar_to(still.get_detector())
          still.set_beam(reference_beam)
          still.set_detector(reference_detector)
      elif params.geometry.mosflm_beam_centre is not None:
        from dxtbx.model.detector_helpers import set_mosflm_beam_centre
        for sweep in sweeps:
          set_mosflm_beam_centre(
            sweep.get_detector(),
            sweep.get_beam(),
            params.geometry.mosflm_beam_centre)
        for still in stills:
          set_mosflm_beam_centre(
            still.get_detector(),
            still.get_beam(),
            params.geometry.mosflm_beam_centre)
      else:
        def override_beam(beam, params):
          if params.wavelength is not None:
            beam.set_wavelength(params.wavelength)
          if params.direction is not None:
            beam.set_direction(params.direction)
        def override_detector(detector, params):
          # need to gather material from multiple phil parameters to set
          frame_hash = { }
          for panel_params in params.panel:
            panel_id = panel_params.id
            panel = detector[panel_params.id]

            if not panel_id in frame_hash:
              frame_hash[panel_id] = {'fast_axis':None,
                                      'slow_axis':None,
                                      'origin':None}

            if panel_params.fast_axis is not None:
              frame_hash[panel_id]['fast_axis'] = panel_params.fast_axis
            if panel_params.slow_axis is not None:
              frame_hash[panel_id]['slow_axis'] = panel_params.slow_axis
            if panel_params.origin is not None:
              frame_hash[panel_id]['origin'] = panel_params.origin

            if panel_params.name is not None:
              panel.set_name(panel_params.name)
            if panel_params.type is not None:
              panel.set_type(panel_params.type)
            if panel_params.pixel_size is not None:
              panel.set_pixel_size(panel_params.pixel_size)
            if panel_params.image_size is not None:
              panel.set_image_size(panel_params.image_size)
            if panel_params.trusted_range is not None:
              panel.set_trusted_range(panel_params.trusted_range)
            if panel_params.thickness is not None:
              panel.set_thickness(panel_params.thickness)
            if panel_params.material is not None:
              panel.set_material(panel_params.material)

          for panel_id in frame_hash:
            fast_axis = frame_hash[panel_id]['fast_axis']
            slow_axis = frame_hash[panel_id]['slow_axis']
            origin = frame_hash[panel_id]['origin']
            # FIXME warn user if not (all none or all not none)
            if (fast_axis is not None and
                slow_axis is not None and
                origin is not None):
              panel = detector[panel_id]
              panel.set_frame(fast_axis, slow_axis, origin)
        def override_goniometer(goniometer, params):
          if params.rotation_axis is not None:
            goniometer.set_rotation_axis(params.rotation_axis)
          if params.fixed_rotation is not None:
            goniometer.set_fixed_rotation(params.fixed_rotation)
          if params.setting_rotation is not None:
            goniometer.set_setting_rotation(params.setting_rotation)
        def override_scan(scan, params):
          if params.image_range is not None:
            most_recent_image_index = scan.get_image_range()[1] - scan.get_image_range()[0]
            scan.set_image_range(params.image_range)
            if params.extrapolate_scan and \
                (params.image_range[1] - params.image_range[0]) > most_recent_image_index:
              exposure_times = scan.get_exposure_times()
              epochs = scan.get_epochs()
              exposure_time = exposure_times[most_recent_image_index]
              epoch_correction = epochs[most_recent_image_index]
              for i in range(most_recent_image_index + 1, \
                  params.image_range[1] - params.image_range[0] + 1):
                exposure_times[i] = exposure_time
                epoch_correction += exposure_time
                epochs[i] = epoch_correction
              scan.set_epochs(epochs)
              scan.set_exposure_times(exposure_times)
          if params.oscillation is not None:
            scan.set_oscillation(params.oscillation)
        for sweep in sweeps:
          override_beam(sweep.get_beam(), params.geometry.beam)
          override_detector(sweep.get_detector(), params.geometry.detector)
          override_goniometer(sweep.get_goniometer(), params.geometry.goniometer)
          override_scan(sweep.get_scan(), params.geometry.scan)
        for still in stills:
          override_beam(still.get_beam(), params.geometry.beam)
          override_detector(still.get_detector(), params.geometry.detector)

      # Print some data block info
      info("-" * 80)
      info("DataBlock %d" % i)
      info("  format: %s" % str(datablock.format_class()))
      info("  num images: %d" % datablock.num_images())
      info("  num sweeps: %d" % len(sweeps))
      info("  num stills: %d" % num_stills)

      # Loop through all the sweeps
      for j, sweep in enumerate(sweeps):
        debug("")
        debug("Sweep %d" % j)
        debug("  Length %d" % len(sweep))
        debug(sweep.get_beam())
        debug(sweep.get_goniometer())
        debug(sweep.get_detector())
        debug(sweep.get_scan())

      # Only allow a single sweep
      if params.input.allow_multiple_sweeps is False:
        if len(sweeps) > 1:
          raise Sorry('''
            More than 1 sweep was found. Two things may be happening here:

            1. There really is more than 1 sweep. If you expected this to be the
               case, set the parameter allow_multiple_sweeps=True. If you don't
               expect this, then check the input to dials.import.

            2. There may be something wrong with your image headers (for example,
               the rotation ranges of each image may not match up). You should
               investigate what went wrong, but you can force dials.import to treat
               your images as a single sweep by using the template=image_####.cbf
               parameter (see help).
          ''')


    # Write the datablock to a JSON or pickle file
    if params.output.datablock:
      info("-" * 80)
      info('Writing datablocks to %s' % params.output.datablock)
      dump = DataBlockDumper(datablocks)
      dump.as_file(params.output.datablock, compact=params.output.compact)

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
