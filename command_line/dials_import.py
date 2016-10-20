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
# LIBTBX_SET_DISPATCHER_NAME dials.import
from __future__ import division
from dxtbx.datablock import DataBlockFactory, DataBlockDumper
from libtbx.utils import Sorry

import libtbx.load_env
import logging
logger = logging.getLogger(libtbx.env.dispatcher_name)

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
containing the reference geometry, by setting the mosflm beam centre or by
setting each variable to be overridden.

Examples::

  dials.import /data/directory-containing-images/

  dials.import image_*.cbf

  dials.import image_1_*.cbf image_2_*.cbf

  dials.import directory/with/images

  dials.import template=image_1_####.cbf

  dials.import directory=directory/with/images

  find . -name "image_*.cbf" | dials.import

  dials.import << EOF
  image_1.cbf
  image_2.cbf
  EOF

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

    as_grid_scan = False
      .type = bool
      .help = "Import as grid scan"

    grid_size = None
      .type = ints(size=2)
      .help = "If importing as a grid scan set the size"
  }

  include scope dials.util.options.geometry_phil_scope

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
''', process_includes=True)


class DataBlockImporter(object):
  '''
  A class to manage the import of the datablocks

  '''

  def __init__(self, params):
    '''
    Init the class

    '''
    self.params = params

  def __call__(self):
    '''
    Import the datablocks

    '''
    from dxtbx.datablock import DataBlockTemplateImporter
    from dxtbx.datablock import DataBlockFactory
    from dials.util.options import flatten_datablocks
    from libtbx.utils import Sorry

    # Get the datablocks
    datablocks = flatten_datablocks(self.params.input.datablock)

    # Check we have some filenames
    if len(datablocks) == 0:

      # Check if a template has been set and print help if not, otherwise try to
      # import the images based on the template input
      if len(self.params.input.template) > 0:
        importer = DataBlockTemplateImporter(
          self.params.input.template,
          max(self.params.verbosity-1, 0))
        datablocks = importer.datablocks
        if len(datablocks) == 0:
          raise Sorry('No datablocks found matching template %s' % self.params.input.template)
      elif len(self.params.input.directory) > 0:
        datablocks = DataBlockFactory.from_filenames(
          self.params.input.directory,
          max(self.params.verbosity-1, 0))
        if len(datablocks) == 0:
          raise Sorry('No datablocks found in directories %s' % self.params.input.directory)
      else:
        raise Sorry('No datablocks found')
    if len(datablocks) > 1:
      raise Sorry("More than 1 datablock found")

    # Return the datablocks
    return datablocks[0]


class ReferenceGeometryUpdater(object):
  '''
  A class to replace beam + detector with a reference

  '''

  def __init__(self, params):
    '''
    Load the reference geometry

    '''
    self.reference = self.load_reference_geometry(params)

  def __call__(self, imageset):
    '''
    Replace with the reference geometry

    '''
    # Check static detector items are the same
    assert self.reference.detector.is_similar_to(
      imageset.get_detector(),
      static_only=True)

    # Set beam and detector
    imageset.set_beam(self.reference.beam)
    imageset.set_detector(self.reference.detector)
    return imageset

  def load_reference_geometry(self, params):
    '''
    Load a reference geoetry file

    '''
    from collections import namedtuple
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
    Reference = namedtuple("Reference", ["detector", "beam"])
    return Reference(detector=reference_detector, beam=reference_beam)


class MosflmBeamCenterUpdater(object):
  '''
  A class to replace geometry with mosflm beam centre

  '''
  def __init__(self, params):
    '''
    Set the params

    '''
    self.params = params

  def __call__(self, imageset):
    '''
    Replace the geometry

    '''
    from dxtbx.model.detector_helpers import set_mosflm_beam_centre
    set_mosflm_beam_centre(
      imageset.get_detector(),
      imageset.get_beam(),
      self.params.geometry.mosflm_beam_centre)
    return imageset


class ManualGeometryUpdater(object):
  '''
  A class to update the geometry manually

  '''
  def __init__(self, params):
    '''
    Save the params

    '''
    self.params = params

  def __call__(self, imageset):
    '''
    Override the parameters

    '''
    from dxtbx.imageset import ImageSweep
    if isinstance(imageset, ImageSweep):
      self.override_beam(
        imageset.get_beam(),
        self.params.geometry.beam)
      self.override_detector(
        imageset.get_detector(),
        self.params.geometry.detector,
        imageset.get_beam(),
        self.params.geometry.beam)
      self.override_goniometer(
        imageset.get_goniometer(),
        self.params.geometry.goniometer)
      self.override_scan(
        imageset.get_scan(),
        self.params.geometry.scan)
    elif not self.params.geometry.scan.convert_stills_to_sweeps:
      for i in range(len(imageset)):
        self.override_beam(
          imageset.get_beam(index=i),
          self.params.geometry.beam)
        self.override_detector(
          imageset.get_detector(index=i),
          self.params.geometry.detector,
          imageset.get_beam(index=i),
          self.params.geometry.beam)
        self.override_goniometer(
          imageset.get_goniometer(index=i),
          self.params.geometry.goniometer)
        self.override_scan(
          imageset.get_scan(index=i),
          self.params.geometry.scan)
    else:
      imageset = self.convert_stills_to_sweep(imageset)
    return imageset

  def convert_stills_to_sweep(self, imageset):
    from dxtbx.model import Scan
    assert self.params.geometry.scan.oscillation is not None
    for i in range(len(imageset)):
      self.override_beam(
        imageset.get_beam(index=i),
        self.params.geometry.beam)
      self.override_detector(
        imageset.get_detector(index=i),
        self.params.geometry.detector)
      self.override_goniometer(
        imageset.get_goniometer(index=i),
        self.params.geometry.goniometer)
    beam = imageset.get_beam(index=0)
    detector = imageset.get_detector(index=0)
    goniometer = imageset.get_goniometer(index=0)
    for i in range(1, len(imageset)):
      assert beam.is_similar_to(
        imageset.get_beam(index=i),
        wavelength_tolerance            = self.params.input.tolerance.beam.wavelength,
        direction_tolerance             = self.params.input.tolerance.beam.direction,
        polarization_normal_tolerance   = self.params.input.tolerance.beam.polarization_normal,
        polarization_fraction_tolerance = self.params.input.tolerance.beam.polarization_fraction)
      assert detector.is_similar_to(
        imageset.get_detector(index=i),
        fast_axis_tolerance = self.params.input.tolerance.detector.fast_axis,
        slow_axis_tolerance = self.params.input.tolerance.detector.slow_axis,
        origin_tolerance    = self.params.input.tolerance.detector.origin)
      assert goniometer.is_similar_to(
        imageset.get_goniometer(index=i),
        rotation_axis_tolerance    = self.params.input.tolerance.goniometer.rotation_axis,
        fixed_rotation_tolerance   = self.params.input.tolerance.goniometer.fixed_rotation,
        setting_rotation_tolerance = self.params.input.tolerance.goniometer.setting_rotation)
    assert beam is not None
    assert detector is not None
    assert goniometer is not None
    if self.params.geometry.scan.image_range is not None:
      image_range = self.params.geometry.scan.image_range
    else:
      image_range = (1, len(imageset))
    oscillation = self.params.geometry.scan.oscillation
    scan = Scan(image_range=image_range, oscillation=oscillation)
    from dxtbx.sweep_filenames import template_regex
    from dxtbx.imageset import ImageSetFactory
    indices      = list(range(image_range[0], image_range[1]+1))
    template = template_regex(imageset.get_path(0))[0]
    if template is None:
      paths = [imageset.get_path(i) for i in range(len(imageset))]
      assert len(set(paths)) == 1
      template = paths[0]
    new_sweep = ImageSetFactory.make_sweep(
      template     = template,
      indices      = indices,
      format_class = imageset.reader().get_format_class(),
      beam         = beam,
      detector     = detector,
      goniometer   = goniometer,
      scan         = scan)
    return new_sweep

  def override_beam(self, beam, params):
    '''
    Override the beam parameters

    '''
    if params.wavelength is not None:
      beam.set_wavelength(params.wavelength)
    if params.direction is not None:
      beam.set_direction(params.direction)

  def override_detector(self, detector, params, beam, beam_params):
    '''
    Override the detector parameters

    '''

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

      # Update the parallax correction
      if (panel_params.material is not None or
          panel_params.thickness is not None or
          beam_params.wavelength is not None):
        from dxtbx.model import ParallaxCorrectedPxMmStrategy
        from cctbx.eltbx import attenuation_coefficient
        table = attenuation_coefficient.get_table(panel.get_material())
        mu = table.mu_at_angstrom(beam.get_wavelength()) / 10.0
        t0 = panel.get_thickness()
        px_mm = ParallaxCorrectedPxMmStrategy(mu, t0)
        panel.set_px_mm_strategy(px_mm)
        panel.set_mu(mu)



    for panel_id in frame_hash:
      fast_axis = frame_hash[panel_id]['fast_axis']
      slow_axis = frame_hash[panel_id]['slow_axis']
      origin = frame_hash[panel_id]['origin']
      if [fast_axis, slow_axis, origin].count(None) == 0:
        panel = detector[panel_id]
        panel.set_frame(fast_axis, slow_axis, origin)
      elif [fast_axis, slow_axis, origin].count(None) != 3:
        raise Sorry("fast_axis, slow_axis and origin must be set together")

  def override_goniometer(self, goniometer, params):
    '''
    Override the goniometer parameters

    '''
    if goniometer is not None:
      if len(params.axis) == 1:
        goniometer.set_rotation_axis(params.axis[0])
      elif len(params.axis) > 1:
        if not hasattr(goniometer, 'get_axes'):
          raise Sorry("Current goniometer is not a multi-axis goniometer")
        if len(goniometer.get_axes()) != len(params.axis):
          raise Sorry("Number of axes must match the current goniometer (%s)"
             %len(goniometer.get_axes()))
        from scitbx.array_family import flex
        goniometer.set_axes(flex.vec3_double(params.axis))
      if len(params.angle) > 1:
        if not hasattr(goniometer, 'get_axes'):
          raise Sorry("Current goniometer is not a multi-axis goniometer")
        if len(goniometer.get_angles()) != len(params.angle):
          raise Sorry("Number of angles must match the current goniometer (%s)"
             %len(goniometer.get_angles()))
        goniometer.set_angles(params.angle)
      if params.fixed_rotation is not None:
        goniometer.set_fixed_rotation(params.fixed_rotation)
      if params.setting_rotation is not None:
        goniometer.set_setting_rotation(params.setting_rotation)

  def override_scan(self, scan, params):
    '''
    Override the scan parameters

    '''
    if scan is not None:
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


class MetaDataUpdater(object):
  '''
  A class to manage updating the datablock metadata

  '''
  def __init__(self, params):
    '''
    Init the class

    '''
    self.params = params

    # Create the geometry updater
    if self.params.input.reference_geometry is not None:
      self.update_geometry = ReferenceGeometryUpdater(self.params)
    elif self.params.geometry.mosflm_beam_centre is not None:
      self.update_geometry = MosflmBeamCenterUpdater(self.params)
    else:
      self.update_geometry = ManualGeometryUpdater(self.params)

  def __call__(self, datablock):
    '''
    Transform the metadata

    '''
    from dxtbx.datablock import DataBlock

    # Import the lookup data
    lookup = self.import_lookup_data(self.params)

    # Convert all to ImageGrid
    if self.params.input.as_grid_scan:
      datablock = self.convert_to_grid_scan(datablock, self.params)

    # Init imageset list
    imageset_list = []

    # Loop through imagesets
    for imageset in datablock.extract_imagesets():

      # Check beam and detector are present
      if imageset.get_beam() == None or imageset.get_detector() == None:
        raise Sorry('''
          Imageset contains no beam or detector model. This means you will be
          unable to process your data.

          Possible causes of this error are:
             - A problem reading the images with one of the dxtbx format classes
             - A lack of header information in the file itself.
        ''')

      # Set the external lookups
      imageset = self.update_lookup(imageset, lookup)

      # Update the geometry
      imageset = self.update_geometry(imageset)

      # Append to new imageset list
      imageset_list.append(imageset)

    # Return the datablock
    return DataBlock(imageset_list)

  def update_lookup(self, imageset, lookup):
    if lookup.size is not None:
      d = imageset.get_detector()
      assert len(lookup.size) == len(d), "Incompatible size"
      for s, p in zip(lookup.size, d):
        assert s == p.get_image_size()[::-1], "Incompatible size"
      if lookup.mask.filename is not None:
        imageset.external_lookup.mask = lookup.mask
      if lookup.gain.filename is not None:
        imageset.external_lookup.gain = lookup.gain
      if lookup.dark.filename is not None:
        imageset.external_lookup.pedestal = lookup.dark
    return imageset

  def import_lookup_data(self, params):
    '''
    Get the lookup data

    '''
    from collections import namedtuple
    from dxtbx.imageset import ExternalLookupItem
    import cPickle as pickle
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
    Lookup = namedtuple("Lookup", ['size', 'mask', 'gain', 'dark'])
    return Lookup(
      size=lookup_size,
      mask=ExternalLookupItem(data=mask, filename=mask_filename),
      gain=ExternalLookupItem(data=gain, filename=gain_filename),
      dark=ExternalLookupItem(data=dark, filename=dark_filename))

  def convert_to_grid_scan(self, datablock, params):
    '''
    Convert the imagesets to grid scans

    '''
    from dxtbx.datablock import DataBlock
    from dxtbx.imageset import ImageGrid
    if params.input.grid_size is None:
      raise Sorry("The input.grid_size parameter is required")
    sweeps = datablock.extract_sweeps()
    stills = datablock.extract_stills()
    imagesets = []
    for iset in sweeps + stills:
      imagesets.append(ImageGrid.from_imageset(iset, params.input.grid_size))
    return DataBlock(imagesets)


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
    from dials.util import log

    # Parse the command line arguments in two passes to set up logging early
    params, options = self.parser.parse_args(show_diff_phil=False, quick_parse=True)

    # Configure logging
    log.config(
      params.verbosity,
      info=params.output.log,
      debug=params.output.debug_log)
    from dials.util.version import dials_version
    logger.info(dials_version())

    # Parse the command line arguments completely
    params, options = self.parser.parse_args(show_diff_phil=False)

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      logger.info('The following parameters have been modified:\n')
      logger.info(diff_phil)

    # Print help if no input
    if len(params.input.datablock) == 0 and not params.input.template:
      self.parser.print_help()
      return

    # Setup the datablock importer
    datablock_importer = DataBlockImporter(params)

    # Setup the metadata updater
    metadata_updater = MetaDataUpdater(params)

    # Get the datablocks
    datablock = metadata_updater(datablock_importer())

    # Extract any sweeps
    sweeps = datablock.extract_sweeps()

    # Extract any stills
    stills = datablock.extract_stills()
    if not stills:
      num_stills = 0
    else:
      num_stills = sum([len(s) for s in stills])

    # Print some data block info - override the output of image range
    # if appropriate
    image_range = params.geometry.scan.image_range

    logger.info("-" * 80)
    logger.info("  format: %s" % str(datablock.format_class()))
    if image_range is None:
      logger.info("  num images: %d" % datablock.num_images())
    else:
      logger.info("  num images: %d" % (image_range[1] - image_range[0] + 1))
    logger.info("  num sweeps: %d" % len(sweeps))
    logger.info("  num stills: %d" % num_stills)

    # Loop through all the sweeps
    for j, sweep in enumerate(sweeps):
      logger.debug("")
      logger.debug("Sweep %d" % j)
      logger.debug("  Length %d" % len(sweep))
      logger.debug(sweep.get_beam())
      logger.debug(sweep.get_goniometer())
      logger.debug(sweep.get_detector())
      logger.debug(sweep.get_scan())

    # Only allow a single sweep
    if params.input.allow_multiple_sweeps is False:
      self.assert_single_sweep(sweeps, params)

    # Write the datablocks to file
    self.write_datablocks([datablock], params)

  def write_datablocks(self, datablocks, params):
    '''
    Output the datablock to file.

    '''
    if params.output.datablock:
      logger.info("-" * 80)
      logger.info('Writing datablocks to %s' % params.output.datablock)
      dump = DataBlockDumper(datablocks)
      dump.as_file(params.output.datablock, compact=params.output.compact)

  def assert_single_sweep(self, sweeps, params):
    '''
    Print an error message if more than 1 sweep
    '''
    if len(sweeps) > 1:

      # Print some info about multiple sweeps
      self.diagnose_multiple_sweeps(sweeps, params)

      # Raise exception
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

  def diagnose_multiple_sweeps(self, sweeps, params):
    '''
    Print a diff between sweeps.

    '''
    logger.info("")
    for i in range(1, len(sweeps)):
      logger.info("=" * 80)
      logger.info("Diff between sweep %d and %d" % (i-1, i))
      logger.info("")
      self.print_sweep_diff(sweeps[i-1], sweeps[i], params)
    logger.info("=" * 80)
    logger.info("")

  def print_sweep_diff(self, sweep1, sweep2, params):
    '''
    Print a diff between sweeps.

    '''
    from dxtbx.datablock import SweepDiff
    diff = SweepDiff(params.input.tolerance)
    diff(sweep1, sweep2)



if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
