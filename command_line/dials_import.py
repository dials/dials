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
from __future__ import absolute_import, division
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

    ignore_unhandled = False
      .type = bool
      .help = "Don't throw exception if some args are unhandled"

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

  include scope dials.util.options.format_phil_scope

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

      format_kwargs = {
        'dynamic_shadowing' : self.params.format.dynamic_shadowing
      }

      # Check if a template has been set and print help if not, otherwise try to
      # import the images based on the template input
      if len(self.params.input.template) > 0:
        importer = DataBlockTemplateImporter(
          self.params.input.template,
          max(self.params.verbosity-1, 0),
          format_kwargs=format_kwargs)
        datablocks = importer.datablocks
        if len(datablocks) == 0:
          raise Sorry('No datablocks found matching template %s' % self.params.input.template)
      elif len(self.params.input.directory) > 0:
        datablocks = DataBlockFactory.from_filenames(
          self.params.input.directory,
          max(self.params.verbosity-1, 0),
          format_kwargs=format_kwargs)
        if len(datablocks) == 0:
          raise Sorry('No datablocks found in directories %s' % self.params.input.directory)
      else:
        raise Sorry('No datablocks found')
    # if len(datablocks) > 1:
    #   raise Sorry("More than 1 datablock found")

    # Return the datablocks
    return datablocks


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
      static_only=True), "Reference detector model does not match input detector model"

    # Set beam and detector
    imageset.set_beam(self.reference.beam)
    imageset.set_detector(self.reference.detector)
    imageset.set_goniometer(self.reference.goniometer)
    return imageset

  def load_reference_geometry(self, params):
    '''
    Load a reference geometry file

    '''
    from collections import namedtuple
    # Load reference geometry
    reference_detector = None
    reference_beam = None
    if params.input.reference_geometry is not None:
      from dxtbx.serialize import load
      experiments, datablock = None, None
      try:
        experiments = load.experiment_list(
          params.input.reference_geometry, check_format=False)
      except Exception, e:
        datablock = load.datablock(params.input.reference_geometry)
      assert experiments or datablock, 'Could not import reference geometry'
      if experiments:
        assert len(experiments.detectors()) >= 1
        assert len(experiments.beams()) >= 1
        if len(experiments.detectors()) > 1:
          raise Sorry('The reference geometry file contains %d detector definitions, but only a single definition is allowed.' % len(experiments.detectors()))
        if len(experiments.beams()) > 1:
          raise Sorry('The reference geometry file contains %d beam definitions, but only a single definition is allowed.' % len(experiments.beams()))
        reference_detector = experiments.detectors()[0]
        reference_beam = experiments.beams()[0]
        reference_goniometer = experiments.goniometers()[0]
      else:
        assert len(datablock) == 1
        imageset = datablock[0].extract_imagesets()[0]
        reference_detector = imageset.get_detector()
        reference_beam = imageset.get_beam()
        reference_goniometer = imageset.get_goniometer()
    Reference = namedtuple("Reference", ["detector", "beam", "goniometer"])
    return Reference(
      detector=reference_detector,
      beam=reference_beam,
      goniometer=reference_goniometer)


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
    from dxtbx.imageset import ImageSet
    from dxtbx.imageset import ImageSweep
    from dxtbx.model import BeamFactory
    from dxtbx.model import DetectorFactory
    from dxtbx.model import GoniometerFactory
    from dxtbx.model import ScanFactory
    from copy import deepcopy
    if self.params.geometry.convert_sweeps_to_stills:
      imageset = ImageSet(data=imageset.data())
    if not isinstance(imageset, ImageSweep):
      if self.params.geometry.convert_stills_to_sweeps:
        imageset = self.convert_stills_to_sweep(imageset)
    if isinstance(imageset, ImageSweep):
      beam = BeamFactory.from_phil(
        self.params.geometry,
        imageset.get_beam())
      detector = DetectorFactory.from_phil(
        self.params.geometry,
        imageset.get_detector(),
        beam)
      goniometer = GoniometerFactory.from_phil(
        self.params.geometry,
        imageset.get_goniometer())
      scan = ScanFactory.from_phil(
        self.params.geometry,
        deepcopy(imageset.get_scan()))
      i0, i1 = scan.get_array_range()
      j0, j1 = imageset.get_scan().get_array_range()
      if i0 < j0 or i1 > j1:
        imageset = self.extrapolate_imageset(
          imageset   = imageset,
          beam       = beam,
          detector   = detector,
          goniometer = goniometer,
          scan       = scan)
      else:
        imageset.set_beam(beam)
        imageset.set_detector(detector)
        imageset.set_goniometer(goniometer)
        imageset.set_scan(scan)
    else:
      for i in range(len(imageset)):
        beam = BeamFactory.from_phil(
          self.params.geometry,
          imageset.get_beam(i))
        detector = DetectorFactory.from_phil(
          self.params.geometry,
          imageset.get_detector(i),
          beam)
        goniometer = GoniometerFactory.from_phil(
          self.params.geometry,
          imageset.get_goniometer(i))
        scan = ScanFactory.from_phil(
          self.params.geometry,
          imageset.get_scan(i))
        imageset.set_beam(beam, i)
        imageset.set_detector(detector, i)
        imageset.set_goniometer(goniometer, i)
        imageset.set_scan(scan, i)
    return imageset

  def extrapolate_imageset(self,
                           imageset=None,
                           beam = None,
                           detector = None,
                           goniometer = None,
                           scan = None):
    from dxtbx.imageset import ImageSetFactory
    first, last = scan.get_image_range()
    sweep = ImageSetFactory.make_sweep(
      template      = imageset.get_template(),
      indices       = range(first, last+1),
      format_class  = imageset.get_format_class(),
      beam          = beam,
      detector      = detector,
      goniometer    = goniometer,
      scan          = scan,
      format_kwargs = imageset.params())
    return sweep

  def convert_stills_to_sweep(self, imageset):
    from dxtbx.model import Scan
    assert self.params.geometry.scan.oscillation is not None
    beam = imageset.get_beam(index=0)
    detector = imageset.get_detector(index=0)
    goniometer = imageset.get_goniometer(index=0)
    for i in range(1, len(imageset)):
      b_i = imageset.get_beam(i)
      d_i = imageset.get_detector(i)
      g_i = imageset.get_goniometer(i)
      assert (beam is None and b_i is None) or beam.is_similar_to(
        imageset.get_beam(index=i),
        wavelength_tolerance            = self.params.input.tolerance.beam.wavelength,
        direction_tolerance             = self.params.input.tolerance.beam.direction,
        polarization_normal_tolerance   = self.params.input.tolerance.beam.polarization_normal,
        polarization_fraction_tolerance = self.params.input.tolerance.beam.polarization_fraction)
      assert (detector is None and d_i is None) or detector.is_similar_to(
        imageset.get_detector(index=i),
        fast_axis_tolerance = self.params.input.tolerance.detector.fast_axis,
        slow_axis_tolerance = self.params.input.tolerance.detector.slow_axis,
        origin_tolerance    = self.params.input.tolerance.detector.origin)
      assert (goniometer is None and g_i is None) or goniometer.is_similar_to(
        imageset.get_goniometer(index=i),
        rotation_axis_tolerance    = self.params.input.tolerance.goniometer.rotation_axis,
        fixed_rotation_tolerance   = self.params.input.tolerance.goniometer.fixed_rotation,
        setting_rotation_tolerance = self.params.input.tolerance.goniometer.setting_rotation)
    oscillation = self.params.geometry.scan.oscillation
    from dxtbx.sweep_filenames import template_regex_from_list
    from dxtbx.imageset import ImageSetFactory
    template, indices = template_regex_from_list(imageset.paths())
    image_range = (min(indices), max(indices))
    assert (image_range[1]+1 - image_range[0]) == len(indices)
    scan = Scan(image_range=image_range, oscillation=oscillation)
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


class MetaDataUpdater(object):
  '''
  A class to manage updating the datablock metadata

  '''
  def __init__(self, params):
    '''
    Init the class

    '''
    from dials.util.options import geometry_phil_scope

    self.params = params

    # Create the geometry updater
    self.update_geometry = []
    update_order = []

    # First add reference geometry is present
    if self.params.input.reference_geometry is not None:
      self.update_geometry.append(ReferenceGeometryUpdater(self.params))
      update_order.append("Reference geometry")

    # Then add manual geometry
    working_phil = geometry_phil_scope.format(self.params)
    diff_phil = geometry_phil_scope.fetch_diff(source=working_phil)
    if diff_phil.as_str() != "":
      self.update_geometry.append(ManualGeometryUpdater(self.params))
      update_order.append("Manual geometry")

    if len(update_order) > 0:
      logger.info("")
      logger.info("Applying input geometry in the following order:")
      for i, item in enumerate(update_order, start=1):
        logger.info("  %d. %s" % (i, item))
      logger.info("")

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

      # Set the external lookups
      imageset = self.update_lookup(imageset, lookup)

      # Update the geometry
      for updater in self.update_geometry:
        imageset = updater(imageset)

      # Check beam and detector are present
      if imageset.get_beam() == None or imageset.get_detector() == None:
        raise Sorry('''
          Imageset contains no beam or detector model. This means you will be
          unable to process your data.

          Possible causes of this error are:
             - A problem reading the images with one of the dxtbx format classes
             - A lack of header information in the file itself.

          You can override this by specifying the metadata as geometry parameters
        ''')

      # Append to new imageset list
      imageset_list.append(imageset)

    # Return the datablock
    return DataBlock(imageset_list)

  def update_lookup(self, imageset, lookup):
    from dxtbx.format.image import ImageBool, ImageDouble
    if lookup.size is not None:
      d = imageset.get_detector()
      assert len(lookup.size) == len(d), "Incompatible size"
      for s, p in zip(lookup.size, d):
        assert s == p.get_image_size()[::-1], "Incompatible size"
      if lookup.mask.filename is not None:
        imageset.external_lookup.mask.filename = lookup.mask.filename
        imageset.external_lookup.mask.data = ImageBool(lookup.mask.data)
      if lookup.gain.filename is not None:
        imageset.external_lookup.gain.filename = lookup.gain.filename
        imageset.external_lookup.gain.data = ImageDouble(lookup.gain.data)
      if lookup.dark.filename is not None:
        imageset.external_lookup.pedestal.filename = lookup.dark.filename
        imageset.external_lookup.pedestal.data = ImageDouble(lookup.dark.data)
    return imageset

  def import_lookup_data(self, params):
    '''
    Get the lookup data

    '''
    from collections import namedtuple
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
    Item = namedtuple("Item", ["data", "filename"])
    return Lookup(
      size=lookup_size,
      mask=Item(data=mask, filename=mask_filename),
      gain=Item(data=gain, filename=gain_filename),
      dark=Item(data=dark, filename=dark_filename))

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
    if params.input.ignore_unhandled:
      params, options, unhandled = self.parser.parse_args(
        show_diff_phil=False,
        return_unhandled=True)
    else:
      params, options = self.parser.parse_args(show_diff_phil=False)
      unhandled = None

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      logger.info('The following parameters have been modified:\n')
      logger.info(diff_phil)

    # Print a warning if something unhandled
    if unhandled is not None and len(unhandled) > 0:
      msg = 'Unable to handle the following arguments:\n'
      msg += '\n'.join(['  %s' % a for a in unhandled])
      msg += '\n'
      logger.warn(msg)

    # Print help if no input
    if (len(params.input.datablock) == 0 and not
        (params.input.template or params.input.directory)):
      self.parser.print_help()
      return

    # Setup the datablock importer
    datablock_importer = DataBlockImporter(params)

    # Setup the metadata updater
    metadata_updater = MetaDataUpdater(params)

    # Extract the datablocks and loop through
    datablocks = map(metadata_updater, datablock_importer())
    for datablock in datablocks:

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
    self.write_datablocks(datablocks, params)

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
