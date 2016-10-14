#!/usr/bin/env python
#
# export.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
from __future__ import division
from libtbx.phil import parse

help_message = '''

This program is used to export the results of dials processing in various
formats.

The output formats currently supported are:

MTZ format exports the files as an unmerged mtz file, ready for input to
downstream programs such as Pointless and Aimless. The required input is an
experiments.json file and an integrated.pickle file.

NXS format exports the files as an NXmx file. The required input is an
experiments.json file and an integrated.pickle file.

MOSFLM format exports the files as an index.mat mosflm-format matrix file and a
mosflm.in file containing basic instructions for input to mosflm. The required
input is an experiments.json file.

XDS format exports an experiments.json file as XDS.INP and XPARM.XDS files. If a
reflection pickle is given it will be exported as a SPOT.XDS file.

Examples::

  # Export to mtz
  dials.export experiments.json integrated.pickle
  dials.export experiments.json integrated.pickle mtz.hklout=integrated.mtz

  # Export to nexus
  dials.export experiments.json integrated.pickle format=nxs
  dials.export experiments.json integrated.pickle format=nxs nxs.hklout=integrated.nxs

  # Export to mosflm
  dials.export experiments.json integrated.pickle format=mosflm

  # Export to xds
  dials.export strong.pickle format=xds
  dials.export indexed.pickle format=xds
  dials.export experiments.json format=xds
  dials.export experiments.json indexed.pickle format=xds

'''

phil_scope = parse('''

  format = *mtz sadabs nxs mosflm xds best xds_ascii
    .type = choice
    .help = "The output file format"

  summation = False
    .type = bool
    .help = "Output summation integrated data (default profile fitted)"

  debug = False
    .type = bool
    .help = "Output additional debugging information"

  mtz {

    ignore_panels = False
      .type = bool
      .help = "Ignore multiple panels / detectors in output"

    include_partials = False
      .type = bool
      .help = "Include partial reflections (scaled) in output"

    keep_partials = False
      .type = bool
      .help = "Keep low partiality reflections"

    min_isigi = -5
      .type = float
      .help = "Exclude reflections with unfeasible values of I/Sig(I)"

    force_static_model = False
      .type = bool
      .help = "Force program to use static model even if scan varying is present"

    hklout = integrated.mtz
      .type = path
      .help = "The output MTZ file"
  }

  sadabs {

    hklout = integrated.sad
      .type = path
      .help = "The output raw sadabs file"
    run = 0
      .type = int
      .help = "Batch number / run number for output file"
    predict = False
      .type = bool
      .help = "Compute centroids with static model, not observations"

  }

  xds_ascii {

    hklout = DIALS.HKL
      .type = path
      .help = "The output raw hkl file"

  }

  nxs {

    hklout = integrated.nxs
      .type = path
      .help = "The output Nexus file"

  }

  mosflm {

    directory = mosflm
      .type = path
      .help = "The output directory for mosflm output"

  }

  xds {

    directory = xds
      .type = path
      .help = "The output directory for xds output"

  }

  best {

    prefix = best
      .type = str
      .help = "The prefix for the output file names for best"
              "(.hkl, .dat and .par files)"

    n_bins = 100
      .type = int(value_min=1)
      .help = "Number of resolution bins for background estimation"

  }

  output {

    log = dials.export.log
      .type = path
      .help = "The log filename"

    debug_log = dials.export.debug.log
      .type = path
      .help = "The debug log filename"

  }
''')


class MTZExporter(object):
  '''
  A class to export stuff in MTZ format

  '''

  def __init__(self, params, experiments, reflections):
    '''
    Initialise the exporter

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables

    '''

    # Check the input
    if len(experiments) == 0:
      raise Sorry('MTZ exporter requires an experiment list')
    if len(reflections) != 1:
      raise Sorry('MTZ exporter requires 1 reflection table')

    # Save the stuff
    self.params = params
    self.experiments = experiments
    self.reflections = reflections[0]

  def export(self):
    '''
    Export the files

    '''
    from dials.util.export_mtz import export_mtz
    m = export_mtz(
      self.reflections,
      self.experiments,
      self.params.mtz.hklout,
      ignore_panels=params.mtz.ignore_panels,
      include_partials=params.mtz.include_partials,
      keep_partials=params.mtz.keep_partials,
      min_isigi=params.mtz.min_isigi,
      force_static_model=params.mtz.force_static_model)
    from cStringIO import StringIO
    summary = StringIO()
    m.show_summary(out=summary)
    info(summary.getvalue())

class SadabsExporter(object):
  '''
  A class to export stuff in HKL format

  '''

  def __init__(self, params, experiments, reflections):
    '''
    Initialise the exporter

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables

    '''

    # Check the input
    if len(experiments) == 0:
      raise Sorry('SADABS exporter requires an experiment list')
    if len(reflections) != 1:
      raise Sorry('SADABS exporter requires 1 reflection table')

    # Save the stuff
    self.params = params
    self.experiments = experiments
    self.reflections = reflections[0]

  def export(self):
    from dials.util.export_sadabs import export_sadabs
    export_sadabs(
      self.reflections,
      self.experiments,
      self.params.sadabs.hklout,
      run=self.params.sadabs.run,
      summation=self.params.summation,
      include_partials=params.mtz.include_partials,
      keep_partials=params.mtz.keep_partials,
      debug=params.debug,
      predict=params.sadabs.predict)

class XDSASCIIExporter(object):
  '''
  A class to export stuff in XDS_ASCII format

  '''

  def __init__(self, params, experiments, reflections):
    '''
    Initialise the exporter

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables

    '''

    # Check the input
    if len(experiments) == 0:
      raise Sorry('XDS_ASCII exporter requires an experiment list')
    if len(reflections) != 1:
      raise Sorry('XDS_ASCII exporter requires 1 reflection table')

    # Save the stuff
    self.params = params
    self.experiments = experiments
    self.reflections = reflections[0]

  def export(self):
    from dials.util.export_xds_ascii import export_xds_ascii
    export_xds_ascii(
      self.reflections,
      self.experiments,
      self.params.xds_ascii.hklout,
      summation=self.params.summation,
      include_partials=params.mtz.include_partials,
      keep_partials=params.mtz.keep_partials)

class NexusExporter(object):
  '''
  A class to export stuff in Nexus format

  '''

  def __init__(self, params, experiments, reflections):
    '''
    Initialise the exporter

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables

    '''

    # Check the input
    if len(experiments) == 0:
      raise Sorry('Nexus exporter requires an experiment list')
    if len(reflections) != 1:
      raise Sorry('Nexus exporter requires 1 reflection table')

    # Save the stuff
    self.params = params
    self.experiments = experiments
    self.reflections = reflections[0]

  def export(self):
    '''
    Export the files

    '''
    from dials.util.nexus import dump
    dump(
      self.experiments,
      self.reflections,
      self.params.nxs.hklout)


class MosflmExporter(object):
  '''
  A class to export stuff in mosflm format

  '''

  def __init__(self, params, experiments, reflections):
    '''
    Initialise the exporter

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables

    '''

    # Check the input
    if len(experiments) == 0:
      raise Sorry('Mosflm exporter requires an experiment list')
    if len(reflections) != 0:
      raise Sorry('Mosflm exporter does need reflection table')

    # Save the stuff
    self.params = params
    self.experiments = experiments

  def export(self):
    '''
    Export the files

    '''
    from dials.util.mosflm import dump
    dump(
      self.experiments,
      self.params.mosflm.directory)


class XDSExporter(object):
  '''
  A class to export stuff in mosflm format

  '''

  def __init__(self, params, experiments, reflections):
    '''
    Initialise the exporter

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables

    '''

    # Check the input
    if len(reflections) > 1:
      raise Sorry('XDS exporter requires 0 or 1 reflection table')

    # Save the stuff
    self.params = params
    self.experiments = experiments
    if len(reflections) == 0:
      self.reflections = reflections
    else:
      self.reflections = reflections[0]

  def export(self):
    '''
    Export the files

    '''
    from dials.util.xds import dump
    dump(
      self.experiments,
      self.reflections,
      self.params.xds.directory)


class BestExporter(object):
  '''
  A class to export stuff in BEST format

  '''

  def __init__(self, params, experiments, reflections):
    '''
    Initialise the exporter

    :param params: The phil parameters
    :param experiments: The experiment list
    :param reflections: The reflection tables

    '''

    # Check the input
    if len(experiments) == 0:
      raise Sorry('BEST exporter requires an experiment list')
    if len(reflections) == 0:
      raise Sorry('BEST exporter require a reflection table')

    # Save the stuff
    self.params = params
    self.experiments = experiments
    self.reflections = reflections

  def export(self):
    '''
    Export the files

    '''
    from dials.util import best

    experiment = self.experiments[0]
    reflections = self.reflections[0]
    imageset = experiment.imageset
    prefix = self.params.best.prefix

    best.write_background_file(
      '%s.dat' %prefix, imageset, n_bins=self.params.best.n_bins)
    best.write_integrated_hkl(prefix, reflections)
    best.write_par_file('%s.par' %prefix, experiment)


if __name__ == '__main__':
  import libtbx.load_env
  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  from dials.util.version import dials_version
  from dials.util import log
  from libtbx.utils import Sorry
  from logging import info

  usage = '%s experiments.json reflections.pickle [options]' % (
              libtbx.env.dispatcher_name)

  # Create the option parser
  parser = OptionParser(
    usage = usage,
    read_experiments=True,
    read_reflections=True,
    phil=phil_scope,
    epilog=help_message)

  # Get the parameters
  params, options = parser.parse_args(show_diff_phil=False)

  # Configure the logging
  log.config(
    info=params.output.log,
    debug=params.output.debug_log)

  # Print the version number
  info(dials_version())

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    info('The following parameters have been modified:\n')
    info(diff_phil)

  # Get the experiments and reflections
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)
  if len(reflections) == 0 and len(experiments) == 0:
    parser.print_help()
    exit(0)

  # Choose the exporter
  if params.format == 'mtz':
    exporter = MTZExporter(params, experiments, reflections)
  elif params.format == 'sadabs':
    exporter = SadabsExporter(params, experiments, reflections)
  elif params.format == 'xds_ascii':
    exporter = XDSASCIIExporter(params, experiments, reflections)
  elif params.format == 'nxs':
    exporter = NexusExporter(params, experiments, reflections)
  elif params.format == 'mosflm':
    exporter = MosflmExporter(params, experiments, reflections)
  elif params.format == 'xds':
    exporter = XDSExporter(params, experiments, reflections)
  elif params.format == 'best':
    exporter = BestExporter(params, experiments, reflections)
  else:
    raise Sorry('Unknown format: %s' % params.format)

  # Export the data
  exporter.export()
