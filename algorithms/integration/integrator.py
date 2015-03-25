#
# integrator.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials_algorithms_integration_integrator_ext import *
from dials.algorithms.integration.processor import Processor3D
from dials.algorithms.integration.processor import ProcessorFlat3D
from dials.algorithms.integration.processor import Processor2D
from dials.algorithms.integration.processor import ProcessorSingle2D
from dials.algorithms.integration.processor import ProcessorStills
from dials.algorithms.integration.processor import job
from iotbx import phil
import libtbx


def generate_phil_scope():
  '''
  Generate the integration phil scope.

  :return: The phil scope

  '''
  import dials.extensions
  from dials.interfaces import BackgroundIface
  from dials.interfaces import CentroidIface

  phil_scope = phil.parse('''

    integration {

      include scope dials.algorithms.integration.processor.phil_scope

      integrator = *auto 3d flat3d 2d single2d stills
        .type = choice
        .help = "The integrator to use."
        .expert_level=3

      profile_fitting = True
        .type = bool
        .help = "Use profile fitting if available"

      filter
        .expert_level = 1
      {
        min_zeta = 0.05
          .help = "Filter the reflections by the value of zeta. A value of less"
                  "than or equal to zero indicates that this will not be used. A"
                  "positive value is used as the minimum permissable value."
          .type = float(value_min=0.0, value_max=1.0)

        max_shoebox_overlap = 1.0
          .type = float(value_min=0.0, value_max=1.0)
          .help = "Filter reflections whose shoeboxes are overlapped by greater"
                  "than the requested amount. Note that this is not the"
                  "percentage of the peak that is overlapped but rather the"
                  "percentage of the shoebox (background and foreground). This"
                  "can be useful when the detector is too close and many"
                  "overlapping reflections are predicted at high resolution"
                  "causing memory issues."

        include scope dials.algorithms.integration.filtering.phil_scope
      }
    }
  ''', process_includes=True)
  main_scope = phil_scope.get_without_substitution("integration")
  assert len(main_scope) == 1
  main_scope = main_scope[0]
  main_scope.adopt_scope(BackgroundIface.phil_scope())
  main_scope.adopt_scope(CentroidIface.phil_scope())
  return phil_scope

# The integration phil scope
phil_scope = generate_phil_scope()

def hist(data, width=80, symbol='#', prefix=''):
  '''
  A utility function to print a histogram of reflections on frames.

  :param data: The data to histogram
  :param width: The number of characters in each line
  :param symbol: The plot symbol
  :param prefix: String to prefix to each line
  :return: The histogram string

  '''
  from collections import defaultdict, Counter
  from math import log10, floor
  assert len(data) > 0, "Need > 0 reflections"
  assert width > 0, "Width should be > 0"
  count = Counter(data)
  count = count.items()
  count.sort()
  frame, count = zip(*count)
  min_frame = min(frame)
  max_frame = max(frame)
  min_count = min(count)
  max_count = max(count)
  assert max_count > 0, "Max should be > 0"
  assert min_count >= 0, "Min should be >= 0"
  if max_frame == 0:
    num_frame_zeros = 1
  else:
    num_frame_zeros = int(floor(log10(max_frame))) + 1
  num_count_zeros = int(floor(log10(max_count))) + 1
  assert num_frame_zeros > 0, "Num should be > 0"
  assert num_count_zeros > 0, "Num should be > 0"
  num_hist = width - (num_frame_zeros + num_count_zeros + 5) - len(prefix)
  assert num_hist > 0, "num_hist should be > 0"
  fmt = '%s%%-%dd [%%-%dd]: %%s' % (prefix, num_frame_zeros, num_count_zeros)
  scale = float(num_hist) / max_count
  return '\n'.join((
    fmt % (key, value, int(value * scale) * symbol)
      for key, value in zip(frame, count)))

def frame_hist(bbox, width=80, symbol='#', prefix=''):
  '''
  A utility function to print a histogram of reflections on frames.

  :param bbox: The bounding boxes
  :param width: The width of each line
  :param symbol: The histogram symbol
  :param prefix: A string to prefix to each line
  :return: The histogram string

  '''
  return hist(
    [z for b in bbox for z in range(b[4], b[5])],
    width=width,
    symbol=symbol,
    prefix=prefix)

def nframes_hist(bbox, width=80, symbol='#', prefix=''):
  '''
  A utility function to print a histogram of number of frames.

  :param bbox: The bounding boxes
  :param width: The width of each line
  :param symbol: The histogram symbol
  :param prefix: A string to prefix to each line
  :return: The histogram string

  '''
  return hist(
    [b[5] - b[4] for b in bbox],
    width=width,
    symbol=symbol,
    prefix=prefix)


class InitializerRot(object):
  '''
  A pre-processing class for oscillation data.

  '''

  def __init__(self,
               experiments,
               profile_model,
               params):
    '''
    Initialise the pre-processor.

    '''
    self.experiments = experiments
    self.profile_model = profile_model
    self.params = params

  def __call__(self, reflections):
    '''
    Do some pre-processing.

    '''
    from dials.algorithms.integration.filtering import MultiPowderRingFilter
    from dials.array_family import flex
    from scitbx.array_family import shared

    # Compute some reflection properties
    reflections.compute_zeta_multi(self.experiments)
    reflections.compute_d(self.experiments)
    reflections.compute_bbox(self.experiments, self.profile_model)

    # Filter the reflections by zeta
    mask = flex.abs(reflections['zeta']) < self.params.integration.filter.min_zeta
    num_ignore = mask.count(True)
    reflections.set_flags(mask, reflections.flags.dont_integrate)

    # Filter the reflections by powder ring
    powder_filter = MultiPowderRingFilter.from_params(
      self.params.integration.filter)
    if powder_filter is not None:
      mask = powder_filter(reflections['d'])
      reflections.set_flags(mask, reflections.flags.in_powder_ring)


class InitializerStills(object):
  '''
  A pre-processing class for stills data.

  '''

  def __init__(self,
               experiments,
               profile_model,
               params):
    '''
    Initialise the pre-processor.

    '''
    self.experiments = experiments
    self.profile_model = profile_model
    self.params = params

  def __call__(self, reflections):
    '''
    Do some pre-processing.

    '''
    from dials.algorithms.integration.filtering import MultiPowderRingFilter
    from dials.array_family import flex

    # Compute some reflection properties
    reflections.compute_d(self.experiments)
    reflections.compute_bbox(self.experiments, self.profile_model)

    # Check the bounding boxes are all 1 frame in width
    z0, z1 = reflections['bbox'].parts()[4:6]
    assert (z1 - z0).all_eq(1), "bbox is invalid"

    # Filter the reflections by powder ring
    powder_filter = MultiPowderRingFilter.from_params(
      self.params.integration.filter)
    if powder_filter is not None:
      mask = powder_filter(reflections['d'])
      reflections.set_flags(mask, reflections.flags.in_powder_ring)


class FinalizerRot(object):
  '''
  A post-processing class for oscillation data.

  '''

  def __init__(self, experiments, profile_model, params):
    '''
    Initialise the post processor.

    '''
    self.experiments = experiments
    self.profile_model = profile_model
    self.params = params

  def __call__(self, reflections):
    '''
    Do some post processing.

    '''

    # Compute the corrections
    reflections.compute_corrections(self.experiments)


class FinalizerStills(object):
  '''
  A post-processing class for stills data.

  '''

  def __init__(self, experiments, profile_model, params):
    '''
    Initialise the post processor.

    '''
    self.experiments = experiments
    self.profile_model = profile_model
    self.params = params

  def __call__(self, reflections):
    '''
    Do some post processing.

    '''
    pass


class ProfileModellerExecutor(Executor):
  '''
  The class to do profile modelling calculations

  '''

  def __init__(self, experiments, profile_model):
    '''
    Initialise the executor

    :param experiments: The experiment list
    :param profile_model: The profile model

    '''
    self.experiments = experiments
    self.profile_model = profile_model
    super(ProfileModellerExecutor, self).__init__()

  def initialize(self, frame0, frame1, reflections):
    '''
    Initialise the processing for a job

    :param frame0: The first frame in the job
    :param frame1: The last frame in the job
    :param reflections: The reflections that will be processed

    '''
    from logging import info

    # Get some info
    EPS = 1e-7
    full_value = (0.997300203937 - EPS)
    fully_recorded = reflections['partiality'] > full_value
    npart = fully_recorded.count(False)
    nfull = fully_recorded.count(True)
    nice = reflections.get_flags(reflections.flags.in_powder_ring).count(True)
    ntot = len(reflections)

    # Write some output
    info(" Beginning modelling job %d" % job.index)
    info("")
    info(" Frames: %d -> %d" % (frame0, frame1))
    info("")
    info(" Number of reflections")
    info("  Partial:     %d" % npart)
    info("  Full:        %d" % nfull)
    info("  In ice ring: %d" % nice)
    info("  Total:       %d" % ntot)
    info("")

    # Print a histogram of reflections on frames
    if frame1 - frame0 > 1:
      info(' The following histogram shows the number of reflections predicted')
      info(' to have all or part of their intensity on each frame.')
      info('')
      info(frame_hist(reflections['bbox'], prefix=' ', symbol='*'))
      info('')

    # Create the modeller
    self.modeller = self.profile_model.modeller(self.experiments)

  def process(self, frame, reflections):
    '''
    Process the data

    :param frame: The frame being processed
    :param reflections: The reflections to process

    '''
    from logging import info

    # Check if pixels are overloaded
    reflections.is_overloaded(self.experiments)

    # Compute the shoebox mask
    reflections.compute_mask(self.experiments, self.profile_model)

    # Process the data
    reflections.compute_background(self.experiments)
    reflections.compute_centroid(self.experiments)
    reflections.compute_summed_intensity()

    # Do the profile modelling
    self.modeller.model(reflections)

    # Print some info
    fmt = ' Modelled % 5d / % 5d reflection profiles on image %d'
    nmod = reflections.get_flags(reflections.flags.used_in_modelling).count(True)
    ntot = len(reflections)
    info(fmt % (nmod, ntot, frame))

  def finalize(self):
    '''
    Finalize the processing

    '''
    pass

  def data(self):
    '''
    :return: the modeller

    '''
    return self.modeller


class IntegratorExecutor(Executor):
  '''
  The class to process the integration data

  '''

  def __init__(self, experiments, profile_model):
    '''
    Initialize the executor

    :param experiments: The experiment list
    :param profile_model: The profile model

    '''
    self.experiments = experiments
    self.profile_model = profile_model
    self.overlaps = None
    super(IntegratorExecutor, self).__init__()

  def initialize(self, frame0, frame1, reflections):
    '''
    Initialize the processing for the job

    :param frame0: The first frame to process
    :param frame1: The last frame to process
    :param reflections: The reflections to process

    '''
    from logging import info

    # Get some info
    EPS = 1e-7
    full_value = (0.997300203937 - EPS)
    fully_recorded = reflections['partiality'] > full_value
    npart = fully_recorded.count(False)
    nfull = fully_recorded.count(True)
    nice = reflections.get_flags(reflections.flags.in_powder_ring).count(True)
    nint = reflections.get_flags(reflections.flags.dont_integrate).count(False)
    ntot = len(reflections)

    # Write some output
    info(" Beginning integration job %d" % job.index)
    info("")
    info(" Frames: %d -> %d" % (frame0, frame1))
    info("")
    info(" Number of reflections")
    info("  Partial:     %d" % npart)
    info("  Full:        %d" % nfull)
    info("  In ice ring: %d" % nice)
    info("  Integrate:   %d" % nint)
    info("  Total:       %d" % ntot)
    info("")

    # Print a histogram of reflections on frames
    if frame1 - frame0 > 1:
      info(' The following histogram shows the number of reflections predicted')
      info(' to have all or part of their intensity on each frame.')
      info('')
      info(frame_hist(reflections['bbox'], prefix=' ', symbol='*'))
      info('')

    # Find any overlaps
    self.overlaps = reflections.find_overlaps(self.experiments)

  def process(self, frame, reflections):
    '''
    Process the reflections on a frame

    :param frame: The frame to process
    :param reflections: The reflections to process

    '''
    from logging import info
    from dials.algorithms.shoebox import MaskCode

    # Check if pixels are overloaded
    reflections.is_overloaded(self.experiments)

    # Compute the shoebox mask
    reflections.compute_mask(self.experiments, self.profile_model)

    # Process the data
    reflections.compute_background(self.experiments)
    reflections.compute_centroid(self.experiments)
    reflections.compute_summed_intensity()
    reflections.compute_fitted_intensity(self.experiments, self.profile_model)

    # Compute the number of background/foreground pixels
    sbox = reflections['shoebox']
    code1 = MaskCode.Valid
    code2 = MaskCode.Background | code1
    code3 = MaskCode.BackgroundUsed | code2
    code4 = MaskCode.Foreground | code1
    reflections['num_pixels.valid'] = sbox.count_mask_values(code1)
    reflections['num_pixels.background'] = sbox.count_mask_values(code2)
    reflections['num_pixels.background_used'] = sbox.count_mask_values(code3)
    reflections['num_pixels.foreground'] = sbox.count_mask_values(code4)

    # Print some info
    fmt = ' Integrated % 5d (sum) + % 5d (prf) /% 5d reflections on image %d'
    nsum = reflections.get_flags(reflections.flags.integrated_sum).count(True)
    nprf = reflections.get_flags(reflections.flags.integrated_prf).count(True)
    ntot = len(reflections)
    info(fmt % (nsum, nprf, ntot, frame))

  def finalize(self):
    '''
    Finalize the processing

    '''
    pass

  def data(self):
    '''
    Return data

    '''
    return None


class Integrator(object):
  '''
  The integrator class

  '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               params):
    '''
    Initialize the integrator

    :param experiments: The experiment list
    :param profile_model: The profile model
    :param reflections: The reflections to process
    :param params: The parameters to use

    '''

    # Save some stuff
    self.experiments = experiments
    self.profile_model = profile_model
    self.reflections = reflections
    self.params = params
    self.profile_model_report = None
    self.integration_report = None

  def integrate(self):
    '''
    Integrate the data

    '''
    from dials.algorithms.integration.report import IntegrationReport
    from dials.algorithms.integration.report import ProfileModelReport
    from dials.algorithms.integration.statistics import modeller_statistics
    from dials.util.command_line import heading
    from logging import info, debug
    from dials.util import pprint

    # Init the report
    self.profile_model_report = None
    self.integration_report = None

    # Heading
    info("=" * 80)
    info("")
    info(heading("Processing reflections"))
    info("")

    # Create summary format
    fmt = (
      ' Processing the following experiments:\n'
      '\n'
      ' Experiments: %d\n'
      ' Beams:       %d\n'
      ' Detectors:   %d\n'
      ' Goniometers: %d\n'
      ' Scans:       %d\n'
      ' Crystals:    %d\n'
      ' Imagesets:   %d\n'
    )

    # Print the summary
    info(fmt % (
      len(self.experiments),
      len(self.experiments.beams()),
      len(self.experiments.detectors()),
      len(self.experiments.goniometers()),
      len(self.experiments.scans()),
      len(self.experiments.crystals()),
      len(self.experiments.imagesets())))

    # Initialize the reflections
    initialize = self.InitializerClass(
      self.experiments,
      self.profile_model,
      self.params)
    initialize(self.reflections)

    # Check if we want to do some profile fitting
    if (self.params.integration.profile_fitting and
        self.profile_model.has_profile_fitting()):

      info("=" * 80)
      info("")
      info(heading("Modelling reflection profiles"))
      info("")

      # Get the selection
      selection = self.reflections.get_flags(
        self.reflections.flags.reference_spot)

      # Get the reference spots
      reference = self.reflections.select(selection)

      # Check if we need to skip
      if len(reference) == 0:
        info("** Skipping profile modelling - no reference profiles given **")

      # Create the data processor
      executor = ProfileModellerExecutor(
        self.experiments,
        self.profile_model)
      processor = self.ProcessorClass(
        self.experiments,
        self.profile_model,
        reference,
        self.params.integration)
      processor.executor = executor

      # Process the reference profiles
      reference, modeller_list, time_info = processor.process()

      # Set the reference spots info
      #self.reflections.set_selected(selection, reference)

      # Finalize the profile model
      assert len(modeller_list) > 0, "No modellers"
      modeller = None
      for index, mod in modeller_list.iteritems():
        if mod is None:
          continue
        if modeller is None:
          modeller = mod
        else:
          modeller.accumulate(mod)
      modeller.finalize()
      self.profile_model.profiles(modeller)

      # Print profiles
      for i in range(len(modeller)):
        m = modeller[i]
        debug("")
        debug("Profiles for experiment %d" % i)
        for j in range(len(m)):
          debug("Profile %d" % j)
          try:
            debug(pprint.profile3d(m.data(j)))
          except Exception:
            debug("** NO PROFILE **")

      # Print out some statistics
      # info("")
      # for summary in modeller_statistics(
      #     self.reflections,
      #     self.experiments,
      #     self.profile_model):
      #   info(str(summary))

      self.profile_model_report = ProfileModelReport(
        self.experiments,
        self.profile_model,
        self.reflections)
      info("")
      info(self.profile_model_report.as_str(prefix=' '))

      # Print the time info
      info("")
      info(str(time_info))
      info("")

    info("=" * 80)
    info("")
    info(heading("Integrating reflections"))
    info("")

    # Create the data processor
    executor = IntegratorExecutor(
      self.experiments,
      self.profile_model)
    processor = self.ProcessorClass(
      self.experiments,
      self.profile_model,
      self.reflections,
      self.params.integration)
    processor.executor = executor

    # Process the reflections
    self.reflections, _, time_info = processor.process()

    # Finalize the reflections
    finalize = self.FinalizerClass(
      self.experiments,
      self.profile_model,
      self.params)
    finalize(self.reflections)

    # Create the integration report
    self.integration_report = IntegrationReport(
      self.experiments,
      self.reflections)
    info("")
    info(self.integration_report.as_str(prefix=' '))

    # Print the time info
    info(str(time_info))
    info("")

    # Return the reflections
    return self.reflections

  def report(self):
    '''
    Return the report of the processing

    '''
    from collections import OrderedDict
    result = OrderedDict()
    if self.profile_model_report is not None:
      result['profile_model'] = self.profile_model_report.as_dict()
    if self.integration_report is not None:
      result['integration'] = self.integration_report.as_dict()
    return result

  def summary(self, block_size, block_size_units):
    ''' Print a summary of the integration stuff. '''
    from libtbx.table_utils import format as table
    from dials.util.command_line import heading
    from logging import info

    # Compute the task table
    if self._experiments.all_stills():
      rows = [["#", "Group", "Frame From", "Frame To"]]
      for i in range(len(self)):
        job = self._manager.job(i)
        group = job.index()
        f0, f1 = job.frames()
        rows.append([str(i), str(group), str(f0), str(f1)])
    elif self._experiments.all_sweeps():
      rows = [["#", "Group", "Frame From", "Frame To", "Angle From", "Angle To"]]
      for i in range(len(self)):
        job = self._manager.job(i)
        group = job.index()
        expr = job.expr()
        f0, f1 = job.frames()
        scan = self._experiments[expr[0]].scan
        p0 = scan.get_angle_from_array_index(f0)
        p1 = scan.get_angle_from_array_index(f1)
        rows.append([str(i), str(group), str(f0), str(f1), str(p0), str(p1)])
    else:
      raise RuntimeError('Experiments must be all sweeps or all stills')
    task_table = table(rows, has_header=True, justify="right", prefix=" ")


class Integrator3D(Integrator):
  '''
  Integrator for 3D algorithms

  '''
  InitializerClass = InitializerRot
  ProcessorClass = Processor3D
  FinalizerClass = FinalizerRot


class IntegratorFlat3D(Integrator):
  '''
  Integrator for flattened 3D algorithms

  '''
  InitializerClass = InitializerRot
  ProcessorClass = ProcessorFlat3D
  FinalizerClass = FinalizerRot


class Integrator2D(Integrator):
  '''
  Integrator for 2D algorithms

  '''
  InitializerClass = InitializerRot
  ProcessorClass = Processor2D
  FinalizerClass = FinalizerRot


class IntegratorSingle2D(Integrator):
  '''
  Integrator for 2D algorithms on a single image

  '''
  InitializerClass = InitializerRot
  ProcessorClass = ProcessorSingle2D
  FinalizerClass = FinalizerRot


class IntegratorStills(Integrator):
  '''
  Integrator for still algorithms

  '''
  InitializerClass = InitializerStills
  ProcessorClass = ProcessorStills
  FinalizerClass = FinalizerStills


class IntegratorFactory(object):
  '''
  A factory for creating integrators.

  '''

  @staticmethod
  def create(params, experiments, profile_model, reflections):
    '''
    Create the integrator from the input configuration.

    :param params: The input phil parameters
    :param experiments: The list of experiments
    :param profile_model: The profile model
    :param reflections: The reflections to integrate
    :return: The integrator class

    '''
    from dials.algorithms.integration.filtering import MultiPowderRingFilter
    from dials.interfaces import BackgroundIface
    from dials.interfaces import CentroidIface
    from dials.array_family import flex
    from libtbx.utils import Abort

    # Check the input
    if len(experiments) != len(profile_model):
      raise RuntimeError(
        'Number of experiments and profile models should be the same')

    # Check each experiment has an imageset
    for exp in experiments:
      if exp.imageset is None:
        raise Abort('''
          One or more experiment does not contain an imageset. Access to the
          image data is crucial for integration.
        ''')

    # Initialise the strategy classes
    BackgroundAlgorithm = BackgroundIface.extension(
      params.integration.background.algorithm)
    CentroidAlgorithm = CentroidIface.extension(
      params.integration.centroid.algorithm)

    # Set the algorithms in the reflection table
    flex.reflection_table._background_algorithm = \
      flex.strategy(BackgroundAlgorithm, params)
    flex.reflection_table._centroid_algorithm = \
      flex.strategy(CentroidAlgorithm, params)

    # Get the classes we need
    if params.integration.integrator == 'auto':
      if experiments.all_stills():
        params.integration.integrator = 'stills'
      else:
        params.integration.integrator = '3d'
    if params.integration.integrator == '3d':
      IntegratorClass = Integrator3D
    elif params.integration.integrator == 'flat3d':
      IntegratorClass = IntegratorFlat3D
    elif params.integration.integrator == '2d':
      IntegratorClass = Integrator2D
    elif params.integration.integrator == 'single2d':
      IntegratorClass = IntegratorSingle2D
    elif params.integration.integrator == 'stills':
      IntegratorClass = IntegratorStills
    else:
      raise RuntimeError("Unknown integration type")

    # Return an instantiation of the class
    return IntegratorClass(experiments, profile_model, reflections, params)
