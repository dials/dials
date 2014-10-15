from __future__ import division
from iotbx import phil

class JobId(object):
  def __init__(self):
    self.value = None
  def __call__(self):
    return self.value

job_id = JobId()


def generate_phil_scope():
  ''' Generate the phil scope. '''
  import dials.extensions
  from dials.interfaces import BackgroundIface
  from dials.interfaces import IntensityIface
  from dials.interfaces import CentroidIface

  phil_scope = phil.parse('''

    integration {

      include scope dials.data.lookup.phil_scope
      include scope dials.data.multiprocessing.phil_scope

      block {
        size = 10
          .type = float
          .help = "The block size in rotation angle (degrees)."
      }

      filter
        .expert_level = 1
      {

        min_zeta = 0.05
          .help = "Filter the reflections by the value of zeta. A value of less"
                  "than or equal to zero indicates that this will not be used. A"
                  "positive value is used as the minimum permissable value."
          .type = float

        ice_rings {
          filter = False
            .type = bool
          unit_cell = 4.498,4.498,7.338,90,90,120
            .type = unit_cell
            .help = "The unit cell to generate d_spacings for ice rings."
          space_group = 194
            .type = space_group
            .help = "The space group used to generate d_spacings for ice rings."
          d_min = 0
            .type = int(value_min=0)
            .help = "The minimum resolution to filter ice rings"
          width = 0.06
            .type = float(value_min=0.0)
            .help = "The width of an ice ring (in d-spacing)."
        }
      }

      debug {
        save_shoeboxes = False
          .type = bool
          .help = "Save shoeboxes after each integration task."
      }
    }
  ''', process_includes=True)
  main_scope = phil_scope.get_without_substitution("integration")
  assert(len(main_scope) == 1)
  main_scope = main_scope[0]
  main_scope.adopt_scope(BackgroundIface.phil_scope())
  main_scope.adopt_scope(IntensityIface.phil_scope())
  main_scope.adopt_scope(CentroidIface.phil_scope())
  return phil_scope

# The integration phil scope
phil_scope = generate_phil_scope()


class Integrator(object):
  ''' Integrator interface class. '''

  def __init__(self, manager, nproc=1, mp_method='multiprocessing'):
    ''' Initialise the integrator.

    The integrator requires a manager class implementing the IntegratorManager
    interface. This class executes all the workers in separate threads and
    accumulates the results to expose to the user.

    Params:
      manager The integration manager
      nproc The number of processors
      mp_method The multiprocessing method

    '''
    self._manager = manager
    self._nproc = nproc
    self._mp_method = mp_method

  def integrate(self):
    ''' Do all the integration tasks.

    Returns
      The integration results

    '''
    from time import time
    from libtbx import easy_mp
    from libtbx.table_utils import format as table
    start_time = time()
    num_proc = len(self._manager)
    if self._nproc > 0:
      num_proc = min(num_proc, self._nproc)
    print ' Using %s with %d processes\n' % (self._mp_method, num_proc)
    if num_proc > 1:
      def process_output(result):
        print result[1]
        self._manager.accumulate(result[0])
      def execute_task(task):
        from cStringIO import StringIO
        import sys
        sys.stdout = StringIO()
        result = task()
        output = sys.stdout.getvalue()
        return result, output
      task_results = easy_mp.parallel_map(
        func=execute_task,
        iterable=list(self._manager.tasks()),
        processes=num_proc,
        callback=process_output,
        method=self._mp_method,
        preserve_order=True,
        preserve_exception_message=True)
      task_results, output = zip(*task_results)
    else:
      task_results = [task() for task in self._manager.tasks()]
      for result in task_results:
        self._manager.accumulate(result)
    assert(self._manager.finished())
    end_time = time()
    rows = [
      ["Read time"        , "%.2f seconds" % (self._manager.read_time)       ],
      ["Extract time"     , "%.2f seconds" % (self._manager.extract_time)    ],
      ["Pre-process time" , "%.2f seconds" % (self._manager.preprocess_time) ],
      ["Process time"     , "%.2f seconds" % (self._manager.process_time)    ],
      ["Post-process time", "%.2f seconds" % (self._manager.postprocess_time)],
      ["Total time"       , "%.2f seconds" % (self._manager.total_time)      ],
      ["User time"        , "%.2f seconds" % (end_time - start_time)         ],
    ]
    print ''
    print table(rows, justify='right', prefix=' ')
    return self._manager.result()


def frame_hist(bbox, width=80, symbol='#', prefix=''):
  ''' A utility function to print a histogram of reflections on frames. '''
  from collections import defaultdict, Counter
  from math import log10, floor
  assert(len(bbox) > 0)
  assert(width > 0)
  count = Counter([z for b in bbox for z in range(b[4], b[5])])
  count = count.items()
  count.sort()
  frame, count = zip(*count)
  min_frame = min(frame)
  max_frame = max(frame)
  min_count = min(count)
  max_count = max(count)
  assert(max_frame > min_frame)
  assert(max_count > min_count)
  assert(min_count >= 0)
  num_frame_zeros = int(floor(log10(max_frame))) + 1
  num_count_zeros = int(floor(log10(max_count))) + 1
  assert(num_frame_zeros > 0)
  assert(num_count_zeros > 0)
  num_hist = width - (num_frame_zeros + num_count_zeros + 5) - len(prefix)
  assert(num_hist) > 0
  fmt = '%s%%-%dd [%%-%dd]: %%s' % (prefix, num_frame_zeros, num_count_zeros)
  scale = float(num_hist) / max_count
  return '\n'.join((
    fmt % (key, value, int(value * scale) * symbol)
      for key, value in zip(frame, count)))


class Result(object):
  ''' A class representing an integration result. '''

  def __init__(self, index, reflections):
    ''' Initialise the data. '''
    self.index = index
    self.reflections = reflections


class Task(object):
  ''' A class to perform an integration task. '''

  def __init__(self,
               index,
               experiments,
               profile_model,
               data,
               job,
               flatten=False,
               save_shoeboxes=False):
    ''' Initialise the task. '''
    self._index = index
    self._experiments = experiments
    self._profile_model = profile_model
    self._data = data
    self._job = job
    self._flatten = flatten
    self._save_shoeboxes = save_shoeboxes
    self._mask = None

  def __call__(self):
    ''' Do the integration. '''
    from dials.array_family import flex
    from time import time
    from dials.util.command_line import heading

    # Set the global process ID
    job_id.value = self._index

    # Print out some info
    EPS = 1e-7
    fully_recorded = self._data['partiality'] > (1.0 - EPS)
    num_partial = fully_recorded.count(False)
    num_full = fully_recorded.count(True)
    num_integrate = self._data.get_flags(self._data.flags.dont_integrate).count(False)
    num_reference = self._data.get_flags(self._data.flags.reference_spot).count(True)
    num_ice_ring = self._data.get_flags(self._data.flags.in_powder_ring).count(True)
    print "=" * 80
    print ""
    print heading("Beginning integration job %d" % self._index)
    print ""
    print " Frames: %d -> %d" % self._job
    print ""
    print " Number of reflections"
    print "  Partial:     %d" % num_partial
    print "  Full:        %d" % num_full
    print "  In ice ring: %d" % num_ice_ring
    print "  Integrate:   %d" % num_integrate
    print "  Reference:   %d" % num_reference
    print "  Total:       %d" % len(self._data)
    print ""
    start_time = time()

    # Get the sub imageset
    imagesets = self._experiments.imagesets()
    assert(len(imagesets) == 1)
    imageset = imagesets[0]
    frame00, frame01 = self._job
    try:
      frame10, frame11 = imageset.get_array_range()
    except Exception:
      frame10, frame11 = (0, len(imageset))
    assert(frame00 < frame01)
    assert(frame10 < frame11)
    assert(frame00 >= frame10)
    assert(frame01 <= frame11)
    index0 = frame00 - frame10
    index1 = index0 + (frame01 - frame00)
    assert(index0 < index1)
    assert(index0 >= 0)
    assert(index1 <= len(imageset))
    imageset = imageset[index0:index1]

    # Print a histogram of reflections on frames
    if frame11 - frame10 > 1:
      print ' The number of reflections recorded on each frame'
      print ''
      print frame_hist(self._data['bbox'], prefix=' ', symbol='*')
      print ''

    # Create the shoeboxes
    self._data["shoebox"] = flex.shoebox(
      self._data["panel"],
      self._data["bbox"],
      allocate=True)

    # Extract the shoeboxes
    read_time, extract_time = self._data.extract_shoeboxes(imageset, self._mask)

    # Optionally flatten the shoeboxes
    if self._flatten:
      self._data['shoebox'].flatten()

    # Process the data
    process_start_time = time()
    self._process()
    process_time = time() - process_start_time

    # Optionally save the shoeboxes
    if self._save_shoeboxes:
      filename = 'shoeboxes_%d.pickle' % self._index
      print 'Saving shoeboxes to %s\n' % filename
      self._data.as_pickle(filename)

    # Delete the shoeboxes
    del self._data["shoebox"]

    # Return the result
    result = Result(self._index, self._data)
    result.read_time = read_time
    result.extract_time = extract_time
    result.process_time = process_time
    result.total_time = time() - start_time
    return result

  def _process(self):
    ''' Process the data. '''
    self._data.compute_mask(self._experiments, self._profile_model)
    self._data.integrate(self._experiments, self._profile_model)
    print ''


class Manager(object):
  ''' An class to manage integration book-keeping '''

  def __init__(self,
               preprocess,
               postprocess,
               experiments,
               profile_model,
               reflections,
               block_size=10,
               block_size_units='degrees',
               flatten=False,
               save_shoeboxes=False):
    ''' Initialise the manager. '''
    from dials.algorithms.integration import IntegrationManagerExecutor
    from dials.algorithms.integration import IntegrationJobCalculator
    imagesets = experiments.imagesets()
    detectors = experiments.detectors()
    scans = experiments.scans()
    assert(len(experiments) == len(profile_model))
    assert(len(imagesets) == 1)
    assert(len(detectors) == 1)
    imageset = imagesets[0]
    array_range = (0, len(imageset))
    if len(scans) != 0:
      assert(len(scans) == 1)
      scan = scans[0]
      if scan is not None:
        assert(len(imageset) == len(scan))
        array_range = scan.get_array_range()
    self._flatten = flatten
    self._save_shoeboxes = save_shoeboxes
    self._preprocess = preprocess
    self._postprocess = postprocess
    self._experiments = experiments
    self._profile_model = profile_model
    self._reflections = reflections
    if block_size_units == 'degrees':
      phi0, dphi = scan.get_oscillation()
      block_size_frames = block_size / dphi
    elif block_size_units == 'frames':
      block_size_frames = block_size
    else:
      raise RuntimeError('Unknown block_size_units = %s' % block_size_units)
    jobcalculator = IntegrationJobCalculator(
      array_range,
      block_size_frames)
    self._preprocess(self._reflections, jobcalculator.jobs())
    self._manager = IntegrationManagerExecutor(
      jobcalculator,
      self._reflections)
    self.read_time = 0
    self.extract_time = 0
    self.preprocess_time = self._preprocess.time
    self.process_time = 0
    self.total_time = 0
    self._print_summary(block_size, block_size_units)

  def task(self, index):
    ''' Get a task. '''
    return Task(
      index,
      self._experiments,
      self._profile_model,
      self._manager.split(index),
      self._manager.job(index),
      self._flatten,
      self._save_shoeboxes)

  def tasks(self):
    ''' Iterate through the tasks. '''
    for i in range(len(self)):
      yield self.task(i)

  def accumulate(self, result):
    ''' Accumulate the results. '''
    self._manager.accumulate(result.index, result.reflections)
    self.read_time += result.read_time
    self.extract_time += result.extract_time
    self.process_time += result.process_time
    self.total_time += result.total_time
    if self.finished():
      self._postprocess(self._manager.data())
      self._print_final_summary(self._manager.data())
      self.postprocess_time = self._postprocess.time

  def result(self):
    ''' Return the result. '''
    return self._manager.data()

  def finished(self):
    ''' Return if all tasks have finished. '''
    return self._manager.finished()

  def __len__(self):
    ''' Return the number of tasks. '''
    return len(self._manager)

  def _print_summary(self, block_size, block_size_units):
    ''' Print a summary of the integration stuff. '''
    from libtbx.table_utils import format as table
    from dials.util.command_line import heading

    # Create a table of integration tasks
    scans = self._experiments.scans()
    assert(len(scans) == 1)
    if scans[0] is None:
      rows = [["#", "Frame From", "Frame To"]]
      for i in range(len(self)):
        f0, f1 = self._manager.job(i)
        rows.append([str(i), str(f0), str(f1)])
    else:
      rows = [["#", "Frame From", "Frame To", "Angle From", "Angle To"]]
      for i in range(len(self)):
        f0, f1 = self._manager.job(i)
        p0 = scans[0].get_angle_from_array_index(f0)
        p1 = scans[0].get_angle_from_array_index(f1)
        rows.append([str(i), str(f0), str(f1), str(p0), str(p1)])
    task_table = table(rows, has_header=True, justify="right", prefix=" ")

    # Create summary format
    summary_format_str = (
      '%s\n'
      '\n'
      + heading('Beginning integration of the following experiments:\n') +
      '\n'
      ' Experiments: %d\n'
      ' Beams:       %d\n'
      ' Detectors:   %d\n'
      ' Goniometers: %d\n'
      ' Scans:       %d\n'
      ' Crystals:    %d\n'
      '\n'
      'Integrating reflections in the following blocks of images:\n'
      '\n'
      ' block_size: %d %s\n'
      '\n'
      '%s\n'
    )

    # Print the summary
    print summary_format_str % (
      '=' * 80,
      len(self._experiments),
      len(self._experiments.beams()),
      len(self._experiments.detectors()),
      len(self._experiments.goniometers()),
      len(self._experiments.scans()),
      len(self._experiments.crystals()),
      block_size,
      block_size_units,
      task_table)

  def _print_final_summary(self, data):
    ''' Print a final summary. '''
    from dials.algorithms.integration import Summary
    print Summary(
      data,
      self._experiments.imagesets()[0].get_array_range(),
      5).as_str()


class PreProcessorOsc(object):
  ''' A pre-processing class for oscillation data. '''

  def __init__(self, experiments, profile_model, min_zeta, partials, **kwargs):
    ''' Initialise the pre-processor. '''
    self.experiments = experiments
    self.profile_model = profile_model
    self.min_zeta = min_zeta
    self.partials = partials
    self.time = 0

  def __call__(self, data, jobs):
    ''' Do some pre-processing. '''
    from dials.array_family import flex
    from dials.util.command_line import heading
    from time import time
    print '=' * 80
    print ''
    print heading('Pre-processing reflections')
    print ''
    st = time()

    # Compute some reflection properties
    data.compute_zeta_multi(self.experiments)
    data.compute_d(self.experiments)
    data.compute_bbox(self.experiments, self.profile_model)

    # Optionally split the reflection table into partials
    if self.partials:
      num_full = len(data)
      data.split_partials()
      num_partial = len(data)
      assert(num_partial >= num_full)
      if (num_partial > num_full):
        print ' Split %d reflections into %d partial reflections\n' % (
          num_full,
          num_partial)
    else:
      num_full = len(data)
      data.split_blocks(jobs)
      num_partial = len(data)
      assert(num_partial >= num_full)
      if (num_partial > num_full):
        num_split = num_partial - num_full
        print ' Split %d reflections overlapping job boundaries\n' % num_split

    # Compute the partiality
    data.compute_partiality(self.experiments, self.profile_model)

    # Filter the reflections by zeta
    mask = flex.abs(data['zeta']) < self.min_zeta
    num_ignore = mask.count(True)
    data.set_flags(mask, data.flags.dont_integrate)
    EPS = 1e-7
    fully_recorded = data['partiality'] > (1.0 - EPS)
    num_partial = fully_recorded.count(False)
    num_full = fully_recorded.count(True)
    num_integrate = data.get_flags(data.flags.dont_integrate).count(False)
    num_reference = data.get_flags(data.flags.reference_spot).count(True)
    num_ice_ring = data.get_flags(data.flags.in_powder_ring).count(True)
    self.time = time() - st
    print ' Number of reflections'
    print '  Partial:     %d' % num_partial
    print '  Full:        %d' % num_full
    print '  In ice ring: %d' % num_ice_ring
    print '  Integrate:   %d' % num_integrate
    print '  Reference:   %d' % num_reference
    print '  Ignore:      %d' % num_ignore
    print '  Total:       %d' % len(data)
    print ''
    print ' Filtered %d reflections by zeta = %0.3f' % (num_ignore, self.min_zeta)
    print ''
    print ' Time taken: %.2f seconds' % self.time
    print ''


class PreProcessorStills(object):
  ''' A pre-processing class for stills data. '''

  def __init__(self, experiments, profile_model, **kwargs):
    ''' Initialise the pre-processor. '''
    self.experiments = experiments
    self.profile_model = profile_model
    self.time = 0

  def __call__(self, data, jobs):
    ''' Do some pre-processing. '''
    from dials.array_family import flex
    from dials.util.command_line import heading
    from time import time
    print '=' * 80
    print ''
    print heading('Pre-processing reflections')
    print ''
    st = time()

    # Compute some reflection properties
    data.compute_d(self.experiments)
    data.compute_bbox(self.experiments, self.profile_model)

    # Check the bounding boxes are all 1 frame in width
    z0, z1 = data['bbox'].parts()[4:6]
    assert((z1 - z0).all_eq(1))

    # Compute the partiality
    data.compute_partiality(self.experiments, self.profile_model)

    # Print out the pre-processing summary
    num_ignore = 0
    EPS = 1e-7
    fully_recorded = data['partiality'] > (1.0 - EPS)
    num_partial = fully_recorded.count(False)
    num_full = fully_recorded.count(True)
    num_integrate = data.get_flags(data.flags.dont_integrate).count(False)
    num_reference = data.get_flags(data.flags.reference_spot).count(True)
    num_ice_ring = data.get_flags(data.flags.in_powder_ring).count(True)
    self.time = time() - st
    print ' Number of reflections'
    print '  Partial:     %d' % num_partial
    print '  Full:        %d' % num_full
    print '  In ice ring: %d' % num_ice_ring
    print '  Integrate:   %d' % num_integrate
    print '  Reference:   %d' % num_reference
    print '  Ignore:      %d' % num_ignore
    print '  Total:       %d' % len(data)
    print ''
    print ' Time taken: %.2f seconds' % self.time
    print ''


class PostProcessorOsc(object):
  ''' A post-processing class for oscillation data. '''

  def __init__(self, experiments):
    ''' Initialise the post processor. '''
    self.experiments = experiments
    self.time = 0

  def __call__(self, data):
    ''' Do some post processing. '''
    from dials.util.command_line import heading
    from time import time
    st = time()
    print '=' * 80
    print ''
    print heading('Post-processing reflections')
    print ''

    # Compute the corrections
    data.compute_corrections(self.experiments)

    print ''
    num = data.get_flags(data.flags.integrated, all=False).count(True)
    self.time = time() - st
    print ' Integrated %d reflections' % num
    print ''
    print ' Time taken: %.2f seconds' % self.time
    print ''


class PostProcessorStills(object):
  ''' A post-processing class for stills data. '''

  def __init__(self, experiments):
    ''' Initialise the post processor. '''
    self.experiments = experiments
    self.time = 0

  def __call__(self, data):
    ''' Do some post processing. '''
    pass


class ManagerOsc(Manager):
  ''' Specialize the manager for oscillation data using the oscillation pre and
  post processors. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=10,
               block_size_units='degrees',
               min_zeta=0.05,
               flatten=False,
               partials=False,
               save_shoeboxes=False,
               **kwargs):
    ''' Initialise the pre-processor, post-processor and manager. '''

    # Ensure we have the correct type of data
    if not experiments.all_sweeps():
      raise RuntimeError('''
        An inappropriate integration algorithm may have been selected!
         Trying to perform rotation integration when not all experiments
         are indicated as rotation experiments.
      ''')

    # Create the pre-processor
    preprocess = PreProcessorOsc(
      experiments,
      profile_model,
      min_zeta=min_zeta,
      partials=partials)

    # Create the post-processor
    postprocess = PostProcessorOsc(experiments)

    # Initialise the manager
    super(ManagerOsc, self).__init__(
      preprocess,
      postprocess,
      experiments,
      profile_model,
      reflections,
      block_size=block_size,
      block_size_units=block_size_units,
      flatten=flatten,
      save_shoeboxes=save_shoeboxes)


class ManagerStills(Manager):
  ''' Specialize the manager for stills data using the stills pre and post
  processors. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               save_shoeboxes=False,
               **kwargs):
    ''' Initialise the pre-processor, post-processor and manager. '''

    # Ensure we have the correct type of data
    if not experiments.all_stills():
      raise RuntimeError('''
        An inappropriate integration algorithm may have been selected!
         Trying to perform stills integration when not all experiments
         are indicated as stills experiments.
      ''')

    # Create the pre-processor
    preprocess = PreProcessorStills(
      experiments,
      profile_model)

    # Create the post-processor
    postprocess = PostProcessorStills(experiments)

    # Initialise the manager
    super(ManagerStills, self).__init__(
      preprocess,
      postprocess,
      experiments,
      profile_model,
      reflections,
      block_size=1,
      block_size_units='frames',
      flatten=False,
      save_shoeboxes=save_shoeboxes)


class Integrator3D(Integrator):
  ''' Top level integrator for 3D integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=1,
               min_zeta=0.05,
               nproc=1,
               mp_method='multiprocessing',
               save_shoeboxes=False):
    ''' Initialise the manager and the integrator. '''

    # Create the integration manager
    manager = ManagerOsc(
      experiments,
      profile_model,
      reflections,
      block_size=block_size,
      min_zeta=min_zeta,
      save_shoeboxes=save_shoeboxes)

    # Initialise the integrator
    super(Integrator3D, self).__init__(manager, nproc, mp_method)


class IntegratorFlat3D(Integrator):
  ''' Top level integrator for flat 2D integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=1,
               min_zeta=0.05,
               nproc=1,
               mp_method='multiprocessing',
               save_shoeboxes=False):
    ''' Initialise the manager and the integrator. '''

    # Create the integration manager
    manager = ManagerOsc(
      experiments,
      profile_model,
      reflections,
      block_size=block_size,
      min_zeta=min_zeta,
      flatten=True,
      save_shoeboxes=save_shoeboxes)

    # Initialise the integrator
    super(IntegratorFlat3D, self).__init__(manager, nproc, mp_method)


class Integrator2D(Integrator):
  ''' Top level integrator for 2D integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=1,
               min_zeta=0.05,
               nproc=1,
               mp_method='multiprocessing',
               save_shoeboxes=False):
    ''' Initialise the manager and the integrator. '''

    # Create the integration manager
    manager = ManagerOsc(
      experiments,
      profile_model,
      reflections,
      block_size=block_size,
      min_zeta=min_zeta,
      partials=True,
      save_shoeboxes=save_shoeboxes)

    # Initialise the integrator
    super(Integrator2D, self).__init__(manager, nproc, mp_method)


class IntegratorSingle2D(Integrator):
  ''' Top level integrator for still image integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=1,
               min_zeta=0.05,
               nproc=1,
               mp_method='multiprocessing',
               save_shoeboxes=False):
    ''' Initialise the manager and the integrator. '''

    # Override block size
    block_size = 1

    # Create the integration manager
    manager = ManagerOsc(
      experiments,
      profile_model,
      reflections,
      block_size=block_size,
      block_size_units='frames',
      min_zeta=min_zeta,
      partials=True,
      save_shoeboxes=save_shoeboxes)

    # Initialise the integrator
    super(IntegratorSingle2D, self).__init__(manager, nproc, mp_method)


class IntegratorStills(Integrator):
  ''' Top level integrator for still image integration. '''

  def __init__(self,
               experiments,
               profile_model,
               reflections,
               block_size=1,
               min_zeta=0.05,
               nproc=1,
               mp_method='multiprocessing',
               save_shoeboxes=False):
    ''' Initialise the manager and the integrator. '''

    # Create the integration manager
    manager = ManagerStills(
      experiments,
      profile_model,
      reflections,
      save_shoeboxes=save_shoeboxes)

    # Initialise the integrator
    super(IntegratorStills, self).__init__(manager, nproc, mp_method)


class IntegratorFactory(object):
  ''' A factory for creating integrators. '''

  @staticmethod
  def create(params, experiments, profile_model, reflections):
    ''' Create the integrator from the input configuration. '''
    from dials.interfaces import IntensityIface
    from dials.interfaces import BackgroundIface
    from dials.interfaces import CentroidIface
    from dials.array_family import flex

    # Check the input
    assert(len(experiments) == len(profile_model))

    # Initialise the strategy classes
    BackgroundAlgorithm = BackgroundIface.extension(
      params.integration.background.algorithm)
    IntensityAlgorithm = IntensityIface.extension(
      params.integration.intensity.algorithm)
    CentroidAlgorithm = CentroidIface.extension(
      params.integration.centroid.algorithm)

    # Set the algorithms in the reflection table
    flex.reflection_table._background_algorithm = \
      flex.strategy(BackgroundAlgorithm, params)
    flex.reflection_table._intensity_algorithm = \
      flex.strategy(IntensityAlgorithm, params)
    flex.reflection_table._centroid_algorithm = \
      flex.strategy(CentroidAlgorithm, params)

    # Get the integrator class
    IntegratorClass = IntegratorFactory.select_integrator(
      IntensityAlgorithm.type(params, experiments))

    # Return an instantiation of the class
    return IntegratorClass(
      experiments=experiments,
      profile_model=profile_model,
      reflections=reflections,
      block_size=params.integration.block.size,
      min_zeta=params.integration.filter.min_zeta,
      nproc=params.integration.mp.nproc,
      mp_method=params.integration.mp.method,
      save_shoeboxes=params.integration.debug.save_shoeboxes)

  @staticmethod
  def select_integrator(integrator_type):
    ''' Select the integrator. '''
    if integrator_type == '3d':
      IntegratorClass = Integrator3D
    elif integrator_type == 'flat3d':
      IntegratorClass = IntegratorFlat3D
    elif integrator_type == '2d':
      IntegratorClass = Integrator2D
    elif integrator_type == 'single2d':
      IntegratorClass = IntegratorSingle2D
    elif integrator_type == 'stills':
      IntegratorClass = IntegratorStills
    else:
      raise RuntimeError("Unknown integration type")
    return IntegratorClass
