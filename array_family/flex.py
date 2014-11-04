#
# flex.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
import boost.python
from dials.model import data
from dials_array_family_flex_ext import *
from cctbx.array_family.flex import *
from cctbx.array_family import flex

# Set the 'real' type to either float or double
if get_real_type() == "float":
  real = flex.float
elif get_real_type() == "double":
  real = flex.double
else:
  raise TypeError('unknown "real" type')


class ImageSummary(object):
  ''' Summary per image. '''

  def __init__(self, data, experiment):
    ''' Compute stats. '''

    # Check some table columns
    assert("bbox" in data)
    assert("partiality" in data)
    assert("intensity.sum.value" in data)
    assert("intensity.sum.variance" in data)

    # Get the array range
    try:
      array_range = experiment.imageset.get_array_range()
    except:
      array_range = (0, len(experiment.imageset))

    # Get arrays for each frame
    data = data.select(data.get_flags(data.flags.integrated, all=False))
    data.split_partials()
    frames = data['bbox'].parts()[4]

    # Create the binner with the bins per image
    binner = Binner(flex.int(range(*array_range)).as_double())

    # Get the bins
    self.bins = binner.bins()

    # Get full and partial counts
    full = data['partiality'] > 0.997300203937
    bin_indexer = binner.indexer(frames.as_double())
    self.full = bin_indexer.sum(full.as_double())
    self.part = bin_indexer.sum((~full).as_double())

    # Get stuff from table for summation
    i_sum_flg = data.get_flags(data.flags.integrated_sum)
    i_sum_val = data['intensity.sum.value'].select(i_sum_flg)
    i_sum_var = data['intensity.sum.variance'].select(i_sum_flg)
    assert(i_sum_var.all_gt(0))
    ios_sum = i_sum_val / flex.sqrt(i_sum_var)
    bin_indexer = binner.indexer(frames.select(i_sum_flg).as_double())
    self.num_sum = bin_indexer.count()
    self.ios_sum = bin_indexer.mean(ios_sum)

    # Get stuff from table for profile fitting
    try:
      i_prf_flg = data.get_flags(data.flags.integrated_prf)
      i_prf_val = data['intensity.prf.value'].select(i_prf_flg)
      i_prf_var = data['intensity.prf.variance'].select(i_prf_flg)
      assert(i_prf_var.all_gt(0))
      ios_prf = i_prf_val / flex.sqrt(i_prf_var)
      bin_indexer = binner.indexer(frames.select(i_prf_flg).as_double())
      self.num_prf = bin_indexer.count()
      self.ios_prf = bin_indexer.mean(ios_prf)
    except RuntimeError:
      self.num_prf = flex.size_t(len(self.bins), 0)
      self.ios_prf = flex.size_t(len(self.bins), 0)

  def __len__(self):
    return len(self.bins)

  def table(self):
    ''' Produce a table of results. '''
    from libtbx.table_utils import format as table
    rows = [["Image",
             "# full",
             "# part",
             "# sum",
             "# prf",
             "<I/sigI>\n (sum)",
             "<I/sigI>\n (prf)"]]
    for i in range(len(self)):
      rows.append([
        '%d' % self.bins[i],
        '%d' % self.full[i],
        '%d' % self.part[i],
        '%d' % self.num_sum[i],
        '%d' % self.num_prf[i],
        '%.1f' % self.ios_sum[i],
        '%.1f' % self.ios_prf[i]])
    return table(rows, has_header=True, justify='right', prefix=' ')


class ResolutionSummary(object):

  def __init__(self, data, experiment, nbins=10):
    from cctbx import miller
    from cctbx import crystal

    # Check some table columns
    assert("d" in data)
    assert("intensity.sum.value" in data)
    assert("intensity.sum.variance" in data)

    # Select integrated reflections
    data = data.select(data.get_flags(data.flags.integrated, all=False))

    # Create the crystal symmetry object
    cs = crystal.symmetry(
      space_group=experiment.crystal.get_space_group(),
      unit_cell=experiment.crystal.get_unit_cell())

    # Create the array of bins
    ms = miller.set(cs, data['miller_index'])
    ms.setup_binner(n_bins=nbins)
    binner = ms.binner()
    brange = list(binner.range_used())
    bins = [binner.bin_d_range(brange[0])[0]]
    for i in brange:
      bins.append(binner.bin_d_range(i)[1])
    bins = flex.double(reversed(bins))

    # Create the binner
    binner = Binner(bins)

    # Get the bins
    self.bins = binner.bins()

    # Get full and partial counts
    full = data['partiality'] > 0.997300203937
    bin_indexer = binner.indexer(data['d'])
    self.num_full = bin_indexer.sum(full.as_double())
    self.num_part = bin_indexer.sum((~full).as_double())

    # Get stuff from table for summation
    i_sum_flg = data.get_flags(data.flags.integrated_sum)
    i_sum_val = data['intensity.sum.value'].select(i_sum_flg)
    i_sum_var = data['intensity.sum.variance'].select(i_sum_flg)
    assert(i_sum_var.all_gt(0))
    ios_sum = i_sum_val / flex.sqrt(i_sum_var)
    bin_indexer = binner.indexer(data['d'].select(i_sum_flg))
    self.num_sum = bin_indexer.count()
    self.ios_sum = bin_indexer.mean(ios_sum)

    # Get stuff from table for profile fitting
    try:
      i_prf_flg = data.get_flags(data.flags.integrated_prf)
      i_prf_val = data['intensity.prf.value'].select(i_prf_flg)
      i_prf_var = data['intensity.prf.variance'].select(i_prf_flg)
      assert(i_prf_var.all_gt(0))
      ios_prf = i_prf_val / flex.sqrt(i_prf_var)
      bin_indexer = binner.indexer(data['d'].select(i_prf_flg))
      self.num_prf = bin_indexer.count()
      self.ios_prf = bin_indexer.mean(ios_prf)
    except Exception:
      self.num_prf = flex.size_t(len(bins)-1, 0)
      self.ios_prf = flex.double(len(bins)-1, 0)

  def __len__(self):
    return len(self.bins)-1

  def table(self):
    from libtbx.table_utils import format as table
    rows = [["d min",
             "d max",
             "# full",
             "# part",
             "# sum",
             "# prf",
             "<I/sigI>\n (sum)",
             "<I/sigI>\n (prf)"]]
    for i in range(len(self)):
      rows.append([
        '%.1f' % self.bins[i],
        '%.1f' % self.bins[i+1],
        '%d'   % self.num_full[i],
        '%d'   % self.num_part[i],
        '%d'   % self.num_sum[i],
        '%d'   % self.num_prf[i],
        '%.1f' % self.ios_sum[i],
        '%.1f' % self.ios_prf[i]])
    return table(rows, has_header=True, justify='right', prefix=' ')


class WholeSummary(object):
  ''' Whole dataset summary. '''

  def __init__(self, data, experiment):
    ''' Compute the results. '''

    # Compute for summation
    flags_sum = data.get_flags(data.flags.integrated_sum)
    I_sum_val = data['intensity.sum.value'].select(flags_sum)
    I_sum_var = data['intensity.sum.variance'].select(flags_sum)
    assert(I_sum_var.all_gt(0))
    self.sum_ios = flex.mean(I_sum_val / flex.sqrt(I_sum_var))

    # Compute for profile fitting
    try:
      flags_prf = data.get_flags(data.flags.integrated_prf)
      I_prf_val = data['intensity.prf.value'].select(flags_prf)
      I_prf_var = data['intensity.prf.variance'].select(flags_prf)
      assert(I_prf_var.all_gt(0))
      self.prf_ios = flex.mean(I_prf_val / flex.sqrt(I_prf_var))
    except Exception:
        self.prf_ios = 0.0

  def table(self):
    ''' Produce a table of results. '''
    from libtbx.table_utils import format as table
    rows = [["<I/sigI>\n (sum)",
             "<I/sigI>\n (prf)"]]
    rows.append([
      '%.1f' % self.sum_ios,
      '%.1f' % self.prf_ios])
    return table(rows, has_header=True, justify='right', prefix=' ')


class Summary(object):
  ''' A class to present a summary of integration results. '''

  def __init__(self, data, experiment):
    ''' Initialise. '''
    self._image_summary = ImageSummary(data, experiment)
    self._resolution_summary = ResolutionSummary(data, experiment)
    self._whole_summary = WholeSummary(data, experiment)

  def __str__(self):
    ''' Return as a string. '''
    from dials.util.command_line import heading
    img_summary = self._image_summary.table()
    res_summary = self._resolution_summary.table()
    who_summary = self._whole_summary.table()
    print '=' * 80
    print ''
    print heading('Summary of integration results')
    print ''
    return (
      '%s\n'
      '\n'
      '%s\n'
      '\n'
      ' Summary of integration results as a function of image number'
      '\n%s\n\n'
      ' Summary of integration results binned by resolution'
      '\n%s\n\n'
      ' Summary of integration results for the whole dataset'
      '\n%s\n'
    ) % ('=' * 80,
         heading('Summary of integration results'),
         img_summary,
         res_summary,
         who_summary)


# class ImageSummaryAux(boost.python.injector, ImageSummary):

#   def table(self, index):
#     from libtbx.table_utils import format as table
#     data = self.data(index)
#     rows = [["Image",
#              "# full",
#              "# part",
#              "<I/sigI>\n (sum)",
#              "<I/sigI>\n (prf)"]]
#     for i in range(len(data)):
#       rows.append([
#         '%d' % data[i].frame,
#         '%d' % data[i].full,
#         '%d' % data[i].part,
#         '%.1f' % data[i].sum_ios,
#         '%.1f' % data[i].prf_ios])
#     return table(rows, has_header=True, justify='right', prefix=' ')

#   def tables(self):
#     t = []
#     for i in range(len(self)):
#       t.append('%s\n%s' % (
#         ' Experiment: %d' % i,
#         self.table(i)))
#     return '\n'.join(t)


# class ResolutionSummaryAux(boost.python.injector, ResolutionSummary):

#   def table(self, index):
#     from libtbx.table_utils import format as table
#     data = self.data(index)
#     rows = [["d",
#              "<I/sigI>\n (sum)",
#              "<I/sigI>\n (prf)"]]
#     for i in range(len(data)):
#       rows.append([
#         '%.1f -> %.1f' % (data[i].d[0], data[i].d[1]),
#         '%.1f' % data[i].sum_ios,
#         '%.1f' % data[i].prf_ios])
#     return table(rows, has_header=True, justify='right', prefix=' ')

#   def tables(self):
#     t = []
#     for i in range(len(self)):
#       t.append('%s\n%s' % (
#         ' Experiment: %d' % i,
#         self.table(i)))
#     return '\n'.join(t)


# class WholeSummaryAux(boost.python.injector, WholeSummary):

#   def table(self, index):
#     from libtbx.table_utils import format as table
#     data = self.data(index)
#     rows = [["<I/sigI>\n (sum)",
#              "<I/sigI>\n (prf)"]]
#     rows.append([
#       '%.1f' % data.sum_ios,
#       '%.1f' % data.prf_ios])
#     return table(rows, has_header=True, justify='right', prefix=' ')

#   def tables(self):
#     t = []
#     for i in range(len(self)):
#       t.append('%s\n%s' % (
#         ' Experiment: %d' % i,
#         self.table(i)))
#     return '\n'.join(t)


# class SummaryAux(boost.python.injector, Summary):

#   def __str__(self):
#     from dials.util.command_line import heading
#     img_summary = self.image_summary().tables()
#     res_summary = self.resolution_summary().tables()
#     who_summary = self.whole_summary().tables()
#     print '=' * 80
#     print ''
#     print heading('Summary of integration results')
#     print ''
#     return (
#       '%s\n'
#       '\n'
#       '%s\n'
#       '\n'
#       ' Summary of integration results as a function of image number'
#       '\n%s\n\n'
#       ' Summary of integration results binned by resolution'
#       '\n%s\n\n'
#       ' Summary of integration results for the whole dataset'
#       '\n%s\n'
#     ) % ('=' * 80,
#          heading('Summary of integration results'),
#          img_summary,
#          res_summary,
#          who_summary)


def strategy(cls, params=None):
  ''' Wrap a class that takes params and experiments as a strategy. '''
  from functools import wraps
  @wraps(cls)
  def call(self, *args):
    return cls(params, *args)
  return call

def default_background_algorithm():
  ''' Get the default background algorithm. '''
  from dials.extensions import SimpleBackgroundExt
  return strategy(SimpleBackgroundExt)

def default_intensity_algorithm():
  ''' Get the default intensity algorithm. '''
  from dials.extensions import SummationIntegrationExt
  return strategy(SummationIntegrationExt)

def default_centroid_algorithm():
  ''' Get the default centroid algorithm. '''
  from dials.extensions import SimpleCentroidExt
  return strategy(SimpleCentroidExt)


class reflection_table_aux(boost.python.injector, reflection_table):
  ''' An injector class to add additional methods to the reflection table. '''

  # Set the default algorithms. These are set as class variables so that if they
  # are changed in the class, all new instances of reflection table will have
  # the modified algorithms. If these are modified on the instance level, then
  # only the instance will have the modified algorithms and new instances will
  # have the defaults
  _background_algorithm = default_background_algorithm()
  _intensity_algorithm = default_intensity_algorithm()
  _centroid_algorithm = default_centroid_algorithm()

  @staticmethod
  def from_predictions(experiment,
                       dmin=None,
                       dmax=None,
                       margin=1,
                       force_static=False):
    ''' Construct a reflection table from predictions. '''
    from dials.algorithms.spot_prediction.reflection_predictor \
      import ReflectionPredictor
    predict = ReflectionPredictor(
      experiment,
      dmin=dmin,
      dmax=dmax,
      margin=margin,
      force_static=force_static)
    return predict()

  @staticmethod
  def from_predictions_multi(experiments,
                             dmin=None,
                             dmax=None,
                             margin=1,
                             force_static=False):
    ''' Construct a reflection table from predictions. '''
    from scitbx.array_family import flex
    result = reflection_table()
    for i, e in enumerate(experiments):
      rlist = reflection_table.from_predictions(
        e,
        dmin=dmin,
        dmax=dmax,
        margin=margin,
        force_static=force_static)
      rlist['id'] = flex.size_t(len(rlist), i)
      result.extend(rlist)
    return result

  @staticmethod
  def from_observations(datablock, params=None):
    ''' Construct a reflection table from observations. '''
    from dials.algorithms.peak_finding.spotfinder_factory \
      import SpotFinderFactory

    # Get the integrator from the input parameters
    print 'Configuring spot finder from input parameters'
    find_spots = SpotFinderFactory.from_parameters(params)

    # Find the spots
    return find_spots(datablock)

  @staticmethod
  def from_pickle(filename):
    ''' Read the reflection table from pickle file. '''
    import cPickle as pickle
    with open(filename, 'rb') as infile:
      result = pickle.load(infile)
      assert(isinstance(result, reflection_table))
      return result

  @staticmethod
  def from_h5(filename):
    ''' Read the reflections table from a HDF5 file. '''
    from dials.util.nexus import NexusFile
    handle = NexusFile(filename, 'r')
    self = handle.get_reflections()
    handle.close()
    return self

  @staticmethod
  def empty_standard(nrows):
    ''' Create an empty table of specified number of rows with most of the
    standard keys'''

    assert nrows > 0
    table = reflection_table(nrows)

    # General properties
    table['flags'] = flex.size_t(nrows, 0)
    table['id'] = flex.double(nrows, 0)
    table['panel'] = flex.size_t(nrows, 0)

    # Predicted properties
    table['miller_index'] = flex.miller_index(nrows)
    table['entering'] = flex.bool(nrows)
    table['s1'] = flex.vec3_double(nrows, (0, 0, 0))
    table['xyzcal.mm'] = flex.vec3_double(nrows, (0, 0, 0))
    table['xyzcal.px'] = flex.vec3_double(nrows, (0, 0, 0))
    #table['ub_matrix'] = flex.mat3_double(nrows, (0, 0, 0, 0, 0, 0, 0, 0, 0))

    # Observed properties
    table['xyzobs.px.value'] = flex.vec3_double(nrows, (0, 0, 0))
    table['xyzobs.px.variance'] = flex.vec3_double(nrows, (0, 0, 0))
    table['xyzobs.mm.value'] = flex.vec3_double(nrows, (0, 0, 0))
    table['xyzobs.mm.variance'] = flex.vec3_double(nrows, (0, 0, 0))
    table['rlp'] = flex.vec3_double(nrows, (0, 0, 0))
    table['intensity.sum.value'] = flex.double(nrows, 0)
    table['intensity.sum.variance'] = flex.double(nrows, 0)
    table['intensity.prf.value'] = flex.double(nrows, 0)
    table['intensity.prf.variance'] = flex.double(nrows, 0)
    table['lp'] = flex.double(nrows, 0)
    table['profile.correlation'] = flex.double(nrows, 0)

    return table

  def as_pickle(self, filename):
    ''' Write the reflection table as a pickle file. '''
    import cPickle as pickle
    with open(filename, 'wb') as outfile:
      pickle.dump(self, outfile, protocol=pickle.HIGHEST_PROTOCOL)

  def as_h5(self, filename):
    ''' Write the reflection table as a HDF5 file. '''
    from dials.util.nexus import NexusFile
    handle = NexusFile(filename, 'w')
    handle.set_reflections(self)
    handle.close()

  def copy(self):
    ''' Copy everything. '''
    from scitbx.array_family import flex
    return self.select(flex.bool(len(self), True))

  def sort(self, name, reverse=False):
    ''' Sort the reflection table by a key. '''
    import __builtin__
    column = self[name]
    indices = __builtin__.sorted(
      range(len(self)),
      key=lambda x: column[x],
      reverse=reverse)
    self.reorder(flex.size_t(indices))

  def match(self, other):
    ''' Match reflections with another set of reflections. '''
    from dials.algorithms.peak_finding.spot_matcher import SpotMatcher
    match = SpotMatcher(max_separation=1)
    oind, sind = match(other, self)
    return sind, oind

  def match_with_reference(self, other):
    ''' Match reflections with another set of reflections. '''
    from dials.util.command_line import Command
    Command.start("Matching reference spots with predicted reflections")
    sind, oind = self.match(other)
    h1 = self.select(sind)['miller_index']
    h2 = other.select(oind)['miller_index']
    mask = (h1 == h2)
    self.set_flags(sind.select(mask), self.flags.reference_spot)
    Command.end("Matched %d reference spots with predicted reflections" %
                mask.count(True))

  #def is_bbox_inside_image_range(self, experiment):
    #''' Check if bbox is within image range. '''
    #from dials.algorithms import filtering
    #assert(len(experiment.detector) == 1)
    #return filtering.is_bbox_outside_image_range(
      #self['bbox'],
      #experiment.detector[0].get_image_size()[::-1],
      #experiment.scan.get_array_range()) != True

  def compute_zeta(self, experiment):
    ''' Compute zeta for each reflection. '''
    from dials.algorithms.profile_model.gaussian_rs import zeta_factor
    m2 = experiment.goniometer.get_rotation_axis()
    s0 = experiment.beam.get_s0()
    self['zeta'] = zeta_factor(m2, s0, self['s1'])
    return self['zeta']

  def compute_zeta_multi(self, experiments):
    ''' Compute zeta for each reflection. '''
    from dials.algorithms.profile_model.gaussian_rs import zeta_factor
    m2 = flex.vec3_double(len(experiments))
    s0 = flex.vec3_double(len(experiments))
    for i, e in enumerate(experiments):
      m2[i] = e.goniometer.get_rotation_axis()
      s0[i] = e.beam.get_s0()
    self['zeta'] = zeta_factor(m2, s0, self['s1'], self['id'])
    return self['zeta']

  def compute_d_single(self, experiment):
    ''' Compute the resolution for each reflection. '''
    from dials.array_family import flex
    uc = flex.unit_cell(1)
    uc[0] = experiment.crystal.get_unit_cell()
    self['d'] = uc.d(self['miller_index'], flex.size_t(len(self), 0))
    return self['d']

  def compute_d(self, experiments):
    ''' Compute the resolution for each reflection. '''
    from dials.array_family import flex
    uc = flex.unit_cell(len(experiments))
    for i, e in enumerate(experiments):
      uc[i] = e.crystal.get_unit_cell()
    self['d'] = uc.d(self['miller_index'], self['id'])
    return self['d']

  def compute_bbox(self, experiments, profile_model, sigma_b_multiplier=2.0):
    ''' Compute the bounding boxes. '''
    from dials.util.command_line import Command
    profile_model.compute_bbox(experiments, self, sigma_b_multiplier)

  def compute_partiality(self, experiments, profile_model):
    ''' Compute the reflection partiality. '''
    profile_model.compute_partiality(experiments, self)

  def compute_background(self, experiments):
    ''' Helper function to compute the background. '''
    self._background_algorithm(experiments).compute_background(self)

  def compute_centroid(self, experiments):
    ''' Helper function to compute the centroid. '''
    self._centroid_algorithm(experiments).compute_centroid(self)

  def compute_intensity(self, experiments, profile_model):
    ''' Helper function to compute the intensity. '''
    self._intensity_algorithm(experiments, profile_model).compute_intensity(self)

  def compute_summed_intensity(self):
    ''' Compute intensity via summation integration. '''
    from dials.algorithms.integration.sum import IntegrationAlgorithm
    algorithm = IntegrationAlgorithm()
    algorithm(self)

  def compute_corrections(self, experiments):
    ''' Helper function to correct the intensity. '''
    from dials.algorithms.integration import Corrections, CorrectionsMulti
    from dials.util.command_line import Command
    Command.start("Calculating lp correction")
    compute = CorrectionsMulti()
    for experiment in experiments:
      compute.append(Corrections(
        experiment.beam,
        experiment.goniometer))
    lp = compute.lp(self['id'], self['s1'])
    self['lp'] = lp
    Command.end("Calculated lp correction")
    return lp

  def integrate(self, experiments, profile_model):
    ''' Helper function to integrate reflections. '''
    self.compute_background(experiments)
    self.compute_centroid(experiments)
    self.compute_summed_intensity()
    self.compute_intensity(experiments, profile_model)

  def compute_mask(self, experiments, profile_model):
    ''' Apply a mask to the shoeboxes. '''
    from dials.algorithms.shoebox import Masker3DProfile
    assert(len(experiments) == 1)
    mask_profiles = Masker3DProfile(experiments, profile_model)
    mask_profiles(self, None)

  def extract_shoeboxes(self, imageset, mask=None):
    ''' Helper function to read a load of shoebox data. '''
    from dials.model.data import make_image
    import sys
    from time import time
    assert("shoebox" in self)
    detector = imageset.get_detector()
    try:
      frame0, frame1 = imageset.get_array_range()
    except Exception:
      frame0, frame1 = (0, len(imageset))
    extractor = ShoeboxExtractor(self, len(detector), frame0, frame1)
    if mask is None:
      image = imageset[0]
      if not isinstance(image, tuple):
        image = (image,)
      mask = []
      for i in range(len(image)):
        tr = detector[i].get_trusted_range()
        mask.append(image[i].as_double() > tr[0])
      mask = tuple(mask)
    sys.stdout.write("Reading images: ")
    read_time = 0
    extract_time = 0
    for i in range(len(imageset)):
      st = time()
      image = imageset[i]
      read_time += time() - st
      if not isinstance(image, tuple):
        image = (image,)
      st = time()
      extractor.next(make_image(image, mask))
      extract_time += time() - st
      sys.stdout.write(".")
      sys.stdout.flush()
      del image
    sys.stdout.write("\n")
    sys.stdout.flush()
    assert(extractor.finished())
    return read_time, extract_time

  def statistics(self, experiments):
    ''' Return some simple statistics. '''
    tables = self.split_by_experiment_id()
    assert(len(tables) == len(experiments))
    summaries = []
    for table, experiment in zip(tables, experiments):
      summaries.append(Summary(table, experiment))
    return summaries
