"""
Module defining methods for filtering reflection tables, which are combined
into functions to perform the relevant filtering on a reflection table, to
produce a filtered reflection table ready for export or further processing.

The set of classes defined in this module have filtering methods implemented as
classmethods/staticmethods, to allow easy use of individual methods. The
different classes are to handle filtering of different intensity types - profile,
scale, sum, profile + sum, etc. All functions and classmethods/staticmethods act on
a reflection table, returning a reflection table that is typically a new object,
due to the use of flex selections.

Functions:
  - filter_reflection_table:
      performs a full filtering algorithm for a given intensity choice
  - sum_partial_reflections:
      combines matching partials, replacing them with a single combined value

  filter_reflection_table takes in the following parameters: min_isigi=float,
  filter_ice_rings=bool, combine_partials=bool, partiality_threshold=float,
  intensity_choice=strings (passed in as a list e.g. ['sum', 'profile'])

Classes:
  - FilteringReductionMethods:
      a collection of staticmethods applicable to any kind of intensity values
  - FilterForExportAlgorithm:
      defines a full, general filtering algorithm for a reflection table
  - PrfIntensityReducer:
      implements methods specific to filtering of profile fitted (prf) intensities
  - SumIntensityReducer
      implements methods specific to filtering of summation (sum) intensities
  - SumAndPrfIntensityReducer
      implements filtering methods when using prf intensities if present, else
      sum (per reflection)
  - ScaleIntensityReducer
      implements filtering methods for intensities output from scaling
  - AllSumPrfScaleIntensityReducer
      implements filtering methods for using all of prf, sum and scale intensities
"""
import logging
import abc
from collections import defaultdict
from cctbx import crystal, miller
from libtbx.utils import Sorry
from libtbx.table_utils import simple_table
from dials.array_family import flex

logger = logging.getLogger('dials')

def filter_reflection_table(reflection_table, intensity_choice, *args, **kwargs):
  """Filter the data and delete unneeded intensity columns."""
  if intensity_choice == ['scale']:
    reducer = ScaleIntensityReducer
  elif intensity_choice == ['sum']:
    reducer = SumIntensityReducer
  elif intensity_choice == ['profile']:
    reducer = PrfIntensityReducer
  elif all([i in intensity_choice for i in ['sum', 'scale', 'profile']]):
    reducer = AllSumPrfScaleIntensityReducer
  elif all([i in intensity_choice for i in ['sum', 'profile']]):
    reducer = SumAndPrfIntensityReducer
  else:
    raise Sorry(("Unrecognised intensity choice for filter_reflection_table,\n"
      "value read: {0}\n"
      "must be either: 'scale', 'profile', 'sum', 'profile sum' or 'profile sum scale'\n"
      "(if parsing from command line, multiple choices passed as e.g. profile+sum"
      ).format(intensity_choice))
  reflection_table = reducer.filter_for_export(reflection_table, *args, **kwargs)
  return reflection_table

def integrated_data_to_filtered_miller_array(reflections, exp_crystal):
  crystal_symmetry = crystal.symmetry(
    unit_cell=exp_crystal.get_unit_cell(),
    space_group=exp_crystal.get_space_group())

  if 'intensity.scale.value' in reflections:
    intensity_choice = ['scale']
    intensity_to_use = 'scale'
  else:
    assert 'intensity.sum.value' in reflections
    intensity_choice = ['sum']
    if 'intensity.prf.value' in reflections:
      intensity_choice.append('profile')
      intensity_to_use = 'prf'
    else:
      intensity_to_use = 'sum'

  reflections = filter_reflection_table(reflections, intensity_choice, min_isigi=-5,
    filter_ice_rings=False, combine_partials=True,
    partiality_threshold=0.2)
  data = reflections['intensity.'+intensity_to_use+'.value']
  variances = reflections['intensity.'+intensity_to_use+'.variance']

  miller_indices = reflections['miller_index']
  assert variances.all_gt(0)
  sigmas = flex.sqrt(variances)

  miller_set = miller.set(
    crystal_symmetry, miller_indices, anomalous_flag=True)
  intensities = miller.array(miller_set, data=data, sigmas=sigmas)
  intensities.set_observation_type_xray_intensity()
  return intensities

class FilteringReductionMethods(object):

  """A collection of methods for filtering. Some internal methods require an
  'intensity' string, which indicates which column to filter on. These
  methods can be called multiple times to filter on multiple intensity
  choices. All methods may reduce the size of the reflection table by deleting
  data."""

  @staticmethod
  def _filter_on_min_isigi(reflection_table, intensity, min_isigi=None):
    if min_isigi:
      selection = (reflection_table['intensity.'+intensity+'.value'] /
        flex.sqrt(reflection_table['intensity.'+intensity+'.variance'])) < min_isigi
      reflection_table.del_selected(selection)
      logger.info('Removing %d reflections with I/Sig(I) < %s' %(
        selection.count(True), min_isigi))
    return reflection_table

  @staticmethod
  def _filter_bad_variances(reflection_table, intensity):
    selection = reflection_table['intensity.'+intensity+'.variance'] <= 0
    if selection.count(True) > 0:
      reflection_table.del_selected(selection)
      logger.info('Removing %d %s reflections with negative variance' % \
            (selection.count(True),  'intensity.'+intensity+'.value'))
    return reflection_table

  @staticmethod
  def calculate_lp_qe_correction_and_filter(reflection_table):
    # FIXME errors in e.g. LP correction need to be propagated here?
    qe = None
    if 'qe' in reflection_table:
      reflection_table = reflection_table.select(reflection_table['qe'] > 0.0)
      qe = reflection_table['qe']
    elif 'dqe' in reflection_table:
      reflection_table = reflection_table.select(reflection_table['dqe'] > 0.0)
      qe = reflection_table['dqe']
    # Now calculate conversion factor
    conversion = flex.double(reflection_table.size(), 1.0)
    if qe:
      conversion /= qe
    if 'lp' in reflection_table:
      conversion *= reflection_table['lp']
    return reflection_table, conversion

  @staticmethod
  def filter_ice_rings(reflection_table):
    selection = reflection_table.get_flags(reflection_table.flags.in_powder_ring)
    reflection_table.del_selected(selection)
    logger.info("Removing %d reflections in ice ring resolutions" %
                selection.count(True))
    return reflection_table

  @staticmethod
  def filter_on_d_min(reflection_table, d_min):
    selection = reflection_table['d'] < d_min
    logger.info("Removing %d reflections below a resolution of %s" %
                (selection.count(True), d_min))
    reflection_table.del_selected(selection)
    return reflection_table

  @staticmethod
  def filter_unassigned_reflections(reflection_table):
    """"Select reflections that are assigned to an experiment (i.e.
    non-negative id). This step will need to be looked at again once UIDS
    are used. This should currently only affect output before the scaling step,
    as scaling assigns an id."""
    reflection_table = reflection_table.select(reflection_table['id'] >= 0)
    if reflection_table.size() == 0:
      raise Sorry("No experiment-assigned reflections were found")
    logger.info('Read %s predicted reflections' % reflection_table.size())
    return reflection_table

  @staticmethod
  def combine_and_filter_partials(reflection_table, partiality_threshold,
    combine_partials=True):
    if 'partiality' in reflection_table:
      reflection_table['fractioncalc'] = reflection_table['partiality']
      if combine_partials and 'partial_id' in reflection_table:
        dataset_ids = set(reflection_table['id'])
        n_datasets = len(dataset_ids)
        if n_datasets > 1:
          total_reflection_table = flex.reflection_table()
          for id_ in dataset_ids:
            single_table = reflection_table.select(reflection_table['id'] == id_)
            total_reflection_table.extend(sum_partial_reflections(single_table))
          reflection_table = total_reflection_table
        else:
          reflection_table = sum_partial_reflections(reflection_table)
        reflection_table['fractioncalc'] = reflection_table['partiality']
      selection = reflection_table['partiality'] < partiality_threshold
      if selection.count(True) > 0:
        reflection_table.del_selected(selection)
        logger.info('Removing %d incomplete reflections' % selection.count(True))
    else:
      reflection_table['fractioncalc'] = flex.double(reflection_table.size(), 1.0)
    return reflection_table


class FilterForExportAlgorithm(FilteringReductionMethods):

  """An abstract class that defines the filter_for_export algorithm and
  abstract methods which must be implemented in a subclass."""

  __metaclass__ = abc.ABCMeta

  allowed_intensities = ['prf', 'scale', 'sum'] #Supported intensities
  # subclasses must define a class attribute intensities, which is a list
  # of a subset of the allowed intensities.

  @classmethod
  def _filter_for_export(cls, reflection_table, min_isigi=None, filter_ice_rings=False,
    combine_partials=True, partiality_threshold=0.99, d_min=None):
    """Designed to be called by subclasses"""
    assert reflection_table.size() > 0, \
      """Empty reflection table given to reduce_data_for_export function"""

    reflection_table = cls.filter_unassigned_reflections(reflection_table)
    reflection_table = cls.reduce_on_intensities(reflection_table)

    if reflection_table.size() == 0:
      raise Sorry('No suitable reflections found for intensity choice: %s' % \
        ' '.join(['intensity.'+i+'.value' for i in cls.intensities]))

    for intensity in cls.allowed_intensities:
      if intensity not in cls.intensities:
        if 'intensity.'+intensity+'.value' in reflection_table:
          del reflection_table['intensity.'+intensity+'.value']
          del reflection_table['intensity.'+intensity+'.variance']
      else:
        msg = ('No intensity.'+intensity+' values found in reflection table')
        assert 'intensity.'+intensity+'.value' in reflection_table, msg
        assert 'intensity.'+intensity+'.variance' in reflection_table, msg

    reflection_table = cls.filter_bad_variances(reflection_table)
    if filter_ice_rings:
      reflection_table = cls.filter_ice_rings(reflection_table)
    if d_min:
      reflection_table = cls.filter_on_d_min(reflection_table, d_min)

    reflection_table = cls.apply_scaling_factors(reflection_table)

    reflection_table = cls.combine_and_filter_partials(
      reflection_table, partiality_threshold, combine_partials)

    if min_isigi: #Select on I/sigI after applying scale factors and combining partials
      reflection_table = cls.filter_on_min_isigi(reflection_table, min_isigi)

    return reflection_table

  @abc.abstractmethod
  def reduce_on_intensities(reflection_table):
    """Reduce the reflection table to contain only the desired reflections
    based on intensity choice."""

  @classmethod
  @abc.abstractmethod
  def filter_bad_variances(cls, reflection_table, intensity):
    """Remove negative/zero variances of the relevant intensities."""

  @classmethod
  @abc.abstractmethod
  def filter_on_min_isigi(cls, reflection_table, min_isigi=None):
    """Filter the given intensities on I/sigI"""

  @abc.abstractmethod
  def apply_scaling_factors(reflection_table):
    """Apply the relevent scaling factors including lp, qde, scale etc."""


class PrfIntensityReducer(FilterForExportAlgorithm):

  """A class to implement methods to reduce prf intensity data and to
  implement filtering for export"""

  intensities = ['prf']

  @classmethod
  def filter_for_export(cls, reflection_table, min_isigi=None,
    filter_ice_rings=False, combine_partials=True, partiality_threshold=0.99,
    d_min=None):
    return cls._filter_for_export(reflection_table, min_isigi, filter_ice_rings,
      combine_partials, partiality_threshold, d_min)

  @classmethod
  def filter_on_min_isigi(cls, reflection_table, min_isigi=None):
    reflection_table = cls._filter_on_min_isigi(reflection_table,
      cls.intensities[0], min_isigi)
    return reflection_table

  @classmethod
  def filter_bad_variances(cls, reflection_table):
    return cls._filter_bad_variances(reflection_table,
      cls.intensities[0])

  @staticmethod
  def reduce_on_intensities(reflection_table):
    """Select profile fitted reflectons and remove bad variances"""
    selection = reflection_table.get_flags(reflection_table.flags.integrated_prf)
    reflection_table = reflection_table.select(selection)
    logger.info("Selected %d profile integrated reflections" % reflection_table.size())
    return reflection_table

  @classmethod
  def apply_scaling_factors(cls, reflection_table):
    if 'partiality' in reflection_table:
      reflection_table = reflection_table.select(reflection_table['partiality'] > 0.0)

    reflection_table, conversion = cls.calculate_lp_qe_correction_and_filter(
      reflection_table)

    reflection_table['intensity.prf.value'] *= conversion
    reflection_table['intensity.prf.variance'] *= conversion * conversion
    return reflection_table


class SumIntensityReducer(FilterForExportAlgorithm):

  """A class to implement methods to reduce sum intensity data and to
  implement filtering for export"""

  intensities = ['sum']

  @classmethod
  def filter_for_export(cls, reflection_table, min_isigi=None,
    filter_ice_rings=False, combine_partials=True, partiality_threshold=0.99,
    d_min=None):
    return cls._filter_for_export(reflection_table, min_isigi, filter_ice_rings,
      combine_partials, partiality_threshold, d_min)

  @classmethod
  def filter_on_min_isigi(cls, reflection_table, min_isigi=None):
    return cls._filter_on_min_isigi(
      reflection_table, cls.intensities[0], min_isigi)

  @classmethod
  def filter_bad_variances(cls, reflection_table):
    return cls._filter_bad_variances(reflection_table, cls.intensities[0])

  @staticmethod
  def reduce_on_intensities(reflection_table):
    """Select integrated summation reflectons and remove bad variances"""
    selection = reflection_table.get_flags(reflection_table.flags.integrated_sum)
    reflection_table = reflection_table.select(selection)
    logger.info("Selected %d summation integrated reflections" %
      reflection_table.size())
    return reflection_table

  @classmethod
  def apply_scaling_factors(cls, reflection_table):

    reflection_table, conversion = cls.calculate_lp_qe_correction_and_filter(
      reflection_table)

    if 'partiality' in reflection_table:
      nonzero_sel = reflection_table['partiality'] > 0.0
      reflection_table = reflection_table.select(nonzero_sel)
      conversion = conversion.select(nonzero_sel)
      conversion /= reflection_table['partiality']

    reflection_table['intensity.sum.value'] *= conversion
    reflection_table['intensity.sum.variance'] *= conversion * conversion
    return reflection_table


class SumAndPrfIntensityReducer(FilterForExportAlgorithm):

  """A class to implement methods to reduce sum and prf intensity data and to
  implement filtering for export. Reflections are kept both a prf and
  sum intensity is defined."""

  intensities = ['sum', 'prf']

  @classmethod
  def filter_for_export(cls, reflection_table, min_isigi=None, filter_ice_rings=False,
      combine_partials=True, partiality_threshold=0.99, d_min=None):
    return cls._filter_for_export(reflection_table, min_isigi, filter_ice_rings,
      combine_partials, partiality_threshold, d_min)

  @classmethod
  def filter_on_min_isigi(cls, reflection_table, min_isigi=None):
    """Remove reflection if either has a IsigI below min_isigi."""
    for intensity in cls.intensities:
      reflection_table = cls._filter_on_min_isigi(reflection_table,
        intensity, min_isigi)
    return reflection_table

  @classmethod
  def filter_bad_variances(cls, reflection_table):
    """Remove reflection if either has a bad variance."""
    for intensity in cls.intensities:
      reflection_table = cls._filter_bad_variances(reflection_table,
        intensity)
    return reflection_table

  @staticmethod
  def reduce_on_intensities(reflection_table):
    """First select the reflections which have successfully been integrated by
    both methods"""
    selection = reflection_table.get_flags(reflection_table.flags.integrated,
      all=True)
    reflection_table = reflection_table.select(selection)
    if reflection_table.size() == 0:
      raise Sorry('No reflections found with both profile and sum intensities,'
        'try selecting a different intensity choice - sum perhaps?')
    logger.info("Selected %d integrated reflections" % reflection_table.size())
    return reflection_table

  @classmethod
  def apply_scaling_factors(cls, reflection_table):

    reflection_table, conversion = cls.calculate_lp_qe_correction_and_filter(
      reflection_table)
    sum_conversion = conversion

    if 'partiality' in reflection_table:
      nonzero_sel = reflection_table['partiality'] > 0.0
      reflection_table = reflection_table.select(nonzero_sel)
      conversion = conversion.select(nonzero_sel)
      sum_conversion = conversion / reflection_table['partiality']

    reflection_table['intensity.sum.value'] *= sum_conversion 
    reflection_table['intensity.sum.variance'] *= sum_conversion * sum_conversion
    reflection_table['intensity.prf.value'] *= conversion
    reflection_table['intensity.prf.variance'] *= conversion * conversion
    return reflection_table


class ScaleIntensityReducer(FilterForExportAlgorithm):

  """A class to implement methods to reduce scale intensity data and to
  implement filtering for export"""

  intensities = ['scale']

  @classmethod
  def filter_for_export(cls, reflection_table, min_isigi=None,
    filter_ice_rings=False, combine_partials=True, partiality_threshold=0.99,
    d_min=None):
    return cls._filter_for_export(reflection_table, min_isigi, filter_ice_rings,
      combine_partials, partiality_threshold, d_min)

  @classmethod
  def filter_on_min_isigi(cls, reflection_table, min_isigi=None):
    return cls._filter_on_min_isigi(
      reflection_table, cls.intensities[0], min_isigi)

  @classmethod
  def filter_bad_variances(cls, reflection_table):
    return cls._filter_bad_variances(reflection_table, cls.intensities[0])

  @staticmethod
  def reduce_on_intensities(reflection_table):
    """Select intensities used for scaling and remove scaling outliers"""
    selection = ~(reflection_table.get_flags(
      reflection_table.flags.bad_for_scaling, all=False))
    outliers = reflection_table.get_flags(reflection_table.flags.outlier_in_scaling)
    reflection_table = reflection_table.select(selection & ~outliers)
    logger.info("Selected %d scaled reflections" % reflection_table.size())
    if not 'inverse_scale_factor' in reflection_table:
      raise Sorry('Inverse scale factor required to output scale intensities.')
    selection = reflection_table['inverse_scale_factor'] <= 0.0
    if selection.count(True) > 0:
      reflection_table.del_selected(selection)
      logger.info("Removed %s reflections with zero or negative inverse scale factors"\
        % selection.count(True))
    return reflection_table

  @staticmethod
  def apply_scaling_factors(reflection_table):
    """Apply the inverse scale factor."""
    if 'partiality' in reflection_table:
      reflection_table = reflection_table.select(
        reflection_table['partiality'] > 0.0)

    if not 'inverse_scale_factor' in reflection_table:
      raise Sorry('Inverse scale factor required to output scale intensities.')
    reflection_table['intensity.scale.value'] /= \
      reflection_table['inverse_scale_factor']
    reflection_table['intensity.scale.variance'] /= \
      reflection_table['inverse_scale_factor']**2
    return reflection_table


class AllSumPrfScaleIntensityReducer(SumAndPrfIntensityReducer):

  """A class to implement methods to reduce data where both prf, sum and scale
  intensities are defined. Only reflections with valid values for all intensity
  types are retained."""

  intensities = ['sum', 'prf', 'scale']

  @staticmethod
  def reduce_on_intensities(reflection_table):
    """Select those with valid reflections for all values."""
    selection = reflection_table.get_flags(reflection_table.flags.integrated,
      all=True)
    reflection_table = reflection_table.select(selection)
    logger.info("Selected %d integrated reflections" % reflection_table.size())
    reflection_table = ScaleIntensityReducer.reduce_on_intensities(reflection_table)
    return reflection_table

  @classmethod
  def apply_scaling_factors(cls, reflection_table):
    reflection_table = SumAndPrfIntensityReducer.apply_scaling_factors(reflection_table)
    reflection_table = ScaleIntensityReducer.apply_scaling_factors(reflection_table)
    return reflection_table

def sum_partial_reflections(reflection_table):
  """Sum partial reflections; weighted sum for summation integration; weighted
  average for profile fitted reflections. N.B. this will report total
  partiality for the summed reflection."""

  intensities = []
  for intensity in ['prf', 'scale', 'sum']:
    if 'intensity.'+intensity+'.value' in reflection_table:
      intensities.append(intensity)

  isel = (reflection_table['partiality'] < 0.99).iselection()
  if len(isel) == 0:
    return reflection_table

  # create map of partial_id to reflections 
  delete = flex.size_t()
  partial_map = defaultdict(list)
  for j in isel:
    partial_map[reflection_table['partial_id'][j]].append(j)

  # now work through this map - get total partiality for every reflection;
  # here only consider reflections with > 1 component;
  partial_ids = []
  for p_id in partial_map:
    if len(partial_map[p_id]) > 1:
      partial_ids.append(p_id)

  header = ["Partial id", "Partiality"]
  for i in intensities:
    header.extend([str(i)+" intensity", str(i)+" variance"])
  rows = []

  # Now loop through 'matched' partials, summing and then deleting before return
  for p_id in partial_ids:
    j = partial_map[p_id]
    for i in j:
      data = [str(p_id), str(reflection_table['partiality'][i])]
      for intensity in intensities:
        data.extend([str(reflection_table['intensity.'+intensity+'.value'][i]),
          str(reflection_table['intensity.'+intensity+'.variance'][i])])
      rows.append(data)

    # do the summing of the partiality values separately to allow looping
    # over multiple times 
    total_partiality = sum([reflection_table['partiality'][i] for i in j])
    if 'prf' in intensities:
      reflection_table = _sum_prf_partials(reflection_table, j)
    if 'sum' in intensities:
      reflection_table = _sum_sum_partials(reflection_table, j)
    if 'scale' in intensities:
      reflection_table = _sum_scale_partials(reflection_table, j)
    # FIXME now that the partials have been summed, should fractioncalc be set
    # to one (except for summation case?)
    reflection_table['partiality'][j[0]] = total_partiality
    delete.extend(flex.size_t(j[1:]))
    data = ['combined '+str(p_id), str(total_partiality)]
    for intensity in intensities:
      data.extend([str(reflection_table['intensity.'+intensity+'.value'][j[0]]),
        str(reflection_table['intensity.'+intensity+'.variance'][j[0]])])
    rows.append(data)
  reflection_table.del_selected(delete)

  logger.debug('\nSummary of combination of partial reflections')
  st = simple_table(rows, header)
  logger.debug(st.format())
  return reflection_table

# FIXME what are the correct weights to use for the different cases? - why
# weighting by (I/sig(I))^2 not just 1/variance for prf. See tests?

def _sum_prf_partials(reflection_table, partials_isel_for_pid):
  """Sum prf partials and set the updated value in the first entry"""
  j = partials_isel_for_pid
  value = reflection_table['intensity.prf.value'][j[0]]
  variance = reflection_table['intensity.prf.variance'][j[0]]
  weight = value * value / variance
  value *= weight
  variance *= weight
  total_weight = weight
  for i in j[1:]:
    _value = reflection_table['intensity.prf.value'][i]
    _variance = reflection_table['intensity.prf.variance'][i]
    _weight = _value * _value / _variance
    value += _weight * _value
    variance += _weight * _variance
    total_weight += _weight
  # now write these back into original reflection
  reflection_table['intensity.prf.value'][j[0]] = value / total_weight
  reflection_table['intensity.prf.variance'][j[0]] = variance / total_weight
  return reflection_table

def _sum_sum_partials(reflection_table, partials_isel_for_pid):
  """Sum sum partials and set the updated value in the first entry"""
  j = partials_isel_for_pid
  value = reflection_table['intensity.sum.value'][j[0]]
  variance = reflection_table['intensity.sum.variance'][j[0]]
  for i in j[1:]:
    value += reflection_table['intensity.sum.value'][i]
    variance += reflection_table['intensity.sum.variance'][i]
  reflection_table['intensity.sum.value'][j[0]] = value
  reflection_table['intensity.sum.variance'][j[0]] = variance
  return reflection_table

def _sum_scale_partials(reflection_table, partials_isel_for_pid):
  """Sum scale partials and set the updated value in the first entry."""
  # Weight scaled intensity partials by 1/variance. See
  # https://en.wikipedia.org/wiki/Weighted_arithmetic_mean, section 
  # 'Dealing with variance'
  j = partials_isel_for_pid
  variance = reflection_table['intensity.scale.variance'][j[0]]
  value = reflection_table['intensity.scale.value'][j[0]] / variance
  total_weight = 1.0/ variance
  for i in j[1:]:
    _variance = reflection_table['intensity.scale.variance'][i]
    value += reflection_table['intensity.scale.value'][i] / _variance
    total_weight += 1.0/_variance
  reflection_table['intensity.scale.value'][j[0]] = value / total_weight
  reflection_table['intensity.scale.variance'][j[0]] = 1.0 / total_weight
  return reflection_table
