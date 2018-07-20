"""
Define classes to handle filtering and reduction of reflection table data from
integrated/scaled pickle files to columns for export.
Intended to replace the 'apply data filters' function of export_mtz.

FilteringReductionMethod defines general filtering methods on reflection tables
as staticmethods, which can be used elsewhere to exporting.

FilterForExportAlgorithm defines a filter_for_export method/algorithm, and
defines the abstract methods which must be defined for individual cases.

Four cases are defined as separate classes for clarity, as there are subtleties
to each case. The main differences are reduce_on_intensities and
apply_scaling_factors. The filter_for_export method calls the namesake in
FilterForExportAlgorithm, while allowing for different parameters to be defined
for each case or even replacement of the algorithm for one case if required.
Methods are implemented as class methods to allow general filtering/reduction
of reflection tables for given intensity types.

Phil parameters needed to operate:
 - filter_ice_rings=False
 - min_isigi=None
 - combine_partials=True; flag to switch on/off partial combining
 - partiality_threshold=0.99; partials are summed (and scaled if summation int)
    and if the total is below the partiality threshold then
 - intensities: *prf, sum, scale - choice, allow choice of both prf+sum, choice of
    which intensities will be exported. scale is for exporting after scaling,
    prf/sum for exporting after integration.
"""
import logging
import abc
from collections import defaultdict
from libtbx.utils import Sorry
from dials.array_family import flex

logger = logging.getLogger('dials')

def filter_for_export(reflection_table, intensity_choice=['scale'],
  *args, **kwargs):
  if intensity_choice == ['scale']:
    return ScaleIntensityReducer.filter_for_export(
      reflection_table, *args, **kwargs)
  elif intensity_choice == ['sum']:
    return SumIntensityReducer.filter_for_export(
      reflection_table, *args, **kwargs)
  elif intensity_choice == ['prf']:
    return PrfIntensityReducer.filter_for_export(
      reflection_table, *args, **kwargs)
  elif intensity_choice == ['sum', 'prf'] or intensity_choice == ['prf', 'sum']:
    return SumAndPrfIntensityReducer.filter_for_export(
      reflection_table, *args, **kwargs)
  else:
    raise Sorry("Unrecognised intensity choice for filter_for_export")

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
      if combine_partials:
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
  def _filter_for_export(cls, reflection_table, min_isigi=None,
    filter_ice_rings=False, combine_partials=True, partiality_threshold=0.99):
    """Designed to be called by subclasses"""
    assert reflection_table.size() > 0, \
      """Empty reflection table given to reduce_data_for_export function"""

    for intensity in cls.allowed_intensities:
      if intensity not in cls.intensities:
        if 'intensity.'+intensity+'.value' in reflection_table:
          del reflection_table['intensity.'+intensity+'.value']
          del reflection_table['intensity.'+intensity+'.variance']

    reflection_table = cls.filter_unassigned_reflections(reflection_table)
    reflection_table = cls.reduce_on_intensities(reflection_table)

    if reflection_table.size() == 0:
      raise Sorry('No suitable reflections found for intensity choice: %s' % \
        ' '.join(['intensity.'+i+'.value' for i in cls.intensities]))

    reflection_table = cls.filter_bad_variances(reflection_table)
    if filter_ice_rings:
      reflection_table = cls.filter_ice_rings(reflection_table)

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
    filter_ice_rings=False, combine_partials=True, partiality_threshold=0.99):
    return cls._filter_for_export(reflection_table, min_isigi, filter_ice_rings,
      combine_partials, partiality_threshold)

  @classmethod
  def filter_on_min_isigi(cls, reflection_table, min_isigi=None):
    reflections = cls._filter_on_min_isigi(reflection_table,
      cls.intensities[0], min_isigi)
    return reflections

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
    filter_ice_rings=False, combine_partials=True, partiality_threshold=0.99):
    return cls._filter_for_export(reflection_table, min_isigi, filter_ice_rings,
      combine_partials, partiality_threshold)

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
  implement filtering for export. Reflections are kept if either a prf or
  sum intensity is defined."""

  intensities = ['sum', 'prf']

  @classmethod
  def filter_for_export(cls, reflection_table, min_isigi=None, filter_ice_rings=False,
      combine_partials=True, partiality_threshold=0.99):
    return cls._filter_for_export(
      reflection_table, min_isigi=None, filter_ice_rings=False,
      combine_partials=True, partiality_threshold=0.99)

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
  # First select the reflections which have successfully been integrated by
  # either method
    selection = reflection_table.get_flags(reflection_table.flags.integrated,
      all=False)
    reflection_table = reflection_table.select(selection)
    logger.info("Selected %d integrated reflections" % reflection_table.size())
    # Want to ensure that sensible intensities are in both columns. So if there is
    # no value for prf or sum for a given reflection, set the value to be equal
    # to the one that does exist. In effect the prf intensity column then becomes
    # "prf intensity if prf fitting successful else summation" and vice versa.
    profile_selection = reflection_table.get_flags(
      reflection_table.flags.integrated_prf)
    if profile_selection.count(True) < reflection_table.size():
      sum_for_prf_sel = ~profile_selection
      sum_int_for_prf = reflection_table['intensity.sum.value'].select(sum_for_prf_sel)
      reflection_table['intensity.prf.value'].set_selected(
        sum_for_prf_sel.iselection(), sum_int_for_prf)
      sum_var_for_prf = reflection_table['intensity.sum.variance'].select(sum_for_prf_sel)
      reflection_table['intensity.prf.variance'].set_selected(
        sum_for_prf_sel.iselection(), sum_var_for_prf)
    summation_selection = reflection_table.get_flags(
    reflection_table.flags.integrated_sum)
    if summation_selection.count(True) < reflection_table.size():
      prf_for_sum_sel = ~summation_selection
      prf_int_for_sum = reflection_table['intensity.prf.value'].select(prf_for_sum_sel)
      reflection_table['intensity.sum.value'].set_selected(
        prf_for_sum_sel.iselection(), prf_int_for_sum)
      prf_var_for_sum = reflection_table['intensity.prf.variance'].select(prf_for_sum_sel)
      reflection_table['intensity.sum.variance'].set_selected(
        prf_for_sum_sel.iselection(), prf_var_for_sum)
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
    filter_ice_rings=False, combine_partials=True, partiality_threshold=0.99):
    return cls._filter_for_export(reflection_table, min_isigi=None,
      filter_ice_rings=False, combine_partials=True, partiality_threshold=0.99)

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
    reflection_table = reflection_table.select(
      reflection_table['inverse_scale_factor'] > 0.0)
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

  # Now loop through 'matched' partials, summing and then deleting before return
  for p_id in partial_ids:
    j = partial_map[p_id]
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
  reflection_table.del_selected(delete)
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
