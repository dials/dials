# coding: utf-8

from __future__ import absolute_import, division, print_function

import logging
import time
from collections import defaultdict
from math import ceil, cos, floor, log, pi, sin, sqrt

from dials.array_family import flex
from dials.util.version import dials_version
from dials.util.filter_reflections import filter_reflection_table
from iotbx import mtz
from libtbx.utils import Sorry
from scitbx import matrix

try:
  from math import isclose
except ImportError:
  # Python 3 backport
  def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
      return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

logger = logging.getLogger(__name__)

def dials_u_to_mosflm(dials_U, uc):
  '''Compute the mosflm U matrix i.e. the U matrix from same UB definition
  as DIALS, but with Busing & Levy B matrix definition.'''

  parameters = uc.parameters()
  dials_B = matrix.sqr(uc.fractionalization_matrix()).transpose()
  dials_UB = dials_U * dials_B

  r_parameters = uc.reciprocal_parameters()

  a = parameters[:3]
  al = [pi * p / 180.0 for p in parameters[3:]]
  b = r_parameters[:3]
  be = [pi * p / 180.0 for p in r_parameters[3:]]

  mosflm_B = matrix.sqr((b[0], b[1] * cos(be[2]), b[2] * cos(be[1]),
                         0, b[1] * sin(be[2]), - b[2] * sin(be[1]) * cos(al[0]),
                         0, 0, 1.0 / a[2]))

  mosflm_U = dials_UB * mosflm_B.inverse()

  return mosflm_U


def _add_batch(mtz, experiment, batch_number, image_number, force_static_model):
  """Add a single image's metadata to an mtz file.

  Returns the batch object.
  """
  assert batch_number > 0

  # Recalculate useful numbers and references here
  wavelength = experiment.beam.get_wavelength()
  # We ignore panels beyond the first one, at the moment
  panel = experiment.detector[0]

  if experiment.goniometer:
    axis = matrix.col(experiment.goniometer.get_rotation_axis())
  else:
    axis = 0.0, 0.0, 0.0

  U = matrix.sqr(experiment.crystal.get_U())
  if experiment.goniometer is not None:
    F = matrix.sqr(experiment.goniometer.get_fixed_rotation())
  else:
    F = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

  # Create the batch object and start configuring it
  o = mtz.add_batch().set_num(batch_number).set_nbsetid(1).set_ncryst(1)
  o.set_time1(0.0).set_time2(0.0).set_title('Batch {}'.format(batch_number))
  o.set_ndet(1).set_theta(flex.float((0.0, 0.0))).set_lbmflg(0)
  o.set_alambd(wavelength).set_delamb(0.0).set_delcor(0.0)
  o.set_divhd(0.0).set_divvd(0.0)

  # FIXME hard-coded assumption on indealized beam vector below... this may be
  # broken when we come to process data from a non-imgCIF frame
  s0n = matrix.col(experiment.beam.get_s0()).normalize().elems
  o.set_so(flex.float(s0n)).set_source(flex.float((0, 0, -1)))

  # these are probably 0, 1 respectively, also flags for how many are set, sd
  o.set_bbfac(0.0).set_bscale(1.0)
  o.set_sdbfac(0.0).set_sdbscale(0.0).set_nbscal(0)

  # unit cell (this is fine) and the what-was-refined-flags FIXME hardcoded

  # take time-varying parameters from the *end of the frame* unlikely to
  # be much different at the end - however only exist if scan-varying
  # refinement was used
  if not force_static_model and experiment.crystal.num_scan_points > 0:
    # Get the index of the image in the sequence e.g. first => 0, second => 1
    image_index = image_number - experiment.image_range[0]
    _unit_cell = experiment.crystal.get_unit_cell_at_scan_point(image_index)
    _U = matrix.sqr(experiment.crystal.get_U_at_scan_point(image_index))
  else:
    _unit_cell = experiment.crystal.get_unit_cell()
    _U = U

  # apply the fixed rotation to this to unify matrix definitions - F * U
  # was what was used in the actual prediction: U appears to be stored
  # as the transpose?! At least is for Mosflm...
  #
  # FIXME Do we need to apply the setting rotation here somehow? i.e. we have
  # the U.B. matrix assuming that the axis is equal to S * axis_datum but
  # here we are just giving the effective axis so at scan angle 0 this will
  # not be correct... FIXME 2 not even sure we can express the stack of
  # matrices S * R * F * U * B in MTZ format?... see [=A=] below
  _U = dials_u_to_mosflm(F * _U, _unit_cell)

  # FIXME need to get what was refined and what was constrained from the
  # crystal model - see https://github.com/dials/dials/issues/355
  o.set_cell(flex.float(_unit_cell.parameters()))
  o.set_lbcell(flex.int((-1, -1, -1, -1, -1, -1)))
  o.set_umat(flex.float(_U.transpose().elems))

  # get the mosaic spread though today it may not actually be set - should
  # this be in the BATCH headers?
  try:
    mosaic = experiment.crystal.get_mosaicity()
  except AttributeError:
    mosaic = 0
  o.set_crydat(flex.float([mosaic, 0.0, 0.0, 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

  o.set_lcrflg(0)
  o.set_datum(flex.float((0.0, 0.0, 0.0)))

  # detector size, distance
  o.set_detlm(flex.float([0.0, panel.get_image_size()[0],
                          0.0, panel.get_image_size()[1],
                          0, 0, 0, 0]))
  o.set_dx(flex.float([panel.get_directed_distance(), 0.0]))

  # goniometer axes and names, and scan axis number, and num axes, missets
  # [=A=] should we be using this to unroll the setting matrix etc?
  o.set_e1(flex.float(axis))
  o.set_e2(flex.float((0.0, 0.0, 0.0)))
  o.set_e3(flex.float((0.0, 0.0, 0.0)))
  o.set_gonlab(flex.std_string(('AXIS', '', '')))
  o.set_jsaxs(1)
  o.set_ngonax(1)
  o.set_phixyz(flex.float((0.0, 0.0, 0.0, 0.0, 0.0, 0.0)))

  # scan ranges, axis
  if experiment.scan:
    phi_start, phi_range = experiment.scan.get_image_oscillation(image_number)
  else:
    phi_start, phi_range = 0.0, 0.0

  o.set_phistt(phi_start)
  o.set_phirange(phi_range)
  o.set_phiend(phi_start + phi_range)
  o.set_scanax(flex.float(axis))

  # number of misorientation angles
  o.set_misflg(0)

  # crystal axis closest to rotation axis (why do I want this?)
  o.set_jumpax(0)

  # type of data - 1; 2D, 2; 3D, 3; Laue
  o.set_ldtype(2)

  return o


def _write_columns(mtz_file, dataset, integrated_data):
  """Write the column definitions AND data for a single dataset."""

  # now create the actual data structures - first keep a track of the columns

  # H K L M/ISYM BATCH I SIGI IPR SIGIPR FRACTIONCALC XDET YDET ROT WIDTH
  # LP MPART FLAG BGPKRATIOS

  # gather the required information for the reflection file

  nref = len(integrated_data['miller_index'])

  # check reflections remain
  if nref == 0:
    raise Sorry('no reflections for export')

  xdet, ydet, zdet = [flex.double(x) for x in integrated_data['xyzobs.px.value'].parts()]

  # compute BATCH values - floor() to get (fortran) image captured within
  #                        +1     because FORTRAN counting; zdet+1=image_index
  #                        +off   because            image_index+o=batch
  batch = (flex.floor(zdet).iround() + 1) + integrated_data["batch_offset"]

  # now add column information...

  # FIXME add DIALS_FLAG which can include e.g. was partial etc.

  type_table = {
    'H': 'H',
    'K': 'H',
    'L': 'H',
    'I': 'J',
    'SIGI': 'Q',
    'IPR': 'J',
    'SIGIPR': 'Q',
    'BG' : 'R',
    'SIGBG' : 'R',
    'XDET': 'R',
    'YDET': 'R',
    'BATCH': 'B',
    'BGPKRATIOS': 'R',
    'WIDTH': 'R',
    'MPART': 'I',
    'M_ISYM': 'Y',
    'FLAG': 'I',
    'LP': 'R',
    'FRACTIONCALC': 'R',
    'ROT': 'R',
    'QE': 'R',
  }

  # derive index columns from original indices with
  #
  # from m.replace_original_index_miller_indices
  #
  # so all that is needed now is to make space for the reflections - fill with
  # zeros...

  mtz_file.adjust_column_array_sizes(nref)
  mtz_file.set_n_reflections(nref)

  # assign H, K, L, M_ISYM space
  for column in 'H', 'K', 'L', 'M_ISYM':
    dataset.add_column(column, type_table[column]).set_values(
      flex.double(nref, 0.0).as_float())

  mtz_file.replace_original_index_miller_indices(integrated_data['miller_index_rebase'])

  dataset.add_column('BATCH', type_table['BATCH']).set_values(
    batch.as_double().as_float())

  #if intensity values used in scaling exist, then just export these as I, SIGI
  if 'intensity.scale.value' in integrated_data:
    I_scaling = integrated_data['intensity.scale.value']
    V_scaling = integrated_data['intensity.scale.variance']
    # Trap negative variances
    assert V_scaling.all_gt(0)
    dataset.add_column('I', type_table['I']).set_values(I_scaling.as_float())
    dataset.add_column('SIGI', type_table['SIGI']).set_values(
      flex.sqrt(V_scaling).as_float())
    dataset.add_column('SCALEUSED', 'R').set_values(
      integrated_data['inverse_scale_factor'].as_float())
    dataset.add_column('SIGSCALEUSED', 'R').set_values(
      flex.sqrt(integrated_data['inverse_scale_factor_variance']).as_float())
  else:
    if 'intensity.prf.value' in integrated_data:
      if 'intensity.sum.value' in integrated_data:
        col_names = ('IPR', 'SIGIPR')
      else:
        col_names = ('I', 'SIGI')
      I_profile = integrated_data['intensity.prf.value']
      V_profile = integrated_data['intensity.prf.variance']
      # Trap negative variances
      assert V_profile.all_gt(0)
      dataset.add_column(col_names[0], type_table['I']).set_values(
        I_profile.as_float())
      dataset.add_column(col_names[1], type_table['SIGI']).set_values(
        flex.sqrt(V_profile).as_float())
    if 'intensity.sum.value' in integrated_data:
      I_sum = integrated_data['intensity.sum.value']
      V_sum = integrated_data['intensity.sum.variance']
      # Trap negative variances
      assert V_sum.all_gt(0)
      dataset.add_column('I', type_table['I']).set_values(I_sum.as_float())
      dataset.add_column('SIGI', type_table['SIGI']).set_values(
        flex.sqrt(V_sum).as_float())
  if ('background.sum.value' in integrated_data and
      'background.sum.variance' in integrated_data):
    bg = integrated_data['background.sum.value']
    varbg = integrated_data['background.sum.variance']
    assert (varbg >= 0).count(False) == 0
    sigbg = flex.sqrt(varbg)
    dataset.add_column('BG', type_table['BG']).set_values(bg.as_float())
    dataset.add_column('SIGBG', type_table['SIGBG']).set_values(sigbg.as_float())

  dataset.add_column('FRACTIONCALC', type_table['FRACTIONCALC']).set_values(
    integrated_data['fractioncalc'].as_float())

  dataset.add_column('XDET', type_table['XDET']).set_values(xdet.as_float())
  dataset.add_column('YDET', type_table['YDET']).set_values(ydet.as_float())
  dataset.add_column('ROT', type_table['ROT']).set_values(integrated_data["ROT"].as_float())
  if 'lp' in integrated_data:
    dataset.add_column('LP', type_table['LP']).set_values(integrated_data['lp'].as_float())
  if 'qe' in integrated_data:
    dataset.add_column('QE', type_table['QE']).set_values(integrated_data['qe'].as_float())
  elif 'dqe' in integrated_data:
    dataset.add_column('QE', type_table['QE']).set_values(integrated_data['dqe'].as_float())
  else:
    dataset.add_column('QE', type_table['QE']).set_values(flex.double(nref, 1.0).as_float())

def _next_epoch(val):
  """Find a reasonably round epoch a small number above an existing one.

  Examples: 130-138     => 140
            139         => 150
            1234        => 1300
            19999-20998 => 21000
  """

  # Find the order of magnitude-1 (minimum: 1 as want no fractional values)
  small_magnitude = 10**max(1, int(floor(log(val, 10))-1))
  # How many units of this we have (float cast for __division__ insensitivity)
  mag_multiple = int(ceil(val / float(small_magnitude)))
  epoch = small_magnitude * mag_multiple
  # If this would give a consecutive number then offset it by a magnitude step
  if epoch <= val + 1:
    epoch = small_magnitude * (mag_multiple+1)
  return epoch


def _calculate_batch_offsets(experiments):
  """Take a list of experiments and resolve and return the batch offsets.

  This is the number added to the image number to give the
  batch number, such that:
  - Each experiment has a unique, nonoverlapping, nonconsecutive range
  - None are zero
  - Image number ranges are kept if at all possible
  """

  experiments_to_shift = []
  existing_ranges = set()
  maximum_batch_number = 0
  batch_offsets = [0]*len(experiments)

  # Handle zeroth shifts and kept ranges
  for i, experiment in enumerate(experiments):
    ilow, ihigh = experiment.image_range
    # Check assumptions
    assert ilow <= ihigh, "Inverted image order!?"
    assert ilow >= 0, "Negative image indices are not expected"
    # Don't emit zero: Causes problems with C/fortran number conversion
    if ilow == 0:
      ilow, ihigh = ilow+1, ihigh+1
    # If we overlap with anything, then process later
    if any( ilow <= high+1 and ihigh >= low-1 for low, high in existing_ranges):
      experiments_to_shift.append((i, experiment))
    else:
      batch_offsets[i] = ilow-experiment.image_range[0]
      existing_ranges.add((ilow, ihigh))
      maximum_batch_number = max(maximum_batch_number, ihigh)

  # Now handle all the experiments that overlapped by pushing them higher
  for i, experiment in experiments_to_shift:
    start_number = _next_epoch(maximum_batch_number)
    range_width = experiment.image_range[1]-experiment.image_range[0]+1
    end_number = start_number + range_width - 1
    batch_offsets[i] = start_number - experiment.image_range[0]
    maximum_batch_number = end_number
    experiment.scan.set_batch_offset(batch_offsets[i])

  return batch_offsets


def export_mtz(integrated_data, experiment_list, params):
  '''Export data from integrated_data corresponding to experiment_list to an
  MTZ file hklout.'''

  # Convert experiment_list to a real python list or else identity assumptions
  # fail like:
  #   assert experiment_list[0] is experiment_list[0]
  # And assumptions about added attributes break
  experiment_list = list(experiment_list)

  # Validate multi-experiment assumptions
  if len(experiment_list) > 1:
    # All experiments should match crystals, or else we need multiple crystals/datasets
    if not all(x.crystal == experiment_list[0].crystal for x in experiment_list[1:]):
      logger.warning("Warning: Experiment crystals differ. Using first experiment crystal for file-level data.")

    # We must match wavelengths (until multiple datasets supported)
    if not all(isclose(x.beam.get_wavelength(), experiment_list[0].beam.get_wavelength(), rel_tol=1e-4) for x in experiment_list[1:]):
      data = [x.beam.get_wavelength() for x in experiment_list]
      raise Sorry("Cannot export multiple experiments with different beam wavelengths ({})".format(data))

  # also only work correctly with one panel (for the moment)
  if any(len(experiment.detector) != 1 for experiment in experiment_list):
    logger.warning('Warning: Ignoring multiple panels in output MTZ')

  # Clean up the data with the passed in options
  integrated_data = filter_reflection_table(integrated_data,
    intensity_choice=params.intensity, 
    partiality_threshold=params.mtz.partiality_threshold,
    combine_partials=params.mtz.combine_partials,
    min_isigi=params.mtz.min_isigi, filter_ice_rings=params.mtz.filter_ice_rings,
    d_min=params.mtz.d_min)

  # Calculate and store the image range for each image
  for experiment in experiment_list:
    # Calculate this once so that we don't have to again
    if experiment.scan:
      experiment.image_range = experiment.scan.get_image_range()
    else:
      experiment.image_range = 1, 1


  batch_offsets = flex.int(
    expt.scan.get_batch_offset() for expt in experiment_list)
  if batch_offsets.all_eq(0):
    # Calculate any offset to the image numbers
    batch_offsets = _calculate_batch_offsets(experiment_list)
  else:
    unique_offsets = set(batch_offsets)
    if len(unique_offsets) != len(batch_offsets):
      import collections
      raise Sorry('Duplicate batch offsets detected: %s' %', '.join(
        str(item) for item, count in collections.Counter(batch_offsets).items()
        if count > 1))

  # Create the mtz file
  mtz_file = mtz.object()
  mtz_file.set_title('from dials.export_mtz')
  date_str = time.strftime('%Y-%m-%d at %H:%M:%S %Z')
  if time.strftime('%Z') != 'GMT':
    date_str += time.strftime('  (%Y-%m-%d at %H:%M:%S %Z)', time.gmtime())
  mtz_file.add_history('From %s, run on %s' % (dials_version(), date_str))

  # FIXME TODO for more than one experiment into an MTZ file:
  #
  # - add an epoch (or recover an epoch) from the scan and add this as an extra
  #   column to the MTZ file for scaling, so we know that the two lattices were
  #   integrated at the same time
  # ✓ decide a sensible BATCH increment to apply to the BATCH value between
  #   experiments and add this

  for experiment_index, experiment in enumerate(experiment_list):
    # Grab our subset of the data
    experiment.data = dict(integrated_data.select(integrated_data["id"] == experiment_index))

    # Do any crystal transformations for the experiment
    cb_op_to_ref = experiment.crystal.get_space_group().info(
        ).change_of_basis_op_to_reference_setting()
    experiment.crystal = experiment.crystal.change_basis(cb_op_to_ref)
    experiment.data["miller_index_rebase"] = cb_op_to_ref.apply(experiment.data["miller_index"])

    s0 = experiment.beam.get_s0()
    s0n = matrix.col(s0).normalize().elems
    logger.debug('Beam vector: %.4f %.4f %.4f' % s0n)

    for i in range(experiment.image_range[0], experiment.image_range[1]+1):
      _add_batch(mtz_file, experiment,
        batch_number=i+experiment.scan.get_batch_offset(),
        image_number=i,
        force_static_model=params.mtz.force_static_model)

    # Create the batch offset array. This gives us an experiment (id)-dependent
    # batch offset to calculate the correct batch from image number.
    experiment.data["batch_offset"] = flex.int(len(experiment.data["id"]),
      experiment.scan.get_batch_offset())

    # Calculate whether we have a ROT value for this experiment, and set the column
    _, _, z = experiment.data['xyzcal.px'].parts()
    if experiment.scan:
      experiment.data["ROT"] = experiment.scan.get_angle_from_array_index(z)
    else:
      experiment.data["ROT"] = z

  # Update the mtz general information now we've processed the experiments
  mtz_file.set_space_group_info(experiment_list[0].crystal.get_space_group().info())
  unit_cell = experiment_list[0].crystal.get_unit_cell()
  mtz_crystal = mtz_file.add_crystal('XTAL', 'DIALS', unit_cell.parameters())
  mtz_dataset = mtz_crystal.add_dataset('FROMDIALS', experiment_list[0].beam.get_wavelength())

  # Combine all of the experiment data columns before writing
  merged_data = {k: v.deep_copy() for k, v in experiment_list[0].data.items()}
  for experiment in experiment_list[1:]:
    for k, v in experiment.data.items():
      merged_data[k].extend(v)
  # ALL columns must be the same length
  assert len(set(len(v) for v in merged_data.values())) == 1, "Column length mismatch"
  assert len(merged_data["id"]) == len(integrated_data["id"]), "Lost rows in split/combine"

  # Write all the data and columns to the mtz file
  _write_columns(mtz_file, mtz_dataset, merged_data)

  logger.info("Saving {} integrated reflections to {}".format(
    len(merged_data['id']), params.mtz.hklout))
  mtz_file.write(params.mtz.hklout)

  return mtz_file
