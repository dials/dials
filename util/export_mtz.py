from __future__ import absolute_import, division

import logging
logger = logging.getLogger(__name__)

def sum_partial_reflections(integrated_data, min_total_partiality=0.5):
  '''Sum partial reflections; weighted sum for summation integration; weighted
  average for profile fitted reflections. N.B. this will report total
  partiality for the summed reflection.'''

  if not 'partiality' in integrated_data:
    return integrated_data

  from dials.array_family import flex
  from collections import defaultdict

  # rest of the columns are pretty well defined - these are all uniform for the
  # reflection so can just ignore other values:
  #
  # d
  # id
  # entering
  # flags
  # lp
  # miller_index
  # panel
  # partial_id
  # s1
  # xyzcal.mm
  # xyzcal.px
  # zeta
  #
  # Then finally these need a weighted average:
  #
  # intensity.prf.value
  # intensity.prf.variance
  # intensity.sum.value
  # intensity.sum.variance
  # partiality
  # bbox - this is not used...
  #
  # now just need to worry about those that I am actually outputting to the MTZ
  # file...

  isel = (integrated_data['partiality'] < 0.99).iselection()

  if len(isel) == 0:
    return integrated_data

  delete = flex.size_t()
  partial_map = defaultdict(list)

  # create map of partial_id to reflections j

  for j in isel:
    partial_map[integrated_data['partial_id'][j]].append(j)

  # now work through this map - get total partiality for every reflection;
  # here only consider reflections with > 1 component; if total partiality
  # less than min_total_partiality discard all parts.

  partial_ids = []

  for p_id in partial_map:
    if len(partial_map[p_id]) > 1:
      partial_ids.append(p_id)

  # work through multipart partials; compute those weighted values I need
  # if total partiality less than min, delete. if summing, delete extra parts

  we_got_profiles = 'intensity.prf.value' in integrated_data
  logger.info('Profile fitted reflections: %s' % we_got_profiles)

  for p_id in partial_ids:
    p_tot = sum([integrated_data['partiality'][j] for j in partial_map[p_id]])
    if p_tot < min_total_partiality:
      for j in partial_map[p_id]:
        delete.append(j)
      continue

    j0 = partial_map[p_id][0]
    jrest = partial_map[p_id][1:]

    # FIXME revisiting this calculation am not sure it is correct - why
    # weighting by (I/sig(I))^2 not just 1/variance?
    if we_got_profiles:
      prf_value = integrated_data['intensity.prf.value'][j0]
      prf_variance = integrated_data['intensity.prf.variance'][j0]
      weight = prf_value * prf_value / prf_variance
      prf_value *= weight
      prf_variance *= weight
      total_weight = weight
    sum_value = integrated_data['intensity.sum.value'][j0]
    sum_variance = integrated_data['intensity.sum.variance'][j0]
    partiality = integrated_data['partiality'][j0]

    # weight profile fitted intensity and variance computed from weights
    # proportional to (I/sig(I))^2; schedule for deletion spare parts

    for j in jrest:
      delete.append(j)

      sum_value += integrated_data['intensity.sum.value'][j]
      sum_variance += integrated_data['intensity.sum.variance'][j]
      partiality += integrated_data['partiality'][j]

      if we_got_profiles:
        _prf_value = integrated_data['intensity.prf.value'][j]
        _prf_variance = integrated_data['intensity.prf.variance'][j]

        _weight = _prf_value * _prf_value / _prf_variance
        prf_value += _weight * _prf_value
        prf_variance += _weight * _prf_variance
        total_weight += _weight

    # now write these back into original reflection
    if we_got_profiles:
      prf_value /= total_weight
      prf_variance /= total_weight
      integrated_data['intensity.prf.value'][j0] = prf_value
      integrated_data['intensity.prf.variance'][j0] = prf_variance
    integrated_data['intensity.sum.value'][j0] = sum_value
    integrated_data['intensity.sum.variance'][j0] = sum_variance
    integrated_data['partiality'][j0] = partiality

  integrated_data.del_selected(delete)

  return integrated_data

def scale_partial_reflections(integrated_data, min_partiality=0.5):
  '''Scale partial reflections (after summation) according to their estimated
  partiality - for profile fitted reflections this will result in no change,
  for summation integrated reflections will be scaled up by 1 / partiality
  with error accordingly scaled. N.B. this will report the scaled up partiality
  for the output reflection.'''

  # assert: in here there will be no multi-part partial reflections

  if not 'partiality' in integrated_data:
    return integrated_data

  isel = (integrated_data['partiality'] < 1.0).iselection()

  if len(isel) == 0:
    return integrated_data

  from dials.array_family import flex
  delete = flex.size_t()

  for j in isel:
    if integrated_data['partiality'][j] < min_partiality:
      delete.append(j)
      continue
    inv_p = 1.0 / integrated_data['partiality'][j]
    integrated_data['intensity.sum.value'][j] *= inv_p
    integrated_data['intensity.sum.variance'][j] *= inv_p
    integrated_data['partiality'][j] *= 1.0

  integrated_data.del_selected(delete)

  return integrated_data

def dials_u_to_mosflm(dials_U, uc):
  '''Compute the mosflm U matrix i.e. the U matrix from same UB definition
  as DIALS, but with Busing & Levy B matrix definition.'''

  from scitbx.matrix import sqr
  from math import sin, cos, pi

  parameters = uc.parameters()
  dials_B = sqr(uc.fractionalization_matrix()).transpose()
  dials_UB = dials_U * dials_B

  r_parameters = uc.reciprocal_parameters()

  a = parameters[:3]
  al = [pi * p / 180.0 for p in parameters[3:]]
  b = r_parameters[:3]
  be = [pi * p / 180.0 for p in r_parameters[3:]]

  mosflm_B = sqr((b[0], b[1] * cos(be[2]), b[2] * cos(be[1]),
                  0, b[1] * sin(be[2]), - b[2] * sin(be[1]) * cos(al[0]),
                  0, 0, 1.0 / a[2]))

  mosflm_U = dials_UB * mosflm_B.inverse()

  return mosflm_U

def export_mtz(integrated_data, experiment_list, hklout, ignore_panels=False,
               include_partials=False, keep_partials=False, min_isigi=None,
               force_static_model=False, filter_ice_rings=False,
               ignore_profile_fitting=False):
  '''Export data from integrated_data corresponding to experiment_list to an
  MTZ file hklout.'''

  from dials.array_family import flex
  from scitbx import matrix

  # for the moment assume (and assert) that we will convert data from exactly
  # one lattice...

  # FIXME allow for more than one experiment in here: this is fine just add
  # multiple MTZ data sets (DIALS1...DIALSN) and multiple batch headers: one
  # range of batches for each experiment

  assert(len(experiment_list) == 1)
  # select reflections that are assigned to an experiment (i.e. non-negative id)
  integrated_data = integrated_data.select(integrated_data['id'] >= 0)
  assert max(integrated_data['id']) == 0
  logger.info('Read %s predicted reflections' % len(integrated_data))

  # Ignore profile fitted
  if ignore_profile_fitting:
    del integrated_data['intensity.prf.value']
    del integrated_data['intensity.prf.variance']

  if 'intensity.prf.variance' in integrated_data:
    selection = integrated_data.get_flags(
      integrated_data.flags.integrated,
      all=True)
  else:
    selection = integrated_data.get_flags(
      integrated_data.flags.integrated_sum)
  integrated_data = integrated_data.select(selection)
  logger.info("Selected %d integrated reflections" % len(integrated_data))

  # check we have reflections left - see #357
  if len(integrated_data) == 0:
    from libtbx.utils import Sorry
    raise Sorry("All reflections excluded based on flags.integrated")

  selection = integrated_data['intensity.sum.variance'] <= 0
  if selection.count(True) > 0:
    integrated_data.del_selected(selection)
    logger.info('Removing %d reflections with negative variance' % \
          selection.count(True))

  if 'intensity.prf.variance' in integrated_data:
    selection = integrated_data['intensity.prf.variance'] <= 0
    if selection.count(True) > 0:
      integrated_data.del_selected(selection)
      logger.info('Removing %d profile reflections with negative variance' % \
            selection.count(True))

  if filter_ice_rings:
    selection = integrated_data.get_flags(integrated_data.flags.in_powder_ring)
    integrated_data.del_selected(selection)
    logger.info("Removing %d reflections in ice ring resolutions" %
                selection.count(True))

  if min_isigi is not None:

    selection = (
      integrated_data['intensity.sum.value']/
      flex.sqrt(integrated_data['intensity.sum.variance'])) < min_isigi
    integrated_data.del_selected(selection)
    logger.info('Removing %d reflections with I/Sig(I) < %s' %(
      selection.count(True), min_isigi))

    if 'intensity.prf.variance' in integrated_data:
      selection = (
        integrated_data['intensity.prf.value']/
        flex.sqrt(integrated_data['intensity.prf.variance'])) < min_isigi
      integrated_data.del_selected(selection)
      logger.info('Removing %d profile reflections with I/Sig(I) < %s' %(
        selection.count(True), min_isigi))

  # FIXME in here work on including partial reflections => at this stage best
  # to split off the partial refections into a different selection & handle
  # gracefully... better to work on a short list as will need to "pop" them &
  # find matching parts to combine.

  if include_partials:
    integrated_data = sum_partial_reflections(integrated_data)
    integrated_data = scale_partial_reflections(integrated_data)

  if 'partiality' in integrated_data:
    selection = integrated_data['partiality'] < 0.99
    if selection.count(True) > 0 and not keep_partials:
      integrated_data.del_selected(selection)
      logger.info('Removing %d incomplete reflections' % \
        selection.count(True))

  # FIXME TODO for more than one experiment into an MTZ file:
  #
  # - add an epoch (or recover an epoch) from the scan and add this as an extra
  #   column to the MTZ file for scaling, so we know that the two lattices were
  #   integrated at the same time
  # - decide a sensible BATCH increment to apply to the BATCH value between
  #   experiments and add this
  #
  # At the moment this is probably enough to be working on.

  experiment = experiment_list[0]

  # also only work with one panel(for the moment)

  if not ignore_panels:
    if len(experiment.detector) != 1:
      from libtbx.utils import Sorry
      raise Sorry('Only a single panel will be considered. To ignore panels '
                  'other than the first for a multi-panel detector, set '
                  'mtz.ignore_panels=true')

  from scitbx import matrix

  if experiment.goniometer:
    axis = matrix.col(experiment.goniometer.get_rotation_axis())
  else:
    axis = 0.0, 0.0, 0.0
  s0 = experiment.beam.get_s0()
  wavelength = experiment.beam.get_wavelength()

  panel = experiment.detector[0]
  origin = matrix.col(panel.get_origin())
  fast = matrix.col(panel.get_fast_axis())
  slow = matrix.col(panel.get_slow_axis())

  pixel_size = panel.get_pixel_size()

  fast *= pixel_size[0]
  slow *= pixel_size[1]

  cb_op_to_ref = experiment.crystal.get_space_group().info(
    ).change_of_basis_op_to_reference_setting()

  experiment.crystal = experiment.crystal.change_basis(cb_op_to_ref)

  U = matrix.sqr(experiment.crystal.get_U())
  if experiment.goniometer is not None:
    F = matrix.sqr(experiment.goniometer.get_fixed_rotation())
  else:
    F = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))
  unit_cell = experiment.crystal.get_unit_cell()

  from iotbx import mtz

  from scitbx.array_family import flex
  from math import floor, sqrt

  m = mtz.object()
  m.set_title('from dials.export_mtz')
  m.set_space_group_info(experiment.crystal.get_space_group().info())

  from dials.util.version import dials_version
  import time
  date_str = time.strftime('%d/%m/%Y at %H:%M:%S', time.gmtime())
  m.add_history('From %s, run on %s' % (dials_version(), date_str))

  if experiment.scan:
    image_range = experiment.scan.get_image_range()
  else:
    image_range = 1, 1

  # pointless (at least) doesn't like batches starting from zero
  b_incr = max(image_range[0], 1)

  s0n = matrix.col(s0).normalize().elems

  logger.info('Beam vector: %.4f %.4f %.4f' % s0n)

  for b in range(image_range[0], image_range[1] + 1):
    o = m.add_batch().set_num(b+b_incr).set_nbsetid(1).set_ncryst(1)
    o.set_time1(0.0).set_time2(0.0).set_title('Batch %d' % (b+b_incr))
    o.set_ndet(1).set_theta(flex.float((0.0, 0.0))).set_lbmflg(0)
    o.set_alambd(wavelength).set_delamb(0.0).set_delcor(0.0)
    o.set_divhd(0.0).set_divvd(0.0)

    # FIXME hard-coded assumption on indealized beam vector below... this may be
    # broken when we come to process data from a non-imgCIF frame
    o.set_so(flex.float(s0n)).set_source(flex.float((0, 0, -1)))

    # these are probably 0, 1 respectively, also flags for how many are set, sd
    o.set_bbfac(0.0).set_bscale(1.0)
    o.set_sdbfac(0.0).set_sdbscale(0.0).set_nbscal(0)

    # unit cell (this is fine) and the what-was-refined-flags FIXME hardcoded

    # take time-varying parameters from the *end of the frame* unlikely to
    # be much different at the end - however only exist if time-varying refinement
    # was used
    if not force_static_model and experiment.crystal.num_scan_points > 0:
      _unit_cell = experiment.crystal.get_unit_cell_at_scan_point(b-image_range[0])
      _U = matrix.sqr(experiment.crystal.get_U_at_scan_point(b-image_range[0]))
    else:
      _unit_cell = unit_cell
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
    mosaic = experiment.crystal.get_mosaicity()
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
      phi_start, phi_range = experiment.scan.get_image_oscillation(b)
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

  # now create the actual data structures - first keep a track of the columns
  # H K L M/ISYM BATCH I SIGI IPR SIGIPR FRACTIONCALC XDET YDET ROT WIDTH
  # LP MPART FLAG BGPKRATIOS

  from cctbx.array_family import flex as cflex # implicit import
  from cctbx.miller import map_to_asu_isym # implicit import

  # gather the required information for the reflection file

  nref = len(integrated_data['miller_index'])
  x_px, y_px, z_px = integrated_data['xyzcal.px'].parts()

  xdet = flex.double(x_px)
  ydet = flex.double(y_px)
  zdet = flex.double(z_px)

  # compute ROT values
  if experiment.scan:
    rot = flex.double([experiment.scan.get_angle_from_image_index(z)
                      for z in zdet])
  else:
    rot = zdet

  # compute BATCH values
  batch = flex.floor(zdet).iround() + 1 + b_incr

  # we're working with full reflections so...
  fractioncalc = flex.double(nref, 1.0)

  # now go for it and make an MTZ file...

  x = m.add_crystal('XTAL', 'DIALS', unit_cell.parameters())
  d = x.add_dataset('FROMDIALS', wavelength)

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
    'DQE': 'R',
  }

  # derive index columns from original indices with
  #
  # from m.replace_original_index_miller_indices
  #
  # so all that is needed now is to make space for the reflections - fill with
  # zeros...

  m.adjust_column_array_sizes(nref)
  m.set_n_reflections(nref)

  # assign H, K, L, M_ISYM space
  for column in 'H', 'K', 'L', 'M_ISYM':
    d.add_column(column, type_table[column]).set_values(
      flex.double(nref, 0.0).as_float())

  m.replace_original_index_miller_indices(cb_op_to_ref.apply(
    integrated_data['miller_index']))

  d.add_column('BATCH', type_table['BATCH']).set_values(
    batch.as_double().as_float())

  if 'lp' in integrated_data:
    lp = integrated_data['lp']
  else:
    lp = flex.double(nref, 1.0)
  if 'dqe' in integrated_data:
    dqe = integrated_data['dqe']
  else:
    dqe = flex.double(nref, 1.0)
  I_profile = None
  V_profile = None
  I_sum = None
  V_sum = None
  # FIXME errors in e.g. LP correction need to be propogated here
  scl = lp / dqe
  if 'intensity.prf.value' in integrated_data:
    I_profile = integrated_data['intensity.prf.value'] * scl
    V_profile = integrated_data['intensity.prf.variance'] * scl * scl
    # Trap negative variances
    assert V_profile.all_gt(0)
    d.add_column('IPR', type_table['I']).set_values(I_profile.as_float())
    d.add_column('SIGIPR', type_table['SIGI']).set_values(
      flex.sqrt(V_profile).as_float())
  if 'intensity.sum.value' in integrated_data:
    I_sum = integrated_data['intensity.sum.value'] * scl
    V_sum = integrated_data['intensity.sum.variance'] * scl * scl
    # Trap negative variances
    assert V_sum.all_gt(0)
    d.add_column('I', type_table['I']).set_values(I_sum.as_float())
    d.add_column('SIGI', type_table['SIGI']).set_values(
      flex.sqrt(V_sum).as_float())
  if ('background.sum.value' in integrated_data and
      'background.sum.variance' in integrated_data):
    bg = integrated_data['background.sum.value']
    varbg = integrated_data['background.sum.variance']
    assert (varbg >= 0).count(False) == 0
    sigbg = flex.sqrt(varbg)
    d.add_column('BG', type_table['BG']).set_values(bg.as_float())
    d.add_column('SIGBG', type_table['SIGBG']).set_values(sigbg.as_float())

  d.add_column('FRACTIONCALC', type_table['FRACTIONCALC']).set_values(
    fractioncalc.as_float())

  d.add_column('XDET', type_table['XDET']).set_values(xdet.as_float())
  d.add_column('YDET', type_table['YDET']).set_values(ydet.as_float())
  d.add_column('ROT', type_table['ROT']).set_values(rot.as_float())
  d.add_column('LP', type_table['LP']).set_values(lp.as_float())
  d.add_column('DQE', type_table['DQE']).set_values(dqe.as_float())

  m.write(hklout)

  return m
