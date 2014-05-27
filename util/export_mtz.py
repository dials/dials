def export_mtz(integrated_data, experiment_list, hklout):
  '''Export data from integrated_data corresponding to experiment_list to an
  MTZ file hklout.'''

  # for the moment assume (and assert) that we will convert data from exactly
  # one lattice...

  # FIXME test to see if the model is scan varying, and if it is (which it can be)
  # write out a different U matrix for each batch... This should be easy. N.B.
  # should use U matrix from middle of batch (probably)

  # FIXME allow for more than one experiment in here: this is fine just add
  # multiple MTZ data sets (DIALS1...DIALSN) and multiple batch headers: one
  # range of batches for each experiment

  assert(len(experiment_list) == 1)
  assert(min(integrated_data['id']) == max(integrated_data['id']) == 0)

  # strip out negative variance reflections: these should not really be there
  # FIXME Doing select on summation results. Should do on profile result if
  # present?
  selection = integrated_data['intensity.sum.variance'] < 0
  if selection.count(True) > 0:
    integrated_data.del_selected(selection)
    print 'Removing %d reflections with negative variance' % selection.count(True)

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

  assert(len(experiment.detector) == 1)

  axis = experiment.goniometer.get_rotation_axis()
  s0 = experiment.beam.get_s0()
  wavelength = experiment.beam.get_wavelength()

  from scitbx import matrix

  panel = experiment.detector[0]
  origin = matrix.col(panel.get_origin())
  fast = matrix.col(panel.get_fast_axis())
  slow = matrix.col(panel.get_slow_axis())

  pixel_size = panel.get_pixel_size()

  fast *= pixel_size[0]
  slow *= pixel_size[1]

  U = experiment.crystal.get_U()
  unit_cell = experiment.crystal.get_unit_cell()
  from iotbx import mtz

  from scitbx.array_family import flex
  from math import floor, sin, sqrt

  m = mtz.object()
  m.set_title('from dials.export_mtz')
  m.set_space_group_info(experiment.crystal.get_space_group().info())

  image_range = experiment.scan.get_image_range()

  for b in range(image_range[0], image_range[1] + 1):
    o = m.add_batch().set_num(b).set_nbsetid(1).set_ncryst(1)
    o.set_time1(0.0).set_time2(0.0).set_title('Batch %d' % b)
    o.set_ndet(1).set_theta(flex.float((0.0, 0.0))).set_lbmflg(0)
    o.set_alambd(wavelength).set_delamb(0.0).set_delcor(0.0)
    o.set_divhd(0.0).set_divvd(0.0)

    # FIXME hard-coded assumption on indealized beam vector below... this may be
    # broken when we come to process data from a non-imgCIF frame
    o.set_so(flex.float(s0)).set_source(flex.float((0, 0, -1)))

    # these are probably 0, 1 respectively, also flags for how many are set, sd
    o.set_bbfac(0.0).set_bscale(1.0)
    o.set_sdbfac(0.0).set_sdbscale(0.0).set_nbscal(0)

    # unit cell (this is fine) and the what-was-refined-flags FIXME hardcoded

    # take time-varying parameters from the *end of the frame* unlikely to
    # be much different at the end - however only exist if time-varying refinement
    # was used
    if experiment.crystal.num_scan_points:
      _unit_cell = experiment.crystal.get_unit_cell_at_scan_point(b)
      _U = experiment.crystal.get_U_at_scan_point(b)
    else:
      _unit_cell = unit_cell
      _U = U

    # FIXME need to get what was refined and what was constrained from the
    # crystal model
    o.set_cell(flex.float(_unit_cell.parameters()))
    o.set_lbcell(flex.int((-1, -1, -1, -1, -1, -1)))
    o.set_umat(flex.float(_U.transpose().elems))

    # get the mosaic spread though today it may not actually be set
    mosaic = experiment.crystal.get_mosaicity()
    o.set_crydat(flex.float([mosaic, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0]))

    o.set_lcrflg(0)
    o.set_datum(flex.float((0.0, 0.0, 0.0)))

    # detector size, distance
    o.set_detlm(flex.float([0.0, panel.get_image_size()[0],
                            0.0, panel.get_image_size()[1],
                            0, 0, 0, 0]))
    o.set_dx(flex.float([panel.get_distance(), 0.0]))

    # goniometer axes and names, and scan axis number, and number of axes, missets
    o.set_e1(flex.float(axis))
    o.set_e2(flex.float((0.0, 0.0, 0.0)))
    o.set_e3(flex.float((0.0, 0.0, 0.0)))
    o.set_gonlab(flex.std_string(('AXIS', '', '')))
    o.set_jsaxs(1)
    o.set_ngonax(1)
    o.set_phixyz(flex.float((0.0, 0.0, 0.0, 0.0, 0.0, 0.0)))

    # scan ranges, axis
    phi_start, phi_range = experiment.scan.get_image_oscillation(b)
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

  from cctbx.array_family import flex as cflex
  from cctbx.miller import map_to_asu_isym

  # gather the required information for the reflection file

  x_px, y_px, z_px = zip(*integrated_data['xyzcal.px'])

  xdet = flex.double(x_px)
  ydet = flex.double(y_px)
  zdet = flex.double(z_px)

  # compute ROT values
  rot = flex.double([experiment.scan.get_angle_from_image_index(z) for z in zdet])

  # compute BATCH values
  batch = flex.floor(zdet).iround() + 1

  # compute M/ISYM
  # mi = integrated_data['miller_index']
  # h, k, l = zip(*integrated_data['miller_index'])
  # isym = flex.int(len(integrated_data['miller_index']))
  # misym = flex.int(len(isym))
  # map_to_asu_isym(experiment.crystal.get_space_group().type(), True,
  #                 integrated_data['miller_index'], isym)

  # for j, i in enumerate(isym):
  #   if i < 0:
  #     mis = (-i + 1) * 2
  #   else:
  #     mis = (i + 1) * 2 - 1
  #   misym[j] = mis

  # we're working with full reflections so...
  fractioncalc = flex.double(len(rot), 1.0)

  # now go for it and make an MTZ file...

  x = m.add_crystal('XTAL', 'DIALS', unit_cell.parameters())
  d = x.add_dataset('FROMDIALS', wavelength)

  # now add column information...

  type_table = {'IPR': 'J', 'BGPKRATIOS': 'R', 'WIDTH': 'R', 'I': 'J',
                'H': 'H', 'K': 'H', 'MPART': 'I', 'L': 'H', 'BATCH': 'B',
                'M_ISYM': 'Y', 'SIGI': 'Q', 'FLAG': 'I', 'XDET': 'R', 'LP': 'R',
                'YDET': 'R', 'SIGIPR': 'Q', 'FRACTIONCALC': 'R', 'ROT': 'R'}

  # OK these are probably the original indices not the reduced indices that the
  # MTZ file should have - now use Richard's new method to add the H, K, L 
  # columns and M_ISYM N.B. now no longer mapping Miller indices to ASU above
  # so these HKL are currently nonsense values...
  # h, k, l = zip(*integrated_data['miller_index'])
  # h, k, l = flex.int(h), flex.int(k), flex.int(l)

  nref = len(integrated_data['miller_index'])
  m.adjust_column_array_sizes(nref)
  m.set_n_reflections(nref)

  # assign H, K, L space
  # d.add_column('H', type_table['H']).set_values(h.as_double().as_float())
  # d.add_column('K', type_table['K']).set_values(k.as_double().as_float())
  # d.add_column('L', type_table['L']).set_values(l.as_double().as_float())
  for column in 'H', 'K', 'L':
    d.add_column(column, type_table[column]).set_values(
      flex.double(nref, 0.0).as_float())
  
  # make M_ISYM space
  d.add_column('M_ISYM', type_table['M_ISYM']).set_values(
    flex.double(nref, 0.0).as_float())
  
  m.replace_original_index_miller_indices(integrated_data['miller_index'])

  #d.add_column('M_ISYM', type_table['M_ISYM']).set_values(
  #  misym.as_double().as_float())
  d.add_column('BATCH', type_table['BATCH']).set_values(
    batch.as_double().as_float())

  lp = integrated_data['lp']
  if 'intensity.prf.value' in integrated_data:
    I = integrated_data['intensity.prf.value'] * lp
    V = integrated_data['intensity.prf.variance'] * lp
  else:
    I = integrated_data['intensity.sum.value'] * lp
    V = integrated_data['intensity.sum.variance'] * lp

  d.add_column('I', type_table['I']).set_values(I.as_float())

  # Trap negative variances

  assert ((V < 0).count(True) == 0)

  d.add_column('SIGI', type_table['SIGI']).set_values(flex.sqrt(V).as_float())

  d.add_column('FRACTIONCALC', type_table['FRACTIONCALC']).set_values(
    fractioncalc.as_float())

  d.add_column('XDET', type_table['XDET']).set_values(xdet.as_float())
  d.add_column('YDET', type_table['YDET']).set_values(ydet.as_float())
  d.add_column('ROT', type_table['ROT']).set_values(rot.as_float())

  m.write(hklout)

  return m
