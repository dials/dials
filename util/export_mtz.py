def export_mtz(integrated_data, experiment_list, hklout):
  '''Export data from integrated_data corresponding to experiment_list to an
  MTZ file hklout.'''

  # for the moment assume (and assert) that we will convert data from exactly
  # one lattice...

  assert(len(experiment_list) == 1)
  assert(min(integrated_data['id']) == max(integrated_data['id']) == 0)

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

  # FIXME this may want to be time-dependent in the future...
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

    # FIXME hard-coded assumption on indealized beam vector below...
    o.set_so(flex.float(s0)).set_source(flex.float((0, 0, 1)))

    # these are probably 0, 1 respectively, also flags for how many are set, sd
    o.set_bbfac(0.0).set_bscale(1.0)
    o.set_sdbfac(0.0).set_sdbscale(0.0).set_nbscal(0)

    # unit cell (this is fine) and the what-was-refined-flags FIXME hardcoded
    o.set_cell(flex.float(unit_cell.parameters()))
    o.set_lbcell(flex.int((-1, -1, -1, 0, 0, 0)))
    o.set_umat(flex.float(U.transpose().elems))

    # sadly we don't record the mosaic spread at the moment (by design)
    mosaic = 0.0
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

  lp = integrated_data['intensity.cor.value'] / \
    integrated_data['intensity.raw.value']

  # compute ROT values
  rot = flex.double([experiment.scan.get_angle_from_image_index(z) for z in zdet])

  # compute BATCH values
  batch = flex.floor(zdet).iround() + 1

  # compute M/ISYM
  mi = integrated_data['miller_index']
  h, k, l = zip(*integrated_data['miller_index'])
  isym = flex.int(len(integrated_data['miller_index']))
  misym = flex.int(len(isym))
  map_to_asu_isym(experiment.crystal.get_space_group().type(), True,
                  integrated_data['miller_index'], isym)

  for j, i in enumerate(isym):
    if i < 0:
      mis = (-i + 1) * 2
    else:
      mis = (i + 1) * 2 - 1
    misym[j] = mis

  # we're working with full reflections so...
  fractioncalc = flex.double(len(isym), 1.0)

  # now go for it and make an MTZ file...

  x = m.add_crystal('XTAL', 'DIALS', unit_cell.parameters())
  d = x.add_dataset('FROMDIALS', wavelength)

  # now add column information...

  type_table = {'IPR': 'J', 'BGPKRATIOS': 'R', 'WIDTH': 'R', 'I': 'J',
                'H': 'H', 'K': 'H', 'MPART': 'I', 'L': 'H', 'BATCH': 'B',
                'M_ISYM': 'Y', 'SIGI': 'Q', 'FLAG': 'I', 'XDET': 'R', 'LP': 'R',
                'YDET': 'R', 'SIGIPR': 'Q', 'FRACTIONCALC': 'R', 'ROT': 'R'}

  h, k, l = zip(*integrated_data['miller_index'])
  h, k, l = flex.int(h), flex.int(k), flex.int(l)

  m.adjust_column_array_sizes(len(h))
  m.set_n_reflections(len(h))

  # assign H, K, L
  d.add_column('H', type_table['H']).set_values(h.as_double().as_float())
  d.add_column('K', type_table['K']).set_values(k.as_double().as_float())
  d.add_column('L', type_table['L']).set_values(l.as_double().as_float())

  d.add_column('M_ISYM', type_table['M_ISYM']).set_values(
    misym.as_double().as_float())
  d.add_column('BATCH', type_table['BATCH']).set_values(
    batch.as_double().as_float())

  d.add_column('I', type_table['I']).set_values(
    integrated_data['intensity.cor.value'].as_float())

  # FIXME properly handle negative variance estimates... like don't make them
  # in the first place...
  d.add_column('SIGI', type_table['SIGI']).set_values(
    flex.sqrt(flex.fabs(integrated_data['intensity.cor.value'])).as_float())

  d.add_column('FRACTIONCALC', type_table['FRACTIONCALC']).set_values(
    fractioncalc.as_float())

  d.add_column('XDET', type_table['XDET']).set_values(xdet.as_float())
  d.add_column('YDET', type_table['YDET']).set_values(ydet.as_float())
  d.add_column('ROT', type_table['ROT']).set_values(rot.as_float())
  d.add_column('LP', type_table['LP']).set_values(lp.as_float())

  m.write(hklout)

  return m
