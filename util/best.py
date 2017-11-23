#!/usr/bin/env python
#
# dials.util.best.py
#
#  Copyright (C) 2016 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import absolute_import, division

import logging
logger = logging.getLogger(__name__)

def write_background_file(file_name, imageset, n_bins):
  from dials.command_line.background import background
  d, I, sig = background(imageset, imageset.indices()[0], n_bins=n_bins)

  logger.info('Saving background file to %s' %file_name)
  with open(file_name, 'wb') as f:
    for d_, I_, sig_ in zip(d, I, sig):
      print >> f, '%10.4f %10.2f %10.2f' %(d_, I_, sig_)

def write_integrated_hkl(prefix, reflections):
  from scitbx.array_family import flex
  expt_ids = reflections['id']
  integrated_sel = reflections.get_flags(reflections.flags.integrated_sum)
  for i_expt in range(flex.max(expt_ids)+1):
    integrated = reflections.select(
      (reflections['id'] == i_expt) & integrated_sel)
    integrated.sort("miller_index")
    h,k,l = integrated['miller_index'].as_vec3_double().parts()
    I = integrated['intensity.sum.value']
    sigI = flex.sqrt(integrated['intensity.sum.variance'])

    suffix = ''
    if flex.max(expt_ids) > 0:
      suffix = '%i' %(i_expt+1)
    file_name = '%s%s.hkl' %(prefix, suffix)
    logger.info('Saving reflections to %s' %file_name)
    with open(file_name, 'wb') as f:
      for i in range(len(integrated)):
        print >> f, '%4.0f %4.0f %4.0f %10.2f %10.2f' %(h[i], k[i], l[i], I[i], sigI[i])

def write_par_file(file_name, experiment):
  from scitbx import matrix
  from dxtbx.model import Crystal
  from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
  from iotbx.mtz.extract_from_symmetry_lib import ccp4_symbol

  imageset = experiment.imageset
  detector = imageset.get_detector()
  goniometer = imageset.get_goniometer()
  beam = imageset.get_beam()
  scan = imageset.get_scan()

  R_to_mosflm = align_reference_frame(
    beam.get_s0(), (1.0, 0.0, 0.0),
    goniometer.get_rotation_axis(), (0.0, 0.0, 1.0))

  cryst = experiment.crystal
  cryst = cryst.change_basis(
    cryst.get_space_group().info()\
      .change_of_basis_op_to_reference_setting())
  A = matrix.sqr(cryst.get_A())
  A_inv = A.inverse()

  real_space_a = R_to_mosflm * A_inv.elems[:3]
  real_space_b = R_to_mosflm * A_inv.elems[3:6]
  real_space_c = R_to_mosflm * A_inv.elems[6:9]

  cryst_mosflm = Crystal(
    real_space_a, real_space_b, real_space_c,
    space_group=cryst.get_space_group())
  A_mosflm = matrix.sqr(cryst_mosflm.get_A())
  U_mosflm = matrix.sqr(cryst_mosflm.get_U())
  B_mosflm = matrix.sqr(cryst_mosflm.get_B())
  UB_mosflm = U_mosflm * B_mosflm
  uc_params = cryst_mosflm.get_unit_cell().parameters()
  assert U_mosflm.is_r3_rotation_matrix(), U_mosflm

  symmetry = cryst_mosflm.get_space_group().type().number()
  beam_centre = tuple(reversed(detector[0].get_beam_centre(beam.get_s0())))
  distance = detector[0].get_directed_distance()
  polarization = R_to_mosflm * matrix.col(beam.get_polarization_normal())
  rotation = matrix.col(goniometer.get_rotation_axis())
  if (rotation.angle(matrix.col(detector[0].get_fast_axis())) <
      rotation.angle(matrix.col(detector[0].get_slow_axis()))):
    direction = 'FAST'
  else:
    direction = 'SLOW'
  rotation = R_to_mosflm * rotation

  # Calculate average spot diameter for SEPARATION parameter
  # http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_parameters.html
  # BEAM_DIVERGENCE=
  # This value is approximately arctan(spot diameter/DETECTOR_DISTANCE)
  import math
  profile = experiment.profile
  spot_diameter = math.tan(profile.delta_b() * math.pi/180) * distance
  spot_diameter_px = spot_diameter * detector[0].get_pixel_size()[0]

  # determine parameters for RASTER keyword
  # http://www.mrc-lmb.cam.ac.uk/harry/cgi-bin/keyword2.cgi?RASTER

  # NXS, NYS (odd integers) define the overall dimensions of the rectangular array of pixels for each spot
  # NXS and NYS are set to twice the spot size plus 5 pixels
  nxs = 2 * int(math.ceil(spot_diameter_px)) + 5
  nys = nxs

  # NRX, NRY are the number of columns or rows of points in the background rim
  # NRX and NRY are set to half the spot size plus 2 pixels
  nrx = int(math.ceil(0.5 * spot_diameter_px)) + 2
  nry = nrx

  # NC the corner background cut-off which corresponds to a half-square of side NC points
  # NC is set to the mean of the spot size in X and Y plus 4
  nc = int(math.ceil(spot_diameter_px)) + 4

  def space_group_symbol(space_group):
    symbol = ccp4_symbol(space_group.info(), lib_name='syminfo.lib',
                         require_at_least_one_lib=False)
    if symbol != 'P 1':
      symbol = symbol.replace(' 1', '')
    symbol = symbol.replace(' ', '')
    return symbol

  logger.info('Saving BEST parameter file to %s' %file_name)
  with open(file_name, 'wb') as f:#
    print >> f, '# parameter file for BEST'
    print >> f, 'TITLE          From DIALS'
    print >> f, 'DETECTOR       PILA'
    print >> f, 'SITE           Not set'
    print >> f, 'DIAMETER       %6.2f' %(max(detector[0].get_image_size()) * detector[0].get_pixel_size()[0])
    print >> f, 'PIXEL          %s' %detector[0].get_pixel_size()[0]
    print >> f, 'ROTAXIS        %4.2f %4.2f %4.2f' %rotation.elems, direction
    print >> f, 'POLAXIS        %4.2f %4.2f %4.2f' %polarization.elems
    print >> f, 'GAIN               1.00' # correct for Pilatus images
    # http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/FAQ#You_said_that_the_XDS_deals_with_high_mosaicity._How_high_mosaicity_is_still_manageable.3F
    # http://journals.iucr.org/d/issues/2012/01/00/wd5161/index.html
    # Transform from XDS defintion of sigma_m to FWHM (MOSFLM mosaicity definition)
    print >> f, 'CMOSAIC            %.2f' %(experiment.profile.sigma_m() * 2.355)
    print >> f, 'PHISTART           %.2f' %scan.get_oscillation_range()[0]
    print >> f, 'PHIWIDTH           %.2f' %scan.get_oscillation()[1]
    print >> f, 'DISTANCE        %7.2f' %distance
    print >> f, 'WAVELENGTH      %.5f' %beam.get_wavelength()
    print >> f, 'POLARISATION    %7.5f' %beam.get_polarization_fraction()
    print >> f, 'SYMMETRY       %s' %space_group_symbol(cryst.get_space_group())
    print >> f, 'UB             %9.6f %9.6f %9.6f' %UB_mosflm[:3]
    print >> f, '               %9.6f %9.6f %9.6f' %UB_mosflm[3:6]
    print >> f, '               %9.6f %9.6f %9.6f' %UB_mosflm[6:]
    print >> f, 'CELL           %8.2f %8.2f %8.2f %6.2f %6.2f %6.2f' %uc_params
    print >> f, 'RASTER           %i %i %i %i %i' %(nxs, nys, nc, nrx, nry)
    print >> f, 'SEPARATION      %.3f  %.3f' %(spot_diameter, spot_diameter)
    print >> f, 'BEAM           %8.3f %8.3f' %beam_centre
    print >> f, '# end of parameter file for BEST'

