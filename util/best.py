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
from __future__ import division

def write_background_file(file_name, imageset, n_bins):
  from dials.command_line.background import background
  d, I, sig = background(imageset, imageset.indices()[0], n_bins=n_bins)

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
    h,k,l = integrated['miller_index'].as_vec3_double().parts()
    I = integrated['intensity.sum.value']
    sigI = flex.sqrt(integrated['intensity.sum.variance'])

    suffix = ''
    if flex.max(expt_ids) > 0:
      suffix = '%i' %(i_expt+1)
    with open('%s%s.hkl' %(prefix, suffix), 'wb') as f:
      for i in range(len(integrated)):
        print >> f, '%4.0f %4.0f %4.0f %10.2f %10.2f' %(h[i], k[i], l[i], I[i], sigI[i])

def write_par_file(file_name, experiment):
  from scitbx import matrix
  from dxtbx.model.crystal import crystal_model
  from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
  from dials.command_line.refine_bravais_settings import short_space_group_name

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
  A = cryst.get_A()
  A_inv = A.inverse()

  real_space_a = R_to_mosflm * A_inv.elems[:3]
  real_space_b = R_to_mosflm * A_inv.elems[3:6]
  real_space_c = R_to_mosflm * A_inv.elems[6:9]

  cryst_mosflm = crystal_model(
    real_space_a, real_space_b, real_space_c,
    space_group=cryst.get_space_group(),
    mosaicity=cryst.get_mosaicity())
  A_mosflm = cryst_mosflm.get_A()
  U_mosflm = cryst_mosflm.get_U()
  B_mosflm = cryst_mosflm.get_B()
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
    print >> f, 'CMOSAIC            %.2f' %experiment.profile.sigma_m()
    print >> f, 'PHISTART           %.2f' %scan.get_oscillation_range()[0]
    print >> f, 'PHIWIDTH           %.2f' %scan.get_oscillation()[1]
    print >> f, 'DISTANCE        %7.2f' %distance
    print >> f, 'WAVELENGTH      %.5f' %beam.get_wavelength()
    print >> f, 'POLARISATION    %7.5f' %beam.get_polarization_fraction()
    print >> f, 'SYMMETRY       %s' %short_space_group_name(cryst.get_space_group())
    print >> f, 'UB             %9.6f %9.6f %9.6f' %UB_mosflm[:3]
    print >> f, '               %9.6f %9.6f %9.6f' %UB_mosflm[3:6]
    print >> f, '               %9.6f %9.6f %9.6f' %UB_mosflm[6:]
    print >> f, 'CELL           %8.2f %8.2f %8.2f %6.2f %6.2f %6.2f' %uc_params
    print >> f, 'RASTER           13  13   7   3   4'
    print >> f, 'SEPARATION      2.960  2.960'
    print >> f, 'BEAM           %8.3f %8.3f' %beam_centre
    print >> f, '# end of parameter file for BEST'

