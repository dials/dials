#!/usr/bin/env python
#
# export_cif.py
#
#  Copyright (C) 2017 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
import logging
logger = logging.getLogger(__name__)

class CIFOutputFile(object):
  '''
  Class to output experiments and reflections as CIF file

  '''

  def __init__(self, filename):
    '''
    Init with the filename

    '''
    self.filename = filename
    self.handle = open(filename, "wb")

  def write(self, experiments, reflections):
    '''
    Write the experiments and reflections to file

    '''
    assert len(experiments) == 1, "Only 1 experiment is handled at the moment"

    # Select reflections
    selection = reflections.get_flags(reflections.flags.integrated, all=True)
    reflections = reflections.select(selection)

    # Write out data block relationships
    # FIXME not sure what correct thing to put in for <dataset name> and <parent dataset name>
    #   data_block_relationships
    #   loop
    #   <dataset name> <parent dataset name>
    self.handle.write("data_block_relationships\n")
    self.handle.write("loop\n")
    self.handle.write("dataset dataset\n")

    # <datablock_name>
    self.handle.write("data_block\n")

    # Hard coding X-ray
    self.handle.write("_pdbx_data_section.type_scattering x-ray\n")
    self.handle.write("_pdbx_data_section.type_merged false\n")
    self.handle.write("_pdbx_data_section.type_scaled false\n")

    # Write the crystal information
    #   <crystal id> <a> <b> <c> <alpha> <beta> <gamma> <wavelength>
    self.handle.write("loop\n")
    crystal = experiments[0].crystal
    wavelength = experiments[0].beam.get_wavelength()
    unit_cell_parameters = {}
    if crystal.num_scan_points > 1:
      for i in range(crystal.num_scan_points):
        a, b, c, alpha, beta, gamma = crystal.get_unit_cell_at_scan_point(i).parameters()
        unit_cell_parameters[i] = (a, b, c, alpha, beta, gamma)
        self.handle.write("%d %f %f %f %f %f %f %f\n" % (
          0, a, b, c, alpha, beta, gamma, wavelength))
    else:
      a, b, c, alpha, beta, gamma = crystal.get_unit_cell().parameters()
      unit_cell_parameters[0] = (a, b, c, alpha, beta, gamma)
      self.handle.write("%d %f %f %f %f %f %f %f\n" % (
        0, a, b, c, alpha, beta, gamma, wavelength))

    # Write the image data
    #  <image id> <image number> <crystal id> <a> <b> <c> <alpha> <beta> <gamma> <phi-image>
    scan = experiments[0].scan
    z0 = scan.get_image_range()[0]
    self.handle.write("loop\n")
    for i in range(len(scan)):
      z = z0 + i
      if crystal.num_scan_points > 1:
        crystal_id = i
      else:
        crystal_id = 0
      a, b, c, alpha, beta, gamma = unit_cell_parameters[crystal_id]
      phi = scan.get_angle_from_image_index(z)
      self.handle.write("%d %d %d %f %f %f %f %f %f %f\n" % (
        i,
        z,
        crystal_id,
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
        phi))

    # Write reflection data
    # FIXME there are three intensiry fields. I've put summation in I and Isum
    #   <reflection id>
    #   <image id start>
    #   <image id end>
    #   <h-original>
    #   <k-original>
    #   <l-original>
    #   <I>
    #   <sigI>
    #   <I-sum>
    #   <sigI-sum>
    #   <I-profile>
    #   <sigI-profile>
    #   <phi-reflection>
    #   <partiality>
    self.handle.write("loop\n")
    for i, r in enumerate(reflections):
      _,_,_,_,z0,z1 = r['bbox']
      h, k, l       = r['miller_index']
      I             = r['intensity.sum.value']
      sigI          = r['intensity.sum.variance']
      Isum          = r['intensity.sum.value']
      sigIsum       = r['intensity.sum.variance']
      Iprf          = r['intensity.prf.value']
      sigIprf       = r['intensity.prf.variance']
      phi           = r['xyzcal.mm'][2]
      partiality    = r['partiality']
      self.handle.write("%d %d %d %d %d %d %f %f %f %f %f %f %f %f\n" % (
        i,
        z0,
        z1,
        h,
        k,
        l,
        I,
        sigI,
        Isum,
        sigIsum,
        Iprf,
        sigIprf,
        phi,
        partiality))

    logger.info("Wrote reflections to %s" % self.filename)
