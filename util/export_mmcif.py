#!/usr/bin/env python
#
# export_mmcif.py
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

class MMCIFOutputFile(object):
  '''
  Class to output experiments and reflections as MMCIF file

  '''

  def __init__(self, filename):
    '''
    Init with the filename

    '''
    import iotbx.cif.model
    self._cif = iotbx.cif.model.cif()
    self.filename = filename

  def write(self, experiments, reflections):
    '''
    Write the experiments and reflections to file

    '''
    import iotbx.cif.model
    assert len(experiments) == 1, "Only 1 experiment is handled at the moment"

    # Select reflections
    selection = reflections.get_flags(reflections.flags.integrated, all=True)
    reflections = reflections.select(selection)

    # Get the cif block
    cif_block = iotbx.cif.model.block()

    # Write out data block relationships
    # FIXME not sure what correct thing to put in for <dataset name> and <parent dataset name>
    #   data_block_relationships
    #   loop
    #   <dataset name> <parent dataset name>
    # FIXME where to input "data_block_relationships"?
    cif_loop = iotbx.cif.model.loop(header=("_dataset_name",
                                            "_parent_dataset_name"))
    cif_loop.add_row(("dataset", "None"))
    cif_block.add_loop(cif_loop)

    # FIXME not sure where to put data_block_<NAME>

    # Hard coding X-ray
    cif_block["_pdbx_data_section.type_scattering"] = "x-ray"
    cif_block["_pdbx_data_section.type_merged"] = False
    cif_block["_pdbx_data_section.type_scaled"] = False

    # Write the crystal information
    #   <crystal id> <a> <b> <c> <alpha> <beta> <gamma> <wavelength>
    cif_loop = iotbx.cif.model.loop(header=("_crystal_id", "_a", "_b", "_c",
                                            "_alpha", "_beta", "_gamma",
                                            "_wavelength"))
    crystal = experiments[0].crystal
    wavelength = experiments[0].beam.get_wavelength()
    unit_cell_parameters = {}
    if crystal.num_scan_points > 1:
      for i in range(crystal.num_scan_points):
        a, b, c, alpha, beta, gamma = crystal.get_unit_cell_at_scan_point(i).parameters()
        unit_cell_parameters[i] = (a, b, c, alpha, beta, gamma)
        cif_loop.add_row((i, a, b, c, alpha, beta, gamma, wavelength))
    else:
      a, b, c, alpha, beta, gamma = crystal.get_unit_cell().parameters()
      unit_cell_parameters[0] = (a, b, c, alpha, beta, gamma)
      cif_loop.add_row((0, a, b, c, alpha, beta, gamma, wavelength))
    #cif_block.add_loop(cif_loop)


    # Write the image data
    #  <image id> <image number> <crystal id> <a> <b> <c> <alpha> <beta> <gamma> <phi-image>
    scan = experiments[0].scan
    z0 = scan.get_image_range()[0]
    cif_loop = iotbx.cif.model.loop(header=("_image_id", "_image_number",
                                        "_crystal_id", "_a", "_b", "_c",
                                        "_alpha", "_beta", "_gamma",
                                        "_phi_image"))
    for i in range(len(scan)):
      z = z0 + i
      if crystal.num_scan_points > 1:
        crystal_id = i
      else:
        crystal_id = 0
      a, b, c, alpha, beta, gamma = unit_cell_parameters[crystal_id]
      phi = scan.get_angle_from_image_index(z)
      cif_loop.add_row((i, z, crystal_id, a, b, c, alpha, beta, gamma, phi))
    cif_block.add_loop(cif_loop)

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
    cif_loop = iotbx.cif.model.loop(header=("_reflection_id", "_image_id_start",
                                            "_image_id_end", "_h_original",
                                            "_k_original", "_l_original", "_I",
                                            "_sigI", "_I_sum", "_sigI_sum",
                                            "_I_profile", "_sigI_profile",
                                            "_phi_reflection", "_partiality"))
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
      cif_loop.add_row((i, z0, z1, h, k, l, I, sigI, Isum, sigIsum, Iprf, sigIprf, phi, partiality))
    cif_block.add_loop(cif_loop)

    # Add the block
    self._cif['dials'] = cif_block

    # Print to file
    print >>open(self.filename, "w"), self._cif

    # Log
    logger.info("Wrote reflections to %s" % self.filename)
