from __future__ import annotations

import os


def test(dials_regression):
    import numpy as np

    from iotbx.xds import integrate_hkl, xparm
    from rstbx.cftbx.coordinate_frame_converter import coordinate_frame_converter

    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.util import ioutil

    # The XDS files to read from
    integrate_filename = os.path.join(
        dials_regression, "data", "sim_mx", "INTEGRATE.HKL"
    )
    gxparm_filename = os.path.join(dials_regression, "data", "sim_mx", "GXPARM.XDS")

    # Read the XDS files
    integrate_handle = integrate_hkl.reader()
    integrate_handle.read_file(integrate_filename)
    gxparm_handle = xparm.reader()
    gxparm_handle.read_file(gxparm_filename)

    # Get the parameters we need from the GXPARM file
    d_min = 1.6
    space_group_type = ioutil.get_space_group_type_from_xparm(gxparm_handle)
    cfc = coordinate_frame_converter(gxparm_filename)
    unit_cell = cfc.get_unit_cell()

    # Generate the indices
    index_generator = IndexGenerator(unit_cell, space_group_type, d_min)
    miller_indices = index_generator.to_array()

    # Get individual generated hkl
    gen_h = [hkl[0] for hkl in miller_indices]
    gen_k = [hkl[1] for hkl in miller_indices]
    gen_l = [hkl[2] for hkl in miller_indices]

    # Get individual xds generated hkl
    xds_h = [hkl[0] for hkl in integrate_handle.hkl]
    xds_k = [hkl[1] for hkl in integrate_handle.hkl]
    xds_l = [hkl[2] for hkl in integrate_handle.hkl]

    # Get min/max generated hkl
    min_gen_h, max_gen_h = np.min(gen_h), np.max(gen_h)
    min_gen_k, max_gen_k = np.min(gen_k), np.max(gen_k)
    min_gen_l, max_gen_l = np.min(gen_l), np.max(gen_l)

    # Get min/max xds generated hkl
    min_xds_h, max_xds_h = np.min(xds_h), np.max(xds_h)
    min_xds_k, max_xds_k = np.min(xds_k), np.max(xds_k)
    min_xds_l, max_xds_l = np.min(xds_l), np.max(xds_l)

    # Ensure we have the whole xds range  in the generated set
    assert min_gen_h <= min_xds_h and max_gen_h >= max_xds_h
    assert min_gen_k <= min_xds_k and max_gen_k >= max_xds_k
    assert min_gen_l <= min_xds_l and max_gen_l >= max_xds_l
