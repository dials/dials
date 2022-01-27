from __future__ import annotations

import logging
import math

import iotbx.phil
from cctbx import sgtbx
from cctbx.crystal import symmetry as crystal_symmetry
from cctbx.miller import set as miller_set
from cctbx.sgtbx import space_group as sgtbx_space_group
from libtbx.utils import format_float_with_standard_uncertainty

from dials.algorithms.symmetry import origin
from dials.array_family import flex
from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.check_indexing_symmetry")

help_message = """

This program can be used to analyse the correlation coefficients between
reflections related by the symmetry operators belonging to the space group of
the input experiment.expt file. It can also check for misindexing of
the diffraction pattern, possibly as a result of an incorrect beam centre.

Examples::

  dials.check_indexing_symmetry indexed.expt indexed.refl \\
    grid=1 symop_threshold=0.7


  dials.check_indexing_symmetry indexed.expt indexed.refl \\
    grid_l=3 symop_threshold=0.7
"""

phil_scope = iotbx.phil.parse(
    """
d_min = 0
  .type = float
  .help = "High resolution limit to use for analysis"
d_max = 0
  .type = float
  .help = "Low resolution limit to use for analysis"
symop_threshold = 0
  .type = float
  .help = "Threshold above which we consider a symmetry operator true"
          "at approximately 95% confidence."
grid = 0
  .type = int
  .help = "Search scope for testing misindexing on h, k, l."
grid_h = 0
  .type = int
  .help = "Search scope for testing misindexing on h."
grid_k = 0
  .type = int
  .help = "Search scope for testing misindexing on k."
grid_l = 0
  .type = int
  .help = "Search scope for testing misindexing on l."
asu = False
  .type = bool
  .help = "Perform search comparing within ASU (assumes input symm)"
normalise = False
  .type = bool
  .help = "Normalise intensities before calculating correlation coefficients."
normalise_bins = 0
  .type = int
  .help = "Number of resolution bins for normalisation"
reference = None
  .type = path
  .help = "Correctly indexed reference set for comparison"
output {
  log = dials.check_indexing_symmetry.log
    .type = str
}
""",
    process_includes=True,
)


def standard_error_of_pearson_cc(corr_coeffs, n):
    """Calculate the standard error of Pearson's correlation coefficient.

    This uses the formula SE = sqrt((1 - r^2) / (n - 2)), which is correct
    only if the data from which the CCs were calculated are normally
    distributed.

    Args:
        corr_coeffs: a flex.double array of correlation coefficient values
        n: a flex.int array of sample sizes used to calculate corr_coeffs

    Returns:
        flex.double: the standard error for values in corr_coeffs
    """

    numerator = 1 - flex.pow2(corr_coeffs)
    denominator = n - 2
    sel = denominator < 1
    denominator.set_selected(sel, 1)
    vals = numerator / denominator.as_double()
    vals.set_selected(sel, float("NaN"))
    return vals


def format_cc_with_standard_error(cc, se, decimal_places=3):
    if math.isnan(se):
        return "?"
    minimum = pow(0.1, decimal_places)
    if se >= minimum:
        cc_str = format_float_with_standard_uncertainty(cc, se, minimum=1.0e-3)
    else:
        cc_str = "{val:0.{dp}f}".format(val=cc, dp=decimal_places)
    return cc_str


def get_symop_correlation_coefficients(miller_array, use_binning=False):
    corr_coeffs = flex.double()
    n_refs = flex.int()
    space_group = miller_array.space_group()
    for smx in space_group.smx():
        reindexed_array = miller_array.change_basis(sgtbx.change_of_basis_op(smx))
        intensity, intensity_rdx = reindexed_array.common_sets(miller_array)
        if use_binning:
            intensity.use_binning_of(miller_array)
            intensity_rdx.use_binning_of(miller_array)
            cc = intensity.correlation(intensity_rdx, use_binning=use_binning)
            corr_coeffs.append(
                flex.mean_weighted(
                    flex.double(i for i in cc.data if i is not None),
                    flex.double(
                        j for i, j in zip(cc.data, cc.binner.counts()) if i is not None
                    ),
                )
            )
        else:
            corr_coeffs.append(
                intensity.correlation(
                    intensity_rdx, use_binning=use_binning
                ).coefficient()
            )
        n_refs.append(intensity.size())
    return corr_coeffs, n_refs


def normalise_intensities(miller_array, n_bins=10):
    miller_array.setup_binner(n_bins=n_bins)
    nomalisations = miller_array.amplitude_quasi_normalisations()
    miller_array = miller_array.customized_copy(
        data=miller_array.data() / nomalisations.data()
    )
    return miller_array


def test_crystal_pointgroup_symmetry(reflections, experiment, params):
    crystal = experiment.crystal

    # in case we pass in reflections from integration
    reflections = reflections.select(reflections["intensity.sum.variance"] > 0)
    reflections = reflections.select(reflections["intensity.sum.value"] > 0)
    original_miller_indices = reflections["miller_index"]

    space_group = crystal.get_space_group()
    unit_cell = crystal.get_unit_cell()

    cs = crystal_symmetry(unit_cell, space_group.type().lookup_symbol())

    ms = miller_set(cs, original_miller_indices)
    ms = ms.array(
        reflections["intensity.sum.value"]
        / flex.sqrt(reflections["intensity.sum.variance"])
    )

    if params.d_min or params.d_max:
        d_spacings = ms.d_spacings().data()
        sel = d_spacings >= params.d_min
        if params.d_max > 0:
            sel = sel & (d_spacings <= params.d_max)
        ms = ms.select(sel)

    if params.normalise:
        if params.normalise_bins:
            ms = normalise_intensities(ms, n_bins=params.normalise_bins)
        else:
            ms = normalise_intensities(ms)

    logger.info("Check symmetry operations on %d reflections:", ms.size())
    logger.info("")
    logger.info("%20s %6s %5s", "Symop", "Nref", "CC")

    true_symops = []

    ccs, n_refs = get_symop_correlation_coefficients(ms)
    ses = standard_error_of_pearson_cc(ccs, n_refs)

    for smx, cc, n_ref, se in zip(space_group.smx(), ccs, n_refs, ses):
        accept = ""
        if params.symop_threshold:
            if (cc - 2.0 * se) > params.symop_threshold:
                true_symops.append(smx)
                accept = "***"
        cc_str = format_cc_with_standard_error(cc, se)
        logger.info("%20s %6d %s %s", smx, n_ref, cc_str, accept)

    if params.symop_threshold:
        sg = sgtbx_space_group()
        for symop in true_symops:
            sg = sg.expand_smx(symop)
        for ltr in space_group.ltr():
            sg = sg.expand_ltr(ltr)
        sg_symbols = sg.match_tabulated_settings()
        logger.info("")
        logger.info(
            "Derived point group from symmetry operations: %s",
            sg_symbols.hermann_mauguin(),
        )
        logger.info("")

    return


def offset_miller_indices(miller_indices, offset):
    return flex.miller_index(
        *[mi.iround() for mi in (miller_indices.as_vec3_double() + offset).parts()]
    )


def get_indexing_offset_correlation_coefficients(
    reflections,
    crystal,
    grid,
    d_min=None,
    d_max=None,
    map_to_asu=False,
    grid_h=0,
    grid_k=0,
    grid_l=0,
    reference=None,
):
    if grid:
        if grid_h == 0:
            grid_h = grid
        if grid_k == 0:
            grid_k = grid
        if grid_l == 0:
            grid_l = grid

    return origin.get_hkl_offset_correlation_coefficients(
        reflections,
        crystal,
        map_to_asu=map_to_asu,
        grid_h=grid_h,
        grid_k=grid_k,
        grid_l=grid_l,
        reference=reference,
    )


def test_P1_crystal_indexing(reflections, experiment, params):
    if not (params.grid or params.grid_h or params.grid_k or params.grid_l):
        return

    logger.info("Checking HKL origin:")
    logger.info("")
    logger.info("dH dK dL %6s %5s", "Nref", "CC")

    if params.reference:
        reference = flex.reflection_table.from_file(params.reference)
    else:
        reference = None

    offsets, ccs, n_refs = get_indexing_offset_correlation_coefficients(
        reflections,
        experiment.crystal,
        grid=params.grid,
        d_min=params.d_min,
        d_max=params.d_max,
        map_to_asu=params.asu,
        grid_h=params.grid_h,
        grid_k=params.grid_k,
        grid_l=params.grid_l,
        reference=reference,
    )

    ses = standard_error_of_pearson_cc(ccs, n_refs)

    for (h, k, l), cc, n, se in zip(offsets, ccs, n_refs, ses):
        if cc > params.symop_threshold or (h == k == l == 0):
            cc_str = format_cc_with_standard_error(cc, se)
            logger.info("%2d %2d %2d %6d %s", h, k, l, n, cc_str)

    logger.info("")

    return


@show_mail_handle_errors()
def run(args=None):
    usage = "dials.check_indexing_symmetry [options] indexed.expt indexed.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=True)

    # Configure the logging
    log.config(logfile=params.output.log)
    logger.info(dials_version())

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    if len(reflections) == 0 or len(experiments) == 0:
        parser.print_help()
        return
    assert len(reflections) == 1
    assert len(experiments) == 1
    experiment = experiments[0]
    reflections = reflections[0]

    # remove reflections with 0, 0, 0 index
    zero = reflections["miller_index"] == (0, 0, 0)
    logger.info("Removing %d unindexed reflections", zero.count(True))
    reflections = reflections.select(~zero)

    h, k, l = reflections["miller_index"].as_vec3_double().parts()

    h = h.iround()
    k = k.iround()
    l = l.iround()

    logger.info("Range on h: %d to %d", flex.min(h), flex.max(h))
    logger.info("Range on k: %d to %d", flex.min(k), flex.max(k))
    logger.info("Range on l: %d to %d", flex.min(l), flex.max(l))

    test_P1_crystal_indexing(reflections, experiment, params)
    test_crystal_pointgroup_symmetry(reflections, experiment, params)


if __name__ == "__main__":
    run()
