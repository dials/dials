# coding: utf-8
"""
Command line script to allow merging and truncating of a dials dataset.
"""
from __future__ import absolute_import, division, print_function

import logging
import sys

from six.moves import cStringIO as StringIO

from libtbx import phil

from dials.algorithms.merging.merge import (
    make_MAD_merged_mtz_file,
    make_merged_mtz_file,
    merge_and_truncate,
)
from dials.util import Sorry, log, show_mail_handle_errors
from dials.util.export_mtz import match_wavelengths
from dials.util.options import OptionParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

help_message = """
Merge scaled dials data.

Examples::

  dials.merge scaled.expt scaled.refl

  dials.merge scaled.expt scaled.refl truncate=False
"""

logger = logging.getLogger("dials")
phil_scope = phil.parse(
    """
assess_space_group = True
    .type = bool
    .help = "Option to assess space group by testing presence of axial reflections"
anomalous = True
    .type = bool
    .help = "Output anomalous as well as mean intensities."
truncate = True
    .type = bool
    .help = "Option to perform truncation on merged data."
d_min = None
    .type = float
    .help = "High resolution limit to apply to the data."
d_max = None
    .type = float
    .help = "Low resolution limit to apply to the data."
combine_partials = True
    .type = bool
    .help = "Combine partials that have the same partial id into one
        reflection, with an updated partiality given by the sum of the
        individual partialities."
partiality_threshold=0.4
    .type = float
    .help = "All reflections with partiality values above the partiality
        threshold will be retained. This is done after any combination of
        partials if applicable."
n_residues = 200
    .type = int
    .help = "Number of residues to use in Wilson scaling"
merging {
    use_internal_variance = False
        .type = bool
    n_bins = 20
        .type = int(value_min=5)
    anomalous = False
        .type = bool
        .help = "Option to control whether reported merging stats are anomalous."
}
reporting {
    wilson_stats = True
        .type = bool
        .help = "Option to turn off reporting of Wilson statistics"
    merging_stats = True
        .type = bool
        .help = "Option to turn off reporting of merging statistics."
}
output {
    log = dials.merge.log
        .type = str
    mtz = merged.mtz
        .type = str
        .help = "Filename to use for mtz output."
    crystal_names = XTAL
        .type = strings
        .help = "Crystal name to be used in MTZ file output (multiple names
            allowed for MAD datasets)"
    project_name = AUTOMATIC
        .type = str
        .help = "Project name to be used in MTZ file output"
    dataset_names = NATIVE
        .type = strings
        .help = "Dataset name to be used in MTZ file output (multiple names
            allowed for MAD datasets)"
}
include scope cctbx.french_wilson.master_phil
""",
    process_includes=True,
)


def merge_data_to_mtz(params, experiments, reflections):
    """Merge data (at each wavelength) and write to an mtz file object."""
    wavelengths = match_wavelengths(experiments)
    if len(wavelengths) > 1:
        logger.info(
            "Multiple wavelengths found: \n%s",
            "\n".join(
                "  Wavlength: %.5f, experiment numbers: %s "
                % (k, ",".join(map(str, v)))
                for k, v in wavelengths.items()
            ),
        )
        return make_MAD_merged_mtz_file(params, experiments, reflections, wavelengths)
    merged_data = merge_and_truncate(params, experiments, reflections)
    return make_merged_mtz_file(*((params, list(wavelengths)[0]) + merged_data))


@show_mail_handle_errors()
def run(args=None):
    """Run the merging from the command-line."""
    usage = """Usage: dials.merge scaled.refl scaled.expt [options]"""

    parser = OptionParser(
        usage=usage,
        read_experiments=True,
        read_reflections=True,
        phil=phil_scope,
        check_format=False,
        epilog=help_message,
    )
    params, options = parser.parse_args(args=args, show_diff_phil=False)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    log.config(verbosity=options.verbose, logfile=params.output.log)
    logger.info(dials_version())

    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    ### Assert that all data have been scaled with dials - should only be
    # able to input one reflection table and experimentlist that are
    # matching and scaled together.

    if len(reflections) != 1:
        raise Sorry(
            """Only data scaled together as a single reflection dataset
can be processed with dials.merge"""
        )

    for k in [
        "intensity.scale.value",
        "intensity.scale.variance",
        "inverse_scale_factor",
    ]:
        if k not in reflections[0]:
            raise Sorry(
                """%s not found in the reflection table.
Only scaled data can be processed with dials.merge"""
                % k
            )

    try:
        mtz_file = merge_data_to_mtz(params, experiments, reflections)
    except ValueError as e:
        raise Sorry(e)

    logger.info("\nWriting reflections to %s", (params.output.mtz))
    out = StringIO()
    mtz_file.show_summary(out=out)
    logger.info(out.getvalue())
    mtz_file.write(params.output.mtz)


if __name__ == "__main__":
    run()
