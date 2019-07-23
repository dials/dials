#!/usr/bin/env python
# coding: utf-8
"""
Command line script to allow merging and truncating of a dials dataset.
"""
from __future__ import absolute_import, division, print_function
import logging
import sys
from cStringIO import StringIO
from dials.util import log, show_mail_on_error, Sorry
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.util.version import dials_version
from dials.util.export_mtz import match_wavelengths
from dials.algorithms.merging.merge import (
    make_MAD_merged_mtz_file,
    make_merged_mtz_file,
    merge_and_truncate,
)
from libtbx import phil


help_message = """Program to merge scaled dials data."""

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
    .help = "Resolution limit to apply to the data."
n_residues = 200
    .type = int
    .help = "Number of residues to use in Wilson scaling"
merging {
    use_internal_variance = False
        .type = bool
    n_bins = 20
        .type = int(value_min=5)
}
reporting {
    wilson_stats = True
        .type = bool
        .help = "Option to turn off reporting of wilson statistics"
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
}
include scope cctbx.french_wilson.master_phil
""",
    process_includes=True,
)


def merge_data_to_mtz(params, experiments, reflections):
    """Merge data (at each wavelength) and write to an mtz file object."""
    if len(experiments) > 1:
        wavelengths = match_wavelengths(experiments)
        if len(wavelengths.keys()) > 1:
            logger.info(
                "Multiple wavelengths found: \n%s",
                "\n".join(
                    "  Wavlength: %.5f, experiment numbers: %s "
                    % (k, ",".join(map(str, v)))
                    for k, v in wavelengths.iteritems()
                ),
            )
            return make_MAD_merged_mtz_file(
                params, experiments, reflections, wavelengths
            )
    merged_data = merge_and_truncate(params, experiments, reflections)
    return make_merged_mtz_file(*merged_data)


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
    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    log.config(verbosity=1, info=params.output.log)
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
            """Only data scaled together in a single reflection data
can be processed with dials.merge"""
        )

    for k in [
        "intensity.scale.value",
        "intensity.scale.variance",
        "inverse_scale_factor",
    ]:
        if k not in reflections[0]:
            raise Sorry(
                """%s not found in the reflection table. Only scaled data can be processed
with dials.merge"""
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
    with show_mail_on_error():
        run()
