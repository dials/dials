"""
Identify connected regions of missing reflections in the asymmetric unit.

This is achieved by first generating the complete set of possible miller indices,
then performing connected components analysis on a graph of nearest neighbours in
the list of missing reflections.

Examples::

  dials.missing_reflections integrated.expt integrated.refl

  dials.missing_reflections scaled.expt scaled.refl min_component_size=10
"""

from __future__ import annotations

import io
import logging
import sys
from typing import List

import libtbx.phil
from cctbx import uctbx

import dials.util.log
from dials.report.analysis import scaled_data_as_miller_array
from dials.util import missing_reflections, tabulate
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections
from dials.util.options import ArgumentParser, flatten_experiments, flatten_reflections
from dials.util.version import dials_version

logger = logging.getLogger("dials.missing_reflections")

phil_scope = libtbx.phil.parse(
    """
    min_component_size = 0
      .type = int(value_min=0)
      .help = "Only show connected regions larger than or equal to this."
    d_min = None
      .type = float(value_min=0)
    d_max = None
      .type = float(value_min=0)
    """
)


@dials.util.show_mail_handle_errors()
def run(args: List[str] = None, phil: libtbx.phil.scope = phil_scope) -> None:
    usage = "dials.missing_reflections [options] scaled.expt scaled.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging.
    dials.util.log.config(options.verbose)
    logger.info(dials_version())

    # Log the difference between the PHIL scope definition and the active PHIL scope,
    # which will include the parsed user inputs.
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    if not experiments or not reflections:
        parser.print_help()
        return
    if len(reflections) != 1 and len(experiments) != len(reflections):
        sys.exit("Number of experiments must equal the number of reflection tables")

    from dials.util.multi_dataset_handling import (
        assign_unique_identifiers,
        parse_multiple_datasets,
    )

    reflections = parse_multiple_datasets(reflections)
    experiments, reflections = assign_unique_identifiers(experiments, reflections)

    if all("inverse_scale_factor" in refl for refl in reflections):
        # Assume all arrays have been scaled
        miller_array = scaled_data_as_miller_array(
            reflections, experiments, anomalous_flag=False
        )
    else:
        # Else get the integrated intensities
        miller_arrays = filtered_arrays_from_experiments_reflections(
            experiments,
            reflections,
        )
        miller_array = miller_arrays[0]
        for ma in miller_arrays[1:]:
            miller_array = miller_array.concatenate(ma)

    if params.d_min or params.d_max:
        miller_array = miller_array.resolution_filter(
            d_min=params.d_min, d_max=params.d_max
        )

    # Print overall summary of input miller array
    s = io.StringIO()
    ma_unique = miller_array.unique_under_symmetry()
    ma_unique.show_comprehensive_summary(f=s)
    logger.info(f"\n{s.getvalue()}")

    # Get the regions of missing reflections
    complete_set, unique_ms = missing_reflections.connected_components(miller_array)
    unique_ms = [ms for ms in unique_ms if ms.size() >= params.min_component_size]

    # Print some output for user
    if len(unique_ms):
        logger.info(
            "The following connected regions of missing reflections have been identified:"
        )
        n_expected = complete_set.size()
        rows = []
        for ms in unique_ms:
            d_max, d_min = (uctbx.d_star_sq_as_d(ds2) for ds2 in ms.min_max_d_star_sq())
            rows.append(
                [
                    ms.size(),
                    f"{100 * ms.size() / n_expected:.1f}",
                    f"{d_max:.2f}-{d_min:.2f}",
                ]
            )
        logger.info(
            tabulate(
                rows, headers=["# reflections", "% missing", "Resolution range (Å)"]
            )
        )
    else:
        logger.info("No connected regions of missing reflections identified")


if __name__ == "__main__":
    run()
