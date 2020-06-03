"""
Identify connected regions of missing reflections in the asymmetric unit.

This is achieved by first generating the complete set of possible miller indices,
then performing connected components analysis on a graph of nearest neighbours in
the list of missing reflections.

Examples::

  dials.missing_reflections integrated.expt integrated.refl

  dials.missing_reflections scaled.expt scaled.refl min_component_size=10
"""

import logging
import sys
from typing import List

import libtbx.phil
from cctbx import uctbx

import dials.util.log
from dials.report.analysis import scaled_data_as_miller_array
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections
from dials.util.options import OptionParser, flatten_experiments, flatten_reflections
from dials.util import missing_reflections
from dials.util.version import dials_version


logger = logging.getLogger("dials.missing_reflections")

phil_scope = libtbx.phil.parse(
    """
    min_component_size = 0
      .type = int(value_min=0)
      .help = "Only show connected regions larger than or equal to this."
    """
)


def run(args=None, phil=phil_scope):  # type: (List[str], libtbx.phil.scope) -> None
    usage = "dials.missing_reflections [options] scaled.expt scaled.refl"

    parser = OptionParser(
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
            experiments, reflections,
        )
        miller_array = miller_arrays[0]
        for ma in miller_arrays[1:]:
            miller_array = miller_array.concatenate(ma)

    # Get the regions of missing reflections
    complete_set, unique_ms = missing_reflections.connected_components(miller_array)

    # Print some output for user
    if len(unique_ms):
        n_expected = complete_set.size()
        unique_ms = [ms for ms in unique_ms if ms.size() >= params.min_component_size]
        for ms in unique_ms:
            d_max, d_min = (uctbx.d_star_sq_as_d(ds2) for ds2 in ms.min_max_d_star_sq())
            logger.info(
                f"{ms.size()} reflections ({100 * ms.size() / n_expected:.1f}%): {d_max:.2f}-{d_min:.2f} Å"
            )
    else:
        logger.info("No connected regions of missing reflections identified")


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()
