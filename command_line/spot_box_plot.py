from __future__ import absolute_import, division, print_function

# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# LIBTBX_SET_DISPATCHER_NAME dev.dials.spot_box_plot

import sys
import numpy

from dials.util.options import OptionParser
from dials.util.options import (
    flatten_reflections,
    flatten_experiments,
    flatten_experiments,
)
import libtbx.load_env

import iotbx.phil

help_message = (
    """

  %s imported_experiments.expt strong.refl

"""
    % libtbx.env.dispatcher_name
)

phil_scope = iotbx.phil.parse(
    """
bin = 1
    .type = int
    .help = "Number of images to bin for each interval"
"""
)


def spot_count_iqr_mean(reflections, experiment, bin_width=1):
    """Compute the spot count IQR and mean across a scan, optionally binning"""

    scan = experiment.scan
    first, last = scan.get_image_range()

    x, y, z = reflections["xyzobs.px.value"].parts()
    z2 = z / bin_width

    bins = [r / bin_width for r in range(first - 1, last, bin_width)]

    # this could just use the z2 values but allowing for e.g. total
    # signal estimation or similar

    spots_per_bin = [reflections.select((z2 > b) & (z2 < (b + 1))).size() for b in bins]

    spots_array = numpy.array(spots_per_bin)

    q25, q75 = numpy.percentile(spots_array, [25, 75])

    return numpy.mean(spots_array), q25, q75


def run(args):
    usage = (
        "%s [options] imported_experiments.expt strong.refl"
        % libtbx.env.dispatcher_name
    )

    parser = OptionParser(
        usage=usage,
        read_reflections=True,
        read_experiments=True,
        phil=phil_scope,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(show_diff_phil=False)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    if not any([reflections, experiments]):
        parser.print_help()
        return

    if len(reflections) != 1:
        sys.exit("Only one reflection list may be passed")

    if len(experiments) != 1:
        sys.exit("Only one experiment list may be passed")

    mean, q25, q75 = spot_count_iqr_mean(reflections[0], experiments[0], params.bin)
    print("%.2f" % (mean / (q75 - q25)))


if __name__ == "__main__":
    run(sys.argv[1:])
