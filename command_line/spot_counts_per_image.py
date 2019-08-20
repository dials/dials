from __future__ import absolute_import, division, print_function

# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import sys

import iotbx.phil
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections, flatten_experiments
from dials.algorithms.spot_finding import per_image_analysis

help_message = """

Reports the number of strong spots and computes an estimate of the resolution
limit for each image, given the results of dials.find_spots. Optionally
generates a plot of the per-image statistics (plot=image.png).

Examples::

  dials.spot_counts_per_image imported.expt strong.refl

  dials.spot_counts_per_image imported.expt strong.refl plot=per_image.png

"""

phil_scope = iotbx.phil.parse(
    """\
resolution_analysis = True
  .type = bool
plot = None
  .type = path
json = None
  .type = path
split_json = False
  .type = bool
joint_json = True
  .type = bool
individual_plots = False
  .type = bool
id = None
  .type = int(value_min=0)
"""
)


def run(args):
    usage = "dials.spot_counts_per_image [options] imported.expt strong.refl"

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

    # FIXME may want to change this to allow many to be passed i.e.
    # from parallel runs
    if len(reflections) != 1:
        sys.exit("Only one reflection list may be passed")
    reflections = reflections[0]
    expts = set(reflections["id"])
    if max(expts) >= len(experiments.imagesets()):
        sys.exit("Unknown experiments in reflection list")

    if params.id is not None:
        reflections = reflections.select(reflections["id"] == params.id)

    all_stats = []
    for j, imageset in enumerate(experiments.imagesets()):
        refl = reflections.select(reflections["id"] == j)
        stats = per_image_analysis.stats_imageset(
            imageset,
            refl,
            resolution_analysis=params.resolution_analysis,
            plot=params.individual_plots,
        )
        all_stats.append(stats)

    # transpose stats
    class empty(object):
        pass

    e = empty()
    for s in all_stats:
        for k in dir(s):
            if k.startswith("_") or k in ["merge", "next"]:
                continue
            if not hasattr(e, k):
                setattr(e, k, [])
            getattr(e, k).extend(getattr(s, k))

    per_image_analysis.print_table(e)
    from libtbx import table_utils

    # FIXME this is now probably nonsense...
    overall_stats = per_image_analysis.stats_single_image(
        imageset, reflections, resolution_analysis=params.resolution_analysis
    )
    rows = [
        ("Overall statistics", ""),
        ("#spots", "%i" % overall_stats.n_spots_total),
        ("#spots_no_ice", "%i" % overall_stats.n_spots_no_ice),
        ("d_min", "%.2f" % overall_stats.estimated_d_min),
        (
            "d_min (distl method 1)",
            "%.2f (%.2f)"
            % (overall_stats.d_min_distl_method_1, overall_stats.noisiness_method_1),
        ),
        (
            "d_min (distl method 2)",
            "%.2f (%.2f)"
            % (overall_stats.d_min_distl_method_1, overall_stats.noisiness_method_1),
        ),
    ]
    print(table_utils.format(rows, has_header=True, prefix="| ", postfix=" |"))

    if params.json:
        import json

        if params.split_json:
            for k in stats.__dict__:
                start, end = params.json.split(".")
                with open("%s_%s.%s" % (start, k, end), "wb") as fp:
                    json.dump(stats.__dict__[k], fp)
        if params.joint_json:
            with open(params.json, "wb") as fp:
                json.dump(stats.__dict__, fp)
    if params.plot:
        import matplotlib

        matplotlib.use("Agg")
        per_image_analysis.plot_stats(stats, filename=params.plot)


if __name__ == "__main__":
    run(sys.argv[1:])
