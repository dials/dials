from __future__ import annotations

import json
import sys

import iotbx.phil

import dials.util
from dials.algorithms.spot_finding import per_image_analysis
from dials.util import tabulate
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

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
id = None
  .type = int(value_min=0)
"""
)


@dials.util.show_mail_handle_errors()
def run(args=None):
    usage = "dials.spot_counts_per_image [options] imported.expt strong.refl"

    parser = ArgumentParser(
        usage=usage,
        read_reflections=True,
        read_experiments=True,
        phil=phil_scope,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=False)
    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    if not reflections and not experiments:
        parser.print_help()
        return

    # FIXME may want to change this to allow many to be passed i.e.
    # from parallel runs
    if len(reflections) != 1:
        sys.exit("Only one reflection list may be passed")
    reflections = reflections[0]

    if "miller_index" in reflections:
        sys.exit("Only unindexed reflections are currently supported")

    if any(experiments.crystals()):
        sys.exit("Only unindexed experiments are currently supported")

    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)

    if params.id is not None:
        reflections = reflections.select(reflections["id"] == params.id)

    all_stats = []
    for i, expt in enumerate(experiments):
        refl = reflections.select(reflections["id"] == i)
        stats = per_image_analysis.stats_per_image(
            expt, refl, resolution_analysis=params.resolution_analysis
        )
        all_stats.append(stats)

    # transpose stats
    summary_table = {}
    for s in all_stats:
        for k, value in s._asdict().items():
            summary_table.setdefault(k, [])
            summary_table[k].extend(value)
    stats = per_image_analysis.StatsMultiImage(**summary_table)
    print(stats)

    overall_stats = per_image_analysis.stats_for_reflection_table(
        reflections, resolution_analysis=params.resolution_analysis
    )
    rows = [
        ("Overall statistics", ""),
        ("#spots", "%i" % overall_stats.n_spots_total),
        ("#spots_no_ice", "%i" % overall_stats.n_spots_no_ice),
        ("d_min", f"{overall_stats.estimated_d_min:.2f}"),
        (
            "d_min (distl method 1)",
            "%.2f (%.2f)"
            % (overall_stats.d_min_distl_method_1, overall_stats.noisiness_method_1),
        ),
        (
            "d_min (distl method 2)",
            "%.2f (%.2f)"
            % (overall_stats.d_min_distl_method_2, overall_stats.noisiness_method_2),
        ),
    ]
    print(tabulate(rows, headers="firstrow"))

    if params.json:
        if params.split_json:
            for k, v in stats._asdict().items():
                start, end = params.json.split(".")
                with open(f"{start}_{k}.{end}", "w") as fp:
                    json.dump(v, fp)
        if params.joint_json:
            with open(params.json, "w") as fp:
                json.dump(stats._asdict(), fp)
    if params.plot:
        import matplotlib

        matplotlib.use("Agg")
        per_image_analysis.plot_stats(stats, filename=params.plot)


if __name__ == "__main__":
    run()
