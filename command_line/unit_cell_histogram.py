from __future__ import absolute_import, division, print_function
import libtbx.phil
from dials.util import log

import logging

logger = logging.getLogger("dials.unit_cell_histogram")

help_message = """

"""

phil_scope = libtbx.phil.parse(
    """
steps_per_angstrom = 20
  .type = int
iqr_ratio = 1.5
  .type = float
""",
    process_includes=True,
)

from scitbx.array_family import flex
from dials.algorithms.clustering.observers import uc_params_from_experiments


def panel_distances_from_experiments(experiments):
    distances = [flex.double() for i in range(len(experiments[0].detector))]
    for expt in experiments:
        for i, panel in enumerate(expt.detector):
            distances[i].append(panel.get_distance())
    return distances


def outlier_selection(uc_params, iqr_ratio=1.5):
    outliers = flex.bool(uc_params[0].size(), False)
    for p in uc_params:
        from scitbx.math import five_number_summary

        min_x, q1_x, med_x, q3_x, max_x = five_number_summary(p)
        logger.info(
            "Five number summary: min %.2f, q1 %.2f, med %.2f, q3 %.2f, max %.2f"
            % (min_x, q1_x, med_x, q3_x, max_x)
        )
        iqr_x = q3_x - q1_x
        if iqr_x < 1e-6:
            continue
        cut_x = iqr_ratio * iqr_x
        outliers.set_selected(p > q3_x + cut_x, True)
        outliers.set_selected(p < q1_x - cut_x, True)
    logger.info("Identified %i unit cell outliers" % outliers.count(True))
    return outliers


def run(args):

    from dials.util.options import OptionParser
    from dials.util.options import flatten_experiments
    import libtbx.load_env

    usage = "%s [options] models.expt" % (libtbx.env.dispatcher_name)

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    from dials.util.version import dials_version

    logger.info(dials_version())

    params, options = parser.parse_args(show_diff_phil=False)
    experiments = flatten_experiments(params.input.experiments)

    if len(experiments) == 0:
        parser.print_help()
        exit(0)

    log.config()

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    uc_params = uc_params_from_experiments(experiments)
    panel_distances = panel_distances_from_experiments(experiments)
    outliers = outlier_selection(uc_params, iqr_ratio=params.iqr_ratio)
    plot_uc_histograms(uc_params, outliers, params.steps_per_angstrom)
    plot_uc_vs_detector_distance(
        uc_params, panel_distances, outliers, params.steps_per_angstrom
    )
    plot_number_of_crystals(experiments)


def plot_uc_histograms(
    uc_params, outliers, steps_per_angstrom=20, plot_name="uc_histograms.png"
):
    from matplotlib import pyplot as plt

    plt.style.use("ggplot")
    f, ax = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))
    a, b, c = uc_params[:3]

    def uc_param_hist2d(p1, p2, ax):
        nbins = 100
        import numpy as np

        H, xedges, yedges = np.histogram2d(p1, p2, bins=nbins)
        H = np.rot90(H)
        H = np.flipud(H)
        Hmasked = np.ma.masked_where(H == 0, H)
        ax.pcolormesh(xedges, yedges, Hmasked)

    uc_param_hist2d(a, b, ax[0][0])
    uc_param_hist2d(b, c, ax[0][1])
    uc_param_hist2d(c, a, ax[0][2])

    for i in range(3):
        mmm = flex.min_max_mean_double(uc_params[i])
        import math

        steps_per_A = steps_per_angstrom
        Amin = math.floor(mmm.min * steps_per_A) / steps_per_A
        Amax = math.floor(mmm.max * steps_per_A) / steps_per_A
        n_slots = max(1, int((Amax - Amin) * steps_per_A))
        if Amin == Amax:
            eps = 0.05
            Amin -= eps
            Amax += eps
        hist = flex.histogram(uc_params[i], Amin, Amax, n_slots=n_slots)
        hist_inliers = flex.histogram(
            uc_params[i].select(~outliers), Amin, Amax, n_slots=n_slots
        )
        ax[1][i].bar(
            hist.slot_centers(),
            hist.slots(),
            align="center",
            width=hist.slot_width(),
            zorder=10,
            color="black",
            edgecolor=None,
            linewidth=0,
        )
        ax[1][i].bar(
            hist_inliers.slot_centers(),
            hist_inliers.slots(),
            align="center",
            width=hist_inliers.slot_width(),
            zorder=10,
            color="red",
            edgecolor=None,
            linewidth=0,
        )
        ax[0][i].set_xlim(ax[1][i].get_xlim())

    ax[0][0].set_ylabel(r"b ($\AA$)")
    ax[0][1].set_ylabel(r"c ($\AA$)")
    ax[0][2].set_ylabel(r"a ($\AA$)")
    ax[1][0].set_xlabel(r"a ($\AA$)")
    ax[1][1].set_xlabel(r"b ($\AA$)")
    ax[1][2].set_xlabel(r"c ($\AA$)")

    f.savefig(plot_name)
    plt.tight_layout()
    plt.close(f)


def plot_uc_vs_detector_distance(
    uc_params,
    panel_distances,
    outliers,
    steps_per_angstrom=20,
    filename="uc_vs_distance.png",
):
    from matplotlib import pyplot as plt

    plt.style.use("ggplot")
    f = plt.figure(figsize=(12, 8))
    ax1 = plt.subplot2grid((2, 3), (0, 0))
    ax2 = plt.subplot2grid((2, 3), (0, 1), sharey=ax1)
    ax3 = plt.subplot2grid((2, 3), (0, 2), sharey=ax1)
    ax4 = plt.subplot2grid((2, 3), (1, 0), colspan=3)
    a, b, c = uc_params[:3]

    def hist2d(p1, p2, ax):
        nbins = 100
        import numpy as np

        H, xedges, yedges = np.histogram2d(p1, p2, bins=nbins)
        H = np.rot90(H)
        H = np.flipud(H)
        Hmasked = np.ma.masked_where(H == 0, H)
        ax.pcolormesh(xedges, yedges, Hmasked)

    hist2d(a, panel_distances[0], ax1)
    hist2d(b, panel_distances[0], ax2)
    hist2d(c, panel_distances[0], ax3)

    mmm = flex.min_max_mean_double(panel_distances[0])
    import math

    steps_per_mm = 20
    Amin = math.floor(mmm.min * steps_per_mm) / steps_per_mm
    Amax = math.floor(mmm.max * steps_per_mm) / steps_per_mm
    n_slots = int((Amax - Amin) * steps_per_mm)
    hist = flex.histogram(panel_distances[0], Amin, Amax, n_slots=n_slots)
    ax4.bar(
        hist.slot_centers(),
        hist.slots(),
        align="center",
        width=hist.slot_width(),
        zorder=10,
        color="red",
        edgecolor=None,
        linewidth=0,
    )

    ax1.set_ylabel("Detector distance (mm)")
    ax1.set_xlabel(r"a ($\AA$)")
    ax2.set_xlabel(r"b ($\AA$)")
    ax3.set_xlabel(r"c ($\AA$)")
    ax4.set_xlabel("Detector distance (mm)")

    f.savefig(filename)
    plt.tight_layout()
    f.clf()


def plot_number_of_crystals(experiments):
    image_to_expts = {}
    for expt in experiments:
        img = expt.imageset.get_image_identifier(0)
        image_to_expts.setdefault(img, [])
        image_to_expts[img].append(expt)

    n_crystals_per_image = flex.int(len(expts) for expts in image_to_expts.values())
    nmax = flex.max(n_crystals_per_image)
    hist = flex.histogram(n_crystals_per_image.as_double(), 0, nmax, n_slots=nmax)
    # hist.show()
    from matplotlib import pyplot as plt

    plt.style.use("ggplot")
    plt.bar(
        hist.slot_centers(),
        hist.slots(),
        align="center",
        width=hist.slot_width(),
        zorder=10,
        color="black",
        edgecolor=None,
    )
    plt.savefig("n_crystals_hist.png")
    plt.clf()


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])
