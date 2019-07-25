"""Creation of 'corrgram' correlation matrix plots"""

from __future__ import absolute_import, division, print_function

import json
import math
import os
import logging

import libtbx.phil

logger = logging.getLogger(__name__)

# PHIL options for corrgrams
phil_str = """
    correlation_plot
      .expert_level = 1
    {
      filename = None
        .type = str
        .help = "The base filename for output of plots of parameter"
                "correlations. A file extension may be added to control"
                "the type of output file, if it is one of matplotlib's"
                "supported types. A JSON file with the same base filename"
                "will also be created, containing the correlation matrix and"
                "column labels for later inspection, replotting etc."

      col_select = None
        .type = strings
        .help = "Specific columns to include in the plots of parameter"
                "correlations, either specifed by parameter name or 0-based"
                "column index. Defaults to all columns."
                "This option is useful when there is a large number of"
                "parameters"

      steps = None
        .type = ints(value_min=0)
        .help = "Steps for which to make correlation plots. By default only"
                "the final step is plotted. Uses zero-based numbering, so"
                "the first step is numbered 0."
    }
"""

phil_scope = libtbx.phil.parse(phil_str)


def corrgram(corrmat, labels):
    """Create a correlation matrix plot or 'corrgram' for the provided
    correlation matrix and row/column labels. Inspired by R's corrplot and
    https://github.com/louridas/corrplot/blob/master/corrplot.py"""

    try:  # is corrmat a scitbx matrix?
        corrmat = corrmat.as_flex_double_matrix()
    except AttributeError:  # assume it is already a flex double matrix
        pass
    assert corrmat.is_square_matrix()

    nr = corrmat.all()[0]
    assert nr == len(labels)

    try:
        import matplotlib

        matplotlib.use("Agg", warn=False)
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
    except ImportError as e:
        msg = "matplotlib modules not available " + str(e)
        logger.info(msg)
        return None

    plt.figure(1)
    ax = plt.subplot(1, 1, 1, aspect="equal")
    clrmap = cm.get_cmap("bwr")

    for x in range(nr):
        for y in range(nr):
            d = corrmat[x, y]
            d_abs = abs(d)
            circ = plt.Circle((x, y), radius=0.9 * math.sqrt(d_abs) / 2)
            circ.set_edgecolor("white")
            # put data into range [0,1] and invert so that 1 == blue and 0 == red
            facecolor = 1 - (0.5 * d + 0.5)
            circ.set_facecolor(clrmap(facecolor))
            ax.add_artist(circ)
    ax.set_xlim(-0.5, nr - 0.5)
    ax.set_ylim(-0.5, nr - 0.5)

    ax.xaxis.tick_top()
    xtickslocs = list(range(len(labels)))
    ax.set_xticks(xtickslocs)
    ax.set_xticklabels(labels, rotation=30, fontsize="small", ha="left")

    ax.invert_yaxis()
    ytickslocs = list(range(len(labels)))
    ax.set_yticks(ytickslocs)
    ax.set_yticklabels(labels, fontsize="small")

    xtickslocs = [e + 0.5 for e in xtickslocs]
    ax.set_xticks(xtickslocs, minor=True)
    ytickslocs = [e + 0.5 for e in ytickslocs]
    ax.set_yticks(ytickslocs, minor=True)
    plt.grid(color="0.8", which="minor", linestyle="-")

    # suppress major tick marks
    ax.tick_params(which="major", width=0)

    # need this otherwise text gets clipped
    plt.tight_layout()

    # FIXME should this also have a colorbar as legend?
    return plt


def create_correlation_plots(refiner, params):
    root, ext = os.path.splitext(params.correlation_plot.filename)
    if not ext:
        ext = ".pdf"

    history = refiner.history
    steps = params.correlation_plot.steps
    if steps is None:
        steps = [history.get_nrows() - 1]

    # extract individual column names or indices
    col_select = params.correlation_plot.col_select

    num_plots = 0
    for step in steps:
        fname_base = root
        if len(steps) > 1:
            fname_base += "_step%02d" % step

        corrmats, labels = refiner.get_parameter_correlation_matrix(step, col_select)
        if [corrmats, labels].count(None) == 0:

            for resid_name, corrmat in corrmats.items():
                plot_fname = fname_base + "_" + resid_name + ext
                plt = corrgram(corrmat, labels)
                if plt is not None:
                    logger.info(
                        "Saving parameter correlation plot to {}".format(plot_fname)
                    )
                    plt.savefig(plot_fname)
                    plt.close()
                    num_plots += 1
            mat_fname = fname_base + ".json"
            with open(mat_fname, "w") as handle:
                for k, corrmat in corrmats.items():
                    corrmats[k] = corrmat.as_scitbx_matrix().as_list_of_lists()
                logger.info(
                    "Saving parameter correlation matrices to {}".format(mat_fname)
                )
                json.dump({"corrmats": corrmats, "labels": labels}, handle)

    if num_plots == 0:
        msg = (
            "Sorry, no parameter correlation plots were produced. Please set "
            "track_parameter_correlation=True to ensure correlations are "
            "tracked, and make sure correlation_plot.col_select is valid."
        )
        logger.info(msg)
