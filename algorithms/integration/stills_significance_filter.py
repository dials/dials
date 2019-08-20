from __future__ import absolute_import, division, print_function

import math

from cctbx.crystal import symmetry
from dials.array_family import flex
from libtbx import table_utils
from libtbx.phil import parse

"""
Class to find a reasonable cutoff for integration based on work in LABELIT.
Bins the reflection table by resolution, then finds the first bin that goes
under a given I/sigI limit.  Cuts the data at that point.
"""

phil_scope = parse(
    """
significance_filter
  .expert_level=1
{
  enable = False
    .type=bool
    .help=If enabled, the significance filter will, for every experiment, find \
          the highest resolution where the I/sigI remains above a certain point\
          (controlled by isigi_cutoff).
  d_min = None
    .type = float
    .help = High resolution cutoff for binning. If None, use the highest \
            resolution reflection as d_min.
  n_bins = 20
    .type = int
    .help = Number of bins to use when examining resolution falloff
  isigi_cutoff = 1.0
    .type = float
    .help = I/sigI cutoff. Reflections in and past the first bin that falls \
            below this cutoff will not be retained
}
"""
)


class SignificanceFilter(object):
    def __init__(self, params):
        self.params = params.significance_filter
        self.best_d_min = None

    def __call__(self, experiments, reflections):
        results = flex.reflection_table()
        table_header = ["", "", "", "I", "IsigI", "N >", "RMSD", "Cutoff"]
        table_header2 = [
            "Bin",
            "Resolution Range",
            "Completeness",
            "",
            "",
            "cutoff",
            "(um)",
            "",
        ]

        for exp_id in range(len(experiments)):
            print("*" * 80)
            print("Significance filtering experiment", exp_id)
            table_data = []
            table_data.append(table_header)
            table_data.append(table_header2)
            experiment = experiments[exp_id]

            # Find the bins for this experiment
            crystal = experiment.crystal
            refls = reflections.select(reflections["id"] == exp_id)
            sym = symmetry(
                unit_cell=crystal.get_unit_cell(), space_group=crystal.get_space_group()
            )
            d = crystal.get_unit_cell().d(refls["miller_index"])
            mset = sym.miller_set(indices=refls["miller_index"], anomalous_flag=False)
            binner = mset.setup_binner(
                n_bins=self.params.n_bins, d_min=self.params.d_min or 0
            )
            acceptable_resolution_bins = []

            # Iterate through the bins, examining I/sigI at each bin
            for i in binner.range_used():
                d_max, d_min = binner.bin_d_range(i)
                if d_max < 0:
                    sel = d > d_min
                else:
                    sel = (d <= d_max) & (d > d_min)
                sel &= refls["intensity.sum.value"] > 0
                bin_refls = refls.select(sel)
                n_refls = len(bin_refls)
                avg_i = (
                    flex.mean(bin_refls["intensity.sum.value"]) if n_refls > 0 else 0
                )
                avg_i_sigi = (
                    flex.mean(
                        bin_refls["intensity.sum.value"]
                        / flex.sqrt(bin_refls["intensity.sum.variance"])
                    )
                    if n_refls > 0
                    else 0
                )
                acceptable_resolution_bins.append(
                    avg_i_sigi >= self.params.isigi_cutoff
                )

                bright_refls = bin_refls.select(
                    (
                        bin_refls["intensity.sum.value"]
                        / flex.sqrt(bin_refls["intensity.sum.variance"])
                    )
                    >= self.params.isigi_cutoff
                )
                n_bright = len(bright_refls)

                rmsd_obs = (
                    1000
                    * math.sqrt(
                        (
                            bright_refls["xyzcal.mm"] - bright_refls["xyzobs.mm.value"]
                        ).sum_sq()
                        / n_bright
                    )
                    if n_bright > 0
                    else 0
                )

                table_row = []
                table_row.append("%3d" % i)
                table_row.append(
                    "%-13s"
                    % binner.bin_legend(
                        i_bin=i,
                        show_bin_number=False,
                        show_bin_range=False,
                        show_d_range=True,
                        show_counts=False,
                    )
                )
                table_row.append(
                    "%13s"
                    % binner.bin_legend(
                        i_bin=i,
                        show_bin_number=False,
                        show_bin_range=False,
                        show_d_range=False,
                        show_counts=True,
                    )
                )

                table_row.append("%.2f" % (avg_i))
                table_row.append("%.2f" % (avg_i_sigi))
                table_row.append("%3d" % n_bright)
                table_row.append("%.2f" % (rmsd_obs))
                table_data.append(table_row)

            # Throw out bins that go back above the cutoff after the first non-passing bin is found
            acceptable_resolution_bins = [
                acceptable_resolution_bins[i]
                for i in range(len(acceptable_resolution_bins))
                if False not in acceptable_resolution_bins[: i + 1]
            ]

            for b, row in zip(acceptable_resolution_bins, table_data[2:]):
                if b:
                    row.append("X")
            print(
                table_utils.format(
                    table_data, has_header=2, justify="center", delim=" "
                )
            )

            # Save the results
            if any(acceptable_resolution_bins):
                best_index = acceptable_resolution_bins.count(True) - 1
                best_row = table_data[best_index + 2]
                d_min = binner.bin_d_range(binner.range_used()[best_index])[1]
                self.best_d_min = d_min
                print("best row:", " ".join(best_row))
                if self.params.enable:
                    results.extend(refls.select(d >= d_min))
            else:
                print("Data didn't pass cutoff")
                self.best_d_min = None
        if self.params.enable:
            return results
        else:
            return reflections
