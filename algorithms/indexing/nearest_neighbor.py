from __future__ import absolute_import, division, print_function

import math


class NeighborAnalysis(object):
    def __init__(
        self,
        reflections,
        step_size=45,
        tolerance=1.5,
        max_height_fraction=0.25,
        percentile=None,
        histogram_binning="linear",
        nn_per_bin=5,
    ):
        self.tolerance = tolerance  # Margin of error for max unit cell estimate
        from scitbx.array_family import flex

        NEAR = 10
        self.NNBIN = nn_per_bin  # target number of neighbors per histogram bin
        self.histogram_binning = histogram_binning

        direct = flex.double()

        if "entering" in reflections:
            entering_flags = reflections["entering"]
        else:
            entering_flags = flex.bool(reflections.size(), True)
        rs_vectors = reflections["rlp"]
        phi_deg = reflections["xyzobs.mm.value"].parts()[2] * (180 / math.pi)

        d_spacings = flex.double()
        # nearest neighbor analysis
        from annlib_ext import AnnAdaptor

        for imageset_id in range(flex.max(reflections["imageset_id"]) + 1):
            sel_imageset = reflections["imageset_id"] == imageset_id
            if sel_imageset.count(True) == 0:
                continue
            phi_min = flex.min(phi_deg.select(sel_imageset))
            phi_max = flex.max(phi_deg.select(sel_imageset))
            d_phi = phi_max - phi_min
            n_steps = max(int(math.ceil(d_phi / step_size)), 1)

            for n in range(n_steps):
                sel_step = (
                    sel_imageset
                    & (phi_deg >= (phi_min + n * step_size))
                    & (phi_deg < (phi_min + (n + 1) * step_size))
                )

                for entering in (True, False):
                    sel_entering = sel_step & (entering_flags == entering)
                    if sel_entering.count(True) == 0:
                        continue

                    query = flex.double()
                    query.extend(rs_vectors.select(sel_entering).as_double())

                    if query.size() == 0:
                        continue

                    IS_adapt = AnnAdaptor(data=query, dim=3, k=1)
                    IS_adapt.query(query)

                    direct.extend(1 / flex.sqrt(IS_adapt.distances))
                    d_spacings.extend(1 / rs_vectors.norms())

        assert (
            len(direct) > NEAR
        ), "Too few spots (%d) for nearest neighbour analysis." % len(direct)

        perm = flex.sort_permutation(direct)
        direct = direct.select(perm)
        d_spacings = d_spacings.select(perm)

        # eliminate nonsensical direct space distances
        sel = direct > 1
        direct = direct.select(sel)
        d_spacings = d_spacings.select(sel)

        if percentile is None:
            # reject top 1% of longest distances to hopefully get rid of any outliers
            n = int(math.floor(0.99 * len(direct)))
            direct = direct[:n]
            d_spacings = d_spacings[:n]

        # determine the most probable nearest neighbor distance (direct space)
        if self.histogram_binning == "log":
            hst = flex.histogram(
                flex.log10(direct), n_slots=int(len(direct) / self.NNBIN)
            )
        else:
            hst = flex.histogram(direct, n_slots=int(len(direct) / self.NNBIN))
        if self.histogram_binning == "log":
            self.slot_start = flex.double(
                [10 ** (s - 0.5 * hst.slot_width()) for s in hst.slot_centers()]
            )
            self.slot_end = flex.double(
                [10 ** (s + 0.5 * hst.slot_width()) for s in hst.slot_centers()]
            )
            self.slot_width = self.slot_end - self.slot_start
        else:
            self.slot_start = hst.slot_centers() - 0.5 * hst.slot_width()
            self.slot_end = hst.slot_centers() + 0.5 * hst.slot_width()
            self.slot_width = hst.slot_width()
        self.relative_frequency = hst.slots().as_double() / self.slot_width
        highest_bin_height = flex.max(self.relative_frequency)

        if percentile is not None:
            # determine the nth-percentile direct-space distance
            perm = flex.sort_permutation(direct, reverse=True)
            self.max_cell = (
                self.tolerance * direct[perm[int((1 - percentile) * len(direct))]]
            )

        else:
            # choose a max cell based on bins above a given fraction of the highest bin height
            # given multiple
            isel = (
                self.relative_frequency.as_double()
                > (max_height_fraction * highest_bin_height)
            ).iselection()
            self.max_cell = (
                self.tolerance * self.slot_end[int(flex.max(isel.as_double()))]
            )

        self.reciprocal_lattice_vectors = rs_vectors
        self.d_spacings = d_spacings
        self.direct = direct
        self.histogram = hst

    def plot_histogram(self, filename="nn_hist.png", figsize=(12, 8)):
        import matplotlib.pyplot as plt

        plt.figure(figsize=figsize)
        plt.bar(
            self.slot_start,
            self.relative_frequency,
            align="center",
            width=self.slot_width,
            color="black",
            edgecolor=None,
        )
        ymin, ymax = plt.ylim()
        if self.histogram_binning == "log":
            ax = plt.gca()
            ax.set_xscale("log")
        plt.vlines(
            self.max_cell / self.tolerance,
            ymin,
            ymax,
            linestyles="--",
            colors="g",
            label="estimated max cell",
        )
        plt.vlines(
            self.max_cell,
            ymin,
            ymax,
            colors="g",
            label="estimated max cell (including tolerance)",
        )
        plt.xlabel("Direct space distance (A)")
        plt.ylabel("Frequency")
        plt.legend(loc="upper left")
        plt.savefig(filename)
        plt.clf()
