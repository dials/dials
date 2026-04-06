from __future__ import annotations

import concurrent.futures
import math

import numpy as np
from scipy.spatial import cKDTree

# Thread count for parallelising the ~160 NN sub-groups. Each group runs
# cKDTree with workers=1; the GIL is released inside scipy's C layer, so
# Python threads genuinely run in parallel. We follow the same CPU_COUNT
# convention used elsewhere in DIALS (e.g. bravais_settings.py) rather than
# hard-coding a thread count.
from dials.util.system import CPU_COUNT

_NN_MAX_WORKERS = CPU_COUNT


def _ckdtree_nn_distances(pts_np):
    """Return 1-NN distances for pts_np (N×3 float64 array).

    scipy.cKDTree includes self-matches, so we query k=2 and take column 1
    (column 0 is always 0.0, the self-match). This is bit-identical to the
    AnnAdaptor result with ANN_ALLOW_SELF_MATCH=ANNfalse and k=1.
    """
    tree = cKDTree(pts_np, leafsize=20)
    dists, _ = tree.query(pts_np, k=2, workers=1)
    return dists[:, 1]  # actual distances in Å⁻¹ (not squared, unlike ANN)


class NeighborAnalysis:
    def __init__(
        self,
        reflections,
        step_size=45,
        tolerance=1.5,
        max_height_fraction=0.25,
        percentile=None,
        histogram_binning="linear",
        nn_per_bin=5,
        convert_reflections_z_to_deg=True,
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

        z = reflections["xyzobs.mm.value"].parts()[2]
        if convert_reflections_z_to_deg:
            z = z * (180 / math.pi)

        d_spacings = flex.double()

        # Collect all sub-groups first, then dispatch to a thread pool.
        # Splitting by (imageset_id, 45° rotation window, entering/exiting)
        # preserves the original grouping semantics exactly.
        group_pts = []  # list of N×3 float64 arrays, one per sub-group
        group_norms = []  # list of flex.double norms, one per sub-group

        for imageset_id in range(flex.max(reflections["imageset_id"]) + 1):
            sel_imageset = reflections["imageset_id"] == imageset_id
            if sel_imageset.count(True) == 0:
                continue
            z_min = flex.min(z.select(sel_imageset))
            z_max = flex.max(z.select(sel_imageset))
            d_z = z_max - z_min
            n_steps = max(int(math.ceil(d_z / step_size)), 1)

            for n in range(n_steps):
                sel_step = (
                    sel_imageset
                    & (z >= (z_min + n * step_size))
                    & (z < (z_min + (n + 1) * step_size))
                )

                for entering in (True, False):
                    sel_entering = sel_step & (entering_flags == entering)
                    if sel_entering.count(True) == 0:
                        continue

                    vecs = rs_vectors.select(sel_entering)
                    pts_np = np.array(vecs.as_double()).reshape(-1, 3)

                    if pts_np.shape[0] == 0:
                        continue

                    group_pts.append(pts_np)
                    group_norms.append(vecs.norms())

        # Parallel NN queries — cKDTree releases the GIL, so threads are effective
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=_NN_MAX_WORKERS
        ) as executor:
            all_dists = list(executor.map(_ckdtree_nn_distances, group_pts))

        for dists, norms in zip(all_dists, group_norms):
            # ANN returned squared distances; cKDTree returns actual distances.
            # Original code: 1 / flex.sqrt(IS_adapt.distances)
            # Here:          1 / dists   (equivalent, bit-identical)
            direct.extend(flex.double(dists.tolist()))
            d_spacings.extend(1 / norms)

        assert len(direct) > NEAR, (
            f"Too few spots ({len(direct)}) for nearest neighbour analysis."
        )

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
            # Vectorised form: avoids a Python loop over ~400k slot centers.
            # Equivalent to the list-comprehension [10**(s ± half_width)] above.
            centers = np.array(hst.slot_centers())
            half_width = 0.5 * hst.slot_width()
            self.slot_start = flex.double(10 ** (centers - half_width))
            self.slot_end = flex.double(10 ** (centers + half_width))
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
