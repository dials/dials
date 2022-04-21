from __future__ import annotations

import collections
import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

import numpy as np
import scipy.spatial.distance as ssd
from scipy.cluster import hierarchy

from cctbx import crystal, uctbx
from cctbx.sgtbx.lattice_symmetry import metric_subgroups
from cctbx.uctbx.determine_unit_cell import NCDist
from dxtbx.util import format_float_with_standard_uncertainty

if TYPE_CHECKING:
    import matplotlib.axes


logger = logging.getLogger(__name__)


class Cluster:
    def __init__(
        self,
        crystal_symmetries: list[crystal.symmetry],
        lattice_ids: Optional[list[int]] = None,
        name="",
    ):
        assert lattice_ids is None or len(crystal_symmetries) == len(lattice_ids)
        self.crystal_symmetries = crystal_symmetries
        self.lattice_ids = lattice_ids
        self.name = name
        self.unit_cells = np.array(
            [cs.unit_cell().parameters() for cs in crystal_symmetries]
        )
        self.median_cell = np.median(self.unit_cells, axis=0).tolist()
        self.cell_std = np.std(self.unit_cells, axis=0).tolist()
        self.mean_cell = np.mean(self.unit_cells, axis=0).tolist()
        self.pg_composition = collections.Counter(
            cs.space_group().type().lookup_symbol() for cs in crystal_symmetries
        )

    def __len__(self):
        return len(self.crystal_symmetries)

    def __str__(self):
        cell_str = ", ".join(
            format_float_with_standard_uncertainty(v, e, minimum=1e-5)
            for (v, e) in zip(self.median_cell, self.cell_std)
        )
        return "\n".join(
            [
                f"Median cell: {cell_str} ",
            ]
        )


@dataclass
class ClusteringResult:
    clusters: list[Cluster]
    dendrogram: Optional[dict] = None

    def __len__(self):
        return len(self.clusters)

    def __str__(self):
        text = [
            "\n{:<16} {:<8} {:<13} {:<13} {:<13} {:<12} {:<12} {:<12}{:<8}".format(
                "Cluster_id",
                "N_xtals",
                "Med_a",
                "Med_b",
                "Med_c",
                "Med_alpha",
                "Med_beta",
                "Med_gamma",
                "Delta(deg)",
            )
        ]
        singletons = []
        for cluster in self.clusters:
            if len(cluster) > 1:
                input_symmetry = crystal.symmetry(
                    unit_cell=uctbx.unit_cell(cluster.median_cell),
                    space_group_symbol="P 1",
                )
                groups = metric_subgroups(
                    input_symmetry, 3.00, enforce_max_delta_for_generated_two_folds=True
                )
                group = groups.result_groups[0]
                uc_params_conv = group["best_subsym"].unit_cell().parameters()
                sorted_pg_comp = sorted(
                    cluster.pg_composition.items(), key=lambda x: -1 * x[1]
                )
                pg_strings = ["{} in {}".format(pg[1], pg[0]) for pg in sorted_pg_comp]
                point_group_string = ", ".join(pg_strings) + "."
                text.append(point_group_string)
                text.append(
                    (
                        "{:<16} {:<8} {:<6.2f}({:<5.2f}) {:<6.2f}({:<5.2f})"
                        " {:<6.2f}({:<5.2f}) {:<6.2f}({:<4.2f}) {:<6.2f}"
                        "({:<4.2f}) {:<6.2f}({:<4.2f})"
                    ).format(
                        cluster.name,
                        len(cluster),
                        cluster.median_cell[0],
                        cluster.cell_std[0],
                        cluster.median_cell[1],
                        cluster.cell_std[1],
                        cluster.median_cell[2],
                        cluster.cell_std[2],
                        cluster.median_cell[3],
                        cluster.cell_std[3],
                        cluster.median_cell[4],
                        cluster.cell_std[4],
                        cluster.median_cell[5],
                        cluster.cell_std[5],
                    )
                )
                text.append(
                    (
                        "{:>24}  {:<6.2f}{:<7} {:<6.2f}{:<7}"
                        " {:<6.2f}{:<7} {:<6.2f}{:<6} {:<6.2f}"
                        "{:<6} {:<6.2f}{:<6}  {:<6.2}"
                    ).format(
                        group["best_subsym"].space_group_info().symbol_and_number(),
                        uc_params_conv[0],
                        "",
                        uc_params_conv[1],
                        "",
                        uc_params_conv[2],
                        "",
                        uc_params_conv[3],
                        "",
                        uc_params_conv[4],
                        "",
                        uc_params_conv[5],
                        "",
                        group["max_angular_difference"],
                    )
                )
            else:
                singletons.append(
                    "".join(
                        [
                            (
                                "{:<14} {:<11.2f} {:<11.2f} {:<11.2f}"
                                "{:<12.1f} {:<12.1f} {:<12.1f}"
                            ).format(
                                list(cluster.pg_composition.keys())[0],
                                cluster.median_cell[0],
                                cluster.median_cell[1],
                                cluster.median_cell[2],
                                cluster.median_cell[3],
                                cluster.median_cell[4],
                                cluster.median_cell[5],
                            ),
                        ]
                    )
                )
        text.append("Standard deviations are in brackets.")
        text.append(
            """Each cluster:
        Input lattice count, with integration Bravais setting space group.
        Cluster median with Niggli cell parameters (std dev in brackets).
        Highest possible metric symmetry and unit cell using LePage (J Appl Cryst 1982, 15:255) method, maximum delta 3deg."""
        )
        n_clusters = len(self.clusters) - len(singletons)
        text.insert(0, f"\n{n_clusters} cluster{'s'[:n_clusters^1]}")
        n_singletons = len(singletons)
        text[:0] = [
            f"{n_singletons} singleton{'s'[:n_singletons^1]}:\n",
            "{:<14} {:<11} {:<11} {:<11}{:<12} {:<12} {:<12}".format(
                "Point group", "a", "b", "c", "alpha", "beta", "gamma"
            ),
            "".join(singletons),
        ]
        return "\n".join(text)


def cluster_unit_cells(
    crystal_symmetries: list[crystal.symmetry],
    lattice_ids: Optional[list[int]] = None,
    threshold: int = 10000,
    ax: Optional["matplotlib.axes.Axes"] = None,
    no_plot: bool = True,
) -> Optional[ClusteringResult]:
    if not lattice_ids:
        lattice_ids = list(range(len(crystal_symmetries)))
    cluster = Cluster(crystal_symmetries, lattice_ids)
    uc = cluster.unit_cells

    a = uc[:, 0] ** 2
    b = uc[:, 1] ** 2
    c = uc[:, 2] ** 2
    d = 2 * uc[:, 1] * uc[:, 2] * np.cos(np.radians(uc[:, 3]))
    e = 2 * uc[:, 0] * uc[:, 2] * np.cos(np.radians(uc[:, 4]))
    f = 2 * uc[:, 0] * uc[:, 1] * np.cos(np.radians(uc[:, 5]))
    g6_cells = np.array([a, b, c, d, e, f]).transpose()

    logger.info(
        "Using Andrews-Bernstein distance from Andrews & Bernstein "
        "J Appl Cryst 47:346 (2014)"
    )
    metric = NCDist
    pair_distances = ssd.pdist(g6_cells, metric=metric)
    if len(pair_distances) > 0:
        logger.info("Distances have been calculated")
        this_linkage = hierarchy.linkage(pair_distances, method="single", metric=metric)
        cluster_ids = hierarchy.fcluster(this_linkage, threshold, criterion="distance")
        logger.debug("Clusters have been calculated")
    else:
        logger.debug("No distances were calculated. Aborting clustering.")
        return None

    # Create an array of sub-cluster objects from the clustering
    sub_clusters: list[Cluster] = []
    for i_cluster in set(cluster_ids):
        sel = np.where(cluster_ids == i_cluster)[0]
        cluster = Cluster(
            [crystal_symmetries[i] for i in sel],
            [lattice_ids[int(i)] for i in sel],
        )
        sub_clusters.append(cluster)

    # Order clusters by size
    sub_clusters = sorted(sub_clusters, key=len)
    for i, cluster in enumerate(sub_clusters):
        cluster.name = f"cluster_{i + 1}"

    dendrogram = hierarchy.dendrogram(
        this_linkage,
        # labels=labels,
        p=200,
        truncate_mode="lastp",  # show only the last p merged clusters
        leaf_font_size=8,
        leaf_rotation=90.0,
        color_threshold=threshold,
        ax=ax,
        no_plot=no_plot,
    )

    return ClusteringResult(clusters=sub_clusters, dendrogram=dendrogram)
