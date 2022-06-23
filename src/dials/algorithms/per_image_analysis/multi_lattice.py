from __future__ import annotations

import math
import time

import matplotlib.pyplot as plt
import numpy as np
import peakutils
from matplotlib import ticker
from numpy.polynomial import Polynomial
from sklearn.neighbors import NearestNeighbors

from dxtbx import flumpy
from dxtbx.model import ExperimentList

from dials.array_family import flex


@ticker.FuncFormatter
def ddv_formatter(d, pos):
    d = 1 / d if d > 0 else 0
    return f"{d:.2f}"


def compute_ddv_histogram(
    reflections: flex.reflection_table, bins: int = 100
) -> tuple(np.ndarray, np.ndarray):
    distances = np.empty((0,))
    k = 200

    for i in set(reflections["id"]):
        refl = reflections.select(reflections["id"] == i)
        n_refl = len(refl)
        if n_refl <= 10:
            continue
        rs_vectors = flumpy.to_numpy(refl["rlp"])
        t0 = time.time()
        nbrs = NearestNeighbors(n_neighbors=min(k, n_refl), algorithm="kd_tree").fit(
            rs_vectors
        )
        dist, indices = nbrs.kneighbors(rs_vectors)
        t1 = time.time()
        print(f"Nearest neighbour search took {t1 - t0:.4f} seconds")
        distances = np.concatenate((distances, dist.flatten()))
    distances = distances[(distances > 0).nonzero()]
    distances = distances[((distances >= 0.001) & (distances <= 0.04)).nonzero()]
    hist, bin_edges = np.histogram(distances, bins=bins)
    return hist, bin_edges


def compute_baseline_fit(hist: np.ndarray, bin_edges: np.ndarray) -> Polynomial:
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2

    base = peakutils.baseline(hist, deg=1)
    fit = Polynomial.fit(bin_centers, base, deg=1)
    print(f"{fit=}")
    return fit


def compute_spherical_cap(
    experiments: ExperimentList, reflections: flex.reflection_table
) -> float:

    max_two_theta = 0

    for i, expt in enumerate(experiments):
        refl = reflections.select(reflections["id"] == i)
        if len(refl) == 0:
            continue
        s1 = flumpy.to_numpy(refl["s1"])
        one_over_lambda = 1 / expt.beam.get_wavelength()
        s0 = np.array(expt.beam.get_s0())
        cos_two_theta = s1.dot(s0) / one_over_lambda**2
        two_theta = np.arccos(cos_two_theta)
        max_two_theta = max(two_theta.max(), max_two_theta)

    spherical_cap_area = (
        2 * math.pi * one_over_lambda**2 * (1 - math.cos(max_two_theta))
    )
    N = len(reflections)
    print(f"max(2θ): {np.degrees(max_two_theta)}°")
    print(f"{spherical_cap_area=} Å⁻²")
    print(f"N refl: {N}")
    return spherical_cap_area


def compute_K(experiments: ExperimentList, reflections: flex.reflection_table) -> float:
    if len(reflections) <= 10:
        return None

    hist, bin_edges = compute_ddv_histogram(reflections)
    if np.all(hist == 0):
        return None
    print(hist)
    baseline_fit = compute_baseline_fit(hist, bin_edges)
    k0 = baseline_fit.coef[0]
    spherical_cap_area = compute_spherical_cap(experiments, reflections)
    N = len(reflections)

    K = k0 * spherical_cap_area / N**2
    print(f"{k0=}, {K=}")
    print("K = k0 * S / N ** 2")
    print(f"K = {k0} Å * {spherical_cap_area} Å⁻² / {N} ** 2")
    print(f"{K=} Å⁻¹")
    # indexes = peakutils.indexes(hist - base, thres=0.5, min_dist=10)

    # print(indexes)

    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    baseline = baseline_fit(bin_centers)

    if 1:
        plt.plot(bin_centers, hist, label="Raw")
        # plt.scatter(bin_centers[indexes], hist[indexes], marker="+", color="red")
        # plt.scatter(bin_centers[indexes], (hist - baseline)[indexes], marker="+", color="red")
        plt.plot(bin_centers, hist - baseline, label="Baseline-subtracted")
        plt.plot(bin_centers, baseline, label="Baseline")
        plt.legend()

        # plt.plot(bin_centers, hist)
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ddv_formatter)
        plt.xlim((0.001, 0.04))
        plt.show()

    return K


import sys

from dials.util import log
from dials.util.options import ArgumentParser, flatten_experiments, flatten_reflections


def run(args=None):
    usage = (
        "dev.dials.nearest_neighbour_analysis [options] datablock.json strong.pickle"
    )

    parser = ArgumentParser(
        usage=usage,
        # phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        # epilog=help_message,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging
    log.config(
        options.verbose,
    )

    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    if len(reflections) == 0:
        raise sys.exit("No reflection lists found in input")
    if len(reflections) > 1:
        assert len(reflections) == len(experiments)
        for i in range(len(reflections)):
            reflections[i]["imageset_id"] = flex.int(len(reflections[i]), i)
            if i > 0:
                reflections[0].extend(reflections[i])

    reflections = reflections[0]
    if "imageset_id" not in reflections:
        reflections["imageset_id"] = reflections["id"]

    for expt in experiments:
        if (
            expt.goniometer is not None
            and expt.scan is not None
            and expt.scan.get_oscillation()[1] == 0
        ):
            expt.goniometer = None
            expt.scan = None

    reflections.centroid_px_to_mm(experiments)
    reflections.map_centroids_to_reciprocal_space(experiments)
    reflections.calculate_entering_flags(experiments)

    # reflections = reflections.select(reflections["id"] == 0)

    compute_K(experiments, reflections)


if __name__ == "__main__":
    run()
