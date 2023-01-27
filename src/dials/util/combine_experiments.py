from __future__ import annotations

import sys
from collections.abc import Sequence
from typing import Iterator, Tuple, TypeVar

import dxtbx.model
import dxtbx.model.compare as compare
from dxtbx.model.experiment_list import (
    BeamComparison,
    DetectorComparison,
    Experiment,
    ExperimentList,
    GoniometerComparison,
)
from libtbx import phil
from scitbx import matrix

from dials.algorithms.clustering.unit_cell import cluster_unit_cells

T = TypeVar("T")


def _split_equal_parts_of_length(a: Sequence[T], n: int) -> Iterator[Sequence[T]]:
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)] for i in range(n))


def find_experiment_in(experiment, all_experiments):
    """Search the phil experiment list and find where an experiment came from.

    :param Experiment experiment: The experiment to search for
    :param all_experiments:       The list of all experiments from phil
    :type  all_experiments:       list[dials.util.phil.FilenameDataWrapper[ExperimentList]]
    :returns:                     The filename and experiment ID
    :rtype:                       (str, int)
    """
    for source in all_experiments:
        try:
            experiment_list = list(source.data)
            index = experiment_list.index(experiment)
            return (source.filename, index)
        except ValueError:
            pass
    raise ValueError("Experiment not found")


class ComparisonError(Exception):
    """Exception to indicate problem with tolerance comparisons"""

    pass


class CombineWithReference:
    def __init__(
        self,
        beam=None,
        goniometer=None,
        scan=None,
        crystal=None,
        detector=None,
        params=None,
    ):

        self.ref_beam = beam
        self.ref_goniometer = goniometer
        self.ref_scan = scan
        self.ref_crystal = crystal
        self.ref_detector = detector
        self.tolerance = None
        self._last_imageset = None
        if params:
            if params.reference_from_experiment.compare_models:
                self.tolerance = params.reference_from_experiment.tolerance
            self.average_detector = params.reference_from_experiment.average_detector
        else:
            self.average_detector = False

    def __call__(self, experiment):
        if self.tolerance:
            compare_beam = BeamComparison(
                wavelength_tolerance=self.tolerance.beam.wavelength,
                direction_tolerance=self.tolerance.beam.direction,
                polarization_normal_tolerance=self.tolerance.beam.polarization_normal,
                polarization_fraction_tolerance=self.tolerance.beam.polarization_fraction,
            )
            compare_detector = DetectorComparison(
                fast_axis_tolerance=self.tolerance.detector.fast_axis,
                slow_axis_tolerance=self.tolerance.detector.slow_axis,
                origin_tolerance=self.tolerance.detector.origin,
            )
            compare_goniometer = GoniometerComparison(
                rotation_axis_tolerance=self.tolerance.goniometer.rotation_axis,
                fixed_rotation_tolerance=self.tolerance.goniometer.fixed_rotation,
                setting_rotation_tolerance=self.tolerance.goniometer.setting_rotation,
            )

        else:
            compare_beam = None
            compare_detector = None
            compare_goniometer = None

        if self.ref_beam:
            if compare_beam:
                if not compare_beam(self.ref_beam, experiment.beam):
                    raise ComparisonError(
                        compare.beam_diff(
                            self.ref_beam,
                            experiment.beam,
                            wavelength_tolerance=self.tolerance.beam.wavelength,
                            direction_tolerance=self.tolerance.beam.direction,
                            polarization_normal_tolerance=self.tolerance.beam.polarization_normal,
                            polarization_fraction_tolerance=self.tolerance.beam.polarization_fraction,
                        )
                    )
            beam = self.ref_beam
        else:
            beam = experiment.beam

        if self.ref_detector and self.average_detector:
            detector = self.ref_detector
        elif self.ref_detector and not self.average_detector:
            if compare_detector:
                if not compare_detector(self.ref_detector, experiment.detector):
                    raise ComparisonError(
                        compare.detector_diff(
                            self.ref_detector,
                            experiment.detector,
                            fast_axis_tolerance=self.tolerance.detector.fast_axis,
                            slow_axis_tolerance=self.tolerance.detector.slow_axis,
                            origin_tolerance=self.tolerance.detector.origin,
                        )
                    )
            detector = self.ref_detector
        else:
            detector = experiment.detector

        if self.ref_goniometer:
            if compare_goniometer:
                if not compare_goniometer(self.ref_goniometer, experiment.goniometer):
                    raise ComparisonError(
                        compare.goniometer_diff(
                            self.ref_goniometer,
                            experiment.goniometer,
                            rotation_axis_tolerance=self.tolerance.goniometer.rotation_axis,
                            fixed_rotation_tolerance=self.tolerance.goniometer.fixed_rotation,
                            setting_rotation_tolerance=self.tolerance.goniometer.setting_rotation,
                        )
                    )
            goniometer = self.ref_goniometer
        else:
            goniometer = experiment.goniometer

        if self.ref_scan:
            scan = self.ref_scan
        else:
            scan = experiment.scan

        if self.ref_crystal:
            crystal = self.ref_crystal
        else:
            crystal = experiment.crystal

        if self._last_imageset == experiment.imageset:
            imageset = self._last_imageset
        else:
            imageset = experiment.imageset
            self._last_imageset = imageset

        return Experiment(
            identifier=experiment.identifier,
            beam=beam,
            detector=detector,
            scan=scan,
            goniometer=goniometer,
            crystal=crystal,
            imageset=imageset,
        )


def do_unit_cell_clustering(experiments, reflections, dendrogram=False, threshold=1000):
    if dendrogram:
        import matplotlib.pyplot as plt

        ax = plt.gca()
    else:
        ax = None

    crystal_symmetries = [expt.crystal.get_crystal_symmetry() for expt in experiments]
    clustering = cluster_unit_cells(
        crystal_symmetries,
        lattice_ids=list(experiments.identifiers()),
        threshold=threshold,
        ax=ax,
        no_plot=not dendrogram,
    )

    if dendrogram:
        ax.set_yscale("symlog", linthresh=1)
        plt.tight_layout()
        plt.show()

    return clustering


def average_detectors(target, panelgroups, depth, average_hierarchy_level=None):
    # Recursive function to do the averaging

    if average_hierarchy_level is None or depth == average_hierarchy_level:
        n = len(panelgroups)
        sum_fast = matrix.col((0.0, 0.0, 0.0))
        sum_slow = matrix.col((0.0, 0.0, 0.0))
        sum_ori = matrix.col((0.0, 0.0, 0.0))

        # Average the d matrix vectors
        for pg in panelgroups:
            sum_fast += matrix.col(pg.get_local_fast_axis())
            sum_slow += matrix.col(pg.get_local_slow_axis())
            sum_ori += matrix.col(pg.get_local_origin())
        sum_fast /= n
        sum_slow /= n
        sum_ori /= n

        # Re-orthagonalize the slow and the fast vectors by rotating around the cross product
        c = sum_fast.cross(sum_slow)
        a = sum_fast.angle(sum_slow, deg=True) / 2
        sum_fast = sum_fast.rotate_around_origin(c, a - 45, deg=True)
        sum_slow = sum_slow.rotate_around_origin(c, -(a - 45), deg=True)

        target.set_local_frame(sum_fast, sum_slow, sum_ori)

    if target.is_group():
        # Recurse
        for i, target_pg in enumerate(target):
            average_detectors(
                target_pg,
                [pg[i] for pg in panelgroups],
                depth + 1,
                average_hierarchy_level,
            )


def parse_ref_models(
    flat_exps: ExperimentList, reference_from_experiment: phil.scope_extract
) -> Tuple[
    None | dxtbx.model.beam,
    None | dxtbx.model.goniometer,
    None | dxtbx.model.scan,
    None | dxtbx.model.crystal,
    None | dxtbx.model.detector,
]:
    ref_beam = reference_from_experiment.beam
    ref_goniometer = reference_from_experiment.goniometer
    ref_scan = reference_from_experiment.scan
    ref_crystal = reference_from_experiment.crystal
    ref_detector = reference_from_experiment.detector

    if ref_beam is not None:
        try:
            ref_beam = flat_exps[ref_beam].beam
        except IndexError:
            sys.exit(f"{ref_beam} is not a valid experiment ID")

    if ref_goniometer is not None:
        try:
            ref_goniometer = flat_exps[ref_goniometer].goniometer
        except IndexError:
            sys.exit(f"{ref_goniometer} is not a valid experiment ID")

    if ref_scan is not None:
        try:
            ref_scan = flat_exps[ref_scan].scan
        except IndexError:
            sys.exit(f"{ref_scan} is not a valid experiment ID")

    if ref_crystal is not None:
        try:
            ref_crystal = flat_exps[ref_crystal].crystal
        except IndexError:
            sys.exit(f"{ref_crystal} is not a valid experiment ID")

    if ref_detector is not None:
        assert not reference_from_experiment.average_detector
        try:
            ref_detector = flat_exps[ref_detector].detector
        except IndexError:
            sys.exit(f"{ref_detector} is not a valid experiment ID")
    elif reference_from_experiment.average_detector:
        # Average all of the detectors together
        ref_detector = flat_exps[0].detector
        average_detectors(
            ref_detector.hierarchy(),
            [e.detector.hierarchy() for e in flat_exps],
            0,
            reference_from_experiment.average_hierarchy_level,
        )
    return (ref_beam, ref_goniometer, ref_scan, ref_crystal, ref_detector)
