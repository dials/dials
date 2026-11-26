"""Tests for scaling time-of-flight data with the tof scaling model."""

from __future__ import annotations

import math

import numpy as np
import pytest

from cctbx import miller, sgtbx
from dxtbx.model import Experiment, ExperimentList
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.util import ersatz_uuid4
from scitbx import matrix

from dials.algorithms.scaling.algorithm import ScalingAlgorithm
from dials.algorithms.scaling.scaling_utilities import (
    align_axis_along_z,
)
from dials.algorithms.scaling.tof.components import PanelScaleComponent
from dials.algorithms.scaling.tof.model import (
    TOFScalingModel,
    check_panel_scales,
    plot_panel_scales,
)
from dials.algorithms.scaling.tof.preparation import (
    LAMBDA_NORM_COLUMN,
    add_tof_crystal_frame_vectors,
    prepare_tof_input,
)
from dials.array_family import flex
from dials.command_line.scale import phil_scope
from dials.report.analysis import reflection_tables_to_batch_dependent_properties
from dials.util.batch_handling import (
    calculate_batch_offsets,
    get_batch_ranges,
)


def _tof_experiments(dials_data, n_runs):
    """ToF experiments, one per run, each with the crystal rotated by a
    different goniometer setting rotation."""
    expt_file = dials_data("isis_sxd_nacl_processed") / "refined.expt"
    experiments = ExperimentList()
    for i in range(n_runs):
        expt = ExperimentListFactory.from_json_file(expt_file, check_format=False)[0]
        expt.identifier = ersatz_uuid4()
        expt.crystal.set_space_group(sgtbx.space_group())
        angle = 2.0 * math.pi * i / n_runs
        setting = matrix.col((0, 1, 0)).axis_and_angle_as_r3_rotation_matrix(angle)
        expt.goniometer.set_setting_rotation(setting)
        experiments.append(expt)
    return experiments


# Simulated relative scales of the 11 SXD panels
PANEL_SCALES = [1.0, 0.8, 1.2, 0.9, 1.1, 1.0, 0.7, 1.3, 1.0, 0.95, 1.05]


def _tof_reflections(experiments, run_scales, seed=0, panel_scales=None):
    """Simulated observations of the same reflections in each run, with each
    run's intensities multiplied by its scale (and by the scale of the panel
    each was recorded on, if given)."""
    rng = np.random.default_rng(seed)
    crystal_symmetry = experiments[0].crystal.get_crystal_symmetry()
    indices = (
        miller.build_set(crystal_symmetry, anomalous_flag=True, d_min=0.6)
        .expand_to_p1()
        .indices()
    )
    true_intensities = rng.exponential(1000.0, len(indices))
    unit_cell = experiments[0].crystal.get_unit_cell()
    tables = []
    for i, (expt, k) in enumerate(zip(experiments, run_scales)):
        sel = rng.random(len(indices)) < 0.7
        n = int(sel.sum())
        table = flex.reflection_table()
        table["miller_index"] = flex.miller_index(
            [h for h, s in zip(indices, sel) if s]
        )
        table["d"] = unit_cell.d(table["miller_index"])
        panels = rng.integers(0, len(expt.detector), n)
        intensity = k * true_intensities[sel]
        if panel_scales is not None:
            intensity *= np.array(panel_scales)[panels]
        variance = intensity + 10.0
        intensity += rng.normal(0.0, np.sqrt(variance))
        table["intensity.sum.value"] = flex.double(intensity)
        table["intensity.sum.variance"] = flex.double(variance)
        wavelength = rng.uniform(0.5, 5.0, n)
        directions = rng.normal(size=(n, 3))
        directions /= np.linalg.norm(directions, axis=1)[:, np.newaxis]
        table["s1"] = flex.vec3_double(directions / wavelength[:, np.newaxis])
        table["wavelength_cal"] = flex.double(wavelength)
        table["partiality"] = flex.double(n, 1.0)
        table[LAMBDA_NORM_COLUMN] = flex.double(n, 1.0)
        table["xyzobs.px.value"] = flex.vec3_double(rng.uniform(1, 1800, size=(n, 3)))
        table["panel"] = flex.size_t(panels.tolist())
        table["id"] = flex.int(n, i)
        table.experiment_identifiers()[i] = expt.identifier
        table.set_flags(flex.bool(n, True), table.flags.integrated_sum)
        tables.append(table)
    return tables


def _params(*options):
    from libtbx.phil import parse

    params = phil_scope.fetch(parse("\n".join(options))).extract()
    params.output.html = None
    params.output.json = None
    return params


def test_model_from_data_and_dict(dials_data):
    params = _params(
        "tof.decay_correction=True",
        "tof.absorption_correction=True",
        "tof.reference_wavelength=1.5",
    )
    (expt,) = _tof_experiments(dials_data, 1)
    model = TOFScalingModel.from_data(params, expt, None)
    assert list(model.components) == ["scale", "decay", "absorption"]
    assert model.components["absorption"].n_params == 24
    assert model.get_shared_components() == ("absorption",)
    assert model.configdict["reference_wavelength"] == 1.5

    copy = TOFScalingModel.from_dict(model.to_dict())
    assert copy.is_scaled
    assert list(copy.components) == list(model.components)
    assert copy.configdict == model.configdict

    params = _params(
        "tof.panel_correction=True", "tof.panel_reference=2", "tof.panel_fixed=5"
    )
    model = TOFScalingModel.from_data(params, expt, None)
    assert list(model.components) == ["scale", "panel"]
    assert model.components["panel"].n_params == 11
    assert model.components["panel"].fixed_panels == [2, 5]
    assert len(model.components["panel"].free_parameters) == 9
    assert model.get_shared_components() == ("panel",)
    copy = TOFScalingModel.from_dict(model.to_dict())
    assert copy.components["panel"].fixed_panels == [2, 5]
    assert copy.configdict == model.configdict

    params = _params(
        "tof.panel_correction=True",
        "tof.panel_reference=0",
        "tof.absorption_correction=True",
        "tof.reference_wavelength=1.5",
    )
    model = TOFScalingModel.from_data(params, expt, None)
    assert list(model.components) == ["scale", "absorption", "panel"]
    assert model.get_shared_components() == ("absorption", "panel")

    model = TOFScalingModel.from_data(_params(), None, None)
    assert list(model.components) == ["scale"]
    assert model.get_shared_components() is None

    params = _params("tof.absorption_correction=True")
    with pytest.raises(ValueError, match="reference wavelength is needed"):
        TOFScalingModel.from_data(params, None, None)


def test_crystal_frame_vectors(dials_data):
    (expt,) = _tof_experiments(dials_data, 1)
    R = matrix.col((0, 1, 0)).axis_and_angle_as_r3_rotation_matrix(0.3)
    expt.goniometer.set_setting_rotation(R)
    table = flex.reflection_table()
    table["s1"] = flex.vec3_double([(0.5, 0.0, 0.0), (0.0, 0.0, 0.0)])
    add_tof_crystal_frame_vectors(table, expt)

    R = matrix.sqr(R)
    s0 = matrix.col(expt.beam.get_sample_to_source_direction())
    expected_s0c = align_axis_along_z(
        (1.0, 0.0, 0.0), flex.vec3_double([R.inverse() * s0])
    )
    expected_s1c = align_axis_along_z(
        (1.0, 0.0, 0.0), flex.vec3_double([R.inverse() * matrix.col((1, 0, 0))])
    )
    assert table["s0c"][0] == pytest.approx(expected_s0c[0])
    assert table["s1c"][0] == pytest.approx(expected_s1c[0])
    # reflections without a diffracted beam vector are left at zero
    assert table["s1c"][1] == (0.0, 0.0, 0.0)


def test_absorption_scaled_by_wavelength(dials_data):
    experiments = _tof_experiments(dials_data, 1)
    (table,) = _tof_reflections(experiments, [1.0])
    add_tof_crystal_frame_vectors(table, experiments[0])
    tables = {}
    for dependence in ("linear", "none"):
        params = _params(
            "tof.absorption_correction=True",
            "tof.reference_wavelength=2.0",
            f"tof.absorption_wavelength_dependence={dependence}",
        )
        model = TOFScalingModel.from_data(params, experiments[0], table)
        model.configure_components(table, experiments[0], params)
        tables[dependence] = model.components["absorption"].data["sph_harm_table"]
    factors = table["wavelength_cal"] / 2.0
    for j in (0, 10, 100):
        linear = tables["linear"].col(j).as_dense_vector()
        none = tables["none"].col(j).as_dense_vector()
        assert list(linear) == pytest.approx(list(none * factors[j]))


def test_prepare_tof_input(dials_data):
    experiments = _tof_experiments(dials_data, 2)
    tables = _tof_reflections(experiments, [1.0, 1.0])

    params, _, _ = prepare_tof_input(_params(), experiments, tables)
    assert params.model == "tof"
    assert "s1c" not in tables[0]

    params, _, _ = prepare_tof_input(
        _params("tof.absorption_correction=True"), experiments, tables
    )
    wavelengths = np.concatenate([np.array(t["wavelength_cal"]) for t in tables])
    assert params.tof.reference_wavelength == pytest.approx(np.median(wavelengths))
    assert "s0c" in tables[0] and "s1c" in tables[1]


@pytest.mark.parametrize(
    "options,n_runs,match",
    [
        (["model=physical"], 2, "physical model cannot be used"),
        (["exclude_images=0:1:100"], 2, "exclude_images cannot be used"),
        (
            ["filtering.method=deltacchalf", "filtering.deltacchalf.mode=image_group"],
            2,
            "mode=image_group cannot be used",
        ),
        (["output.merged_mtz=merged.mtz"], 2, "merged_mtz is not supported"),
        ([], 1, "needs more than one run"),
    ],
)
def test_prepare_tof_input_errors(dials_data, options, n_runs, match):
    experiments = _tof_experiments(dials_data, n_runs)
    tables = _tof_reflections(experiments, [1.0] * n_runs)
    with pytest.raises(ValueError, match=match):
        prepare_tof_input(_params(*options), experiments, tables)


def test_prepare_tof_input_requires_lambda_norm(dials_data):
    experiments = _tof_experiments(dials_data, 2)
    tables = _tof_reflections(experiments, [1.0, 1.0])
    del tables[1][LAMBDA_NORM_COLUMN]
    with pytest.raises(ValueError, match="dials.lambda_norm"):
        prepare_tof_input(_params(), experiments, tables)
    prepare_tof_input(_params("tof.require_lambda_norm=False"), experiments, tables)


def test_prepare_tof_input_mixed_types(dials_data):
    experiments = _tof_experiments(dials_data, 1)
    experiments.append(Experiment())
    with pytest.raises(ValueError, match="cannot be scaled together"):
        prepare_tof_input(_params(), experiments, [None, None])


@pytest.mark.parametrize(
    "options",
    [
        [],
        ["tof.decay_correction=True", "tof.absorption_correction=True"],
        ["model=KB"],
    ],
)
def test_scale_tof_runs(dials_data, tmp_path, monkeypatch, options):
    """Scaling recovers the relative scale of each run."""
    monkeypatch.chdir(tmp_path)
    run_scales = [1.0, 1.3, 0.8, 1.1]
    experiments = _tof_experiments(dials_data, len(run_scales))
    tables = _tof_reflections(experiments, run_scales)
    algorithm = ScalingAlgorithm(_params(*options), experiments, tables)
    algorithm.run()

    models = [e.scaling_model for e in algorithm.experiments]
    expected_id = "KB" if "model=KB" in options else "tof"
    assert all(m.id_ == expected_id for m in models)
    scales = [m.components["scale"].parameters[0] for m in models]
    for scale, k in zip(scales, run_scales):
        assert scale / scales[0] == pytest.approx(k / run_scales[0], rel=0.02)
    if "tof.absorption_correction=True" in options:
        surfaces = [list(m.components["absorption"].parameters) for m in models]
        assert all(s == surfaces[0] for s in surfaces)
    assert algorithm.merging_statistics_result.overall.cc_one_half > 0.95

    # The report has one batch per run
    reflection_tables_to_batch_dependent_properties(
        algorithm.reflections, algorithm.experiments, algorithm.scaled_miller_array
    )
    offsets = calculate_batch_offsets(algorithm.experiments)
    batch_ranges = get_batch_ranges(algorithm.experiments, offsets)
    for table, batch_range in zip(algorithm.reflections, batch_ranges):
        assert batch_range[0] == batch_range[1]
        assert set(table["batch"]) == {batch_range[0]}


def test_panel_scale_component():
    component = PanelScaleComponent(flex.double([1.0, 2.0, 3.0, 4.0]), fixed_panels=[1])
    component.data = {"panel": flex.size_t([0, 1, 2, 3, 3])}
    component.update_reflection_data()
    assert list(component.free_parameters) == [1.0, 3.0, 4.0]
    scales, derivatives = component.calculate_scales_and_derivatives()
    assert list(scales) == [1.0, 2.0, 3.0, 4.0, 4.0]
    assert derivatives.n_cols == 3
    dense = [[derivatives[i, j] for j in range(3)] for i in range(5)]
    assert dense == [[1, 0, 0], [0, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 1]]

    component.free_parameters = flex.double([5.0, 6.0, 7.0])
    assert list(component.parameters) == [5.0, 2.0, 6.0, 7.0]
    component.free_parameter_esds = flex.double([0.1, 0.2, 0.3])
    assert list(component.parameter_esds) == [0.1, 0.0, 0.2, 0.3]

    component.update_reflection_data(block_selections=[flex.size_t([3, 1])])
    assert list(component.calculate_scales()) == [7.0, 2.0]


def test_prepare_panel_options(dials_data):
    experiments = _tof_experiments(dials_data, 2)
    tables = _tof_reflections(experiments, [1.0, 1.0])
    for table in tables:
        sel = (table["panel"] == 3) | (table["panel"] == 7)
        table["panel"].set_selected(sel, 4)
    params, _, _ = prepare_tof_input(
        _params("tof.panel_correction=True"), experiments, tables
    )
    assert params.tof.panel_reference == 4
    assert params.tof.panel_fixed == [3, 7]

    params = _params("tof.panel_correction=True", "tof.panel_reference=3")
    with pytest.raises(ValueError, match="reference panel 3 has no observations"):
        prepare_tof_input(params, experiments, tables)


@pytest.mark.parametrize("absorption", [False, True])
def test_scale_tof_panels(dials_data, tmp_path, monkeypatch, absorption):
    """Scaling recovers the relative scales of the runs and the panels, also
    with an absorption surface shared across runs as well as the panel scales."""
    monkeypatch.chdir(tmp_path)
    run_scales = [1.0, 1.3, 0.8]
    experiments = _tof_experiments(dials_data, len(run_scales))
    tables = _tof_reflections(experiments, run_scales, panel_scales=PANEL_SCALES)
    params = _params(
        "tof.panel_correction=True",
        "tof.panel_reference=0",
        f"tof.absorption_correction={absorption}",
    )
    algorithm = ScalingAlgorithm(params, experiments, tables)
    algorithm.run()

    models = [e.scaling_model for e in algorithm.experiments]
    panels = [list(m.components["panel"].parameters) for m in models]
    assert all(p == panels[0] for p in panels)
    if absorption:
        surfaces = [list(m.components["absorption"].parameters) for m in models]
        assert all(s == surfaces[0] for s in surfaces)
    assert panels[0][0] == 1.0
    assert panels[0] == pytest.approx(PANEL_SCALES, rel=0.05)
    scales = [m.components["scale"].parameters[0] for m in models]
    # Each run samples different crystal-frame directions, so a free absorption
    # surface can take up a little of the run scales.
    rel = 0.03 if absorption else 0.02
    for scale, k in zip(scales, run_scales):
        assert scale / scales[0] == pytest.approx(k / run_scales[0], rel=rel)
    plot = plot_panel_scales(models[0])["panel_scales"]
    assert plot["data"][0]["y"] == pytest.approx(panels[0], abs=1e-5)
    assert check_panel_scales(algorithm.experiments) == []
    models[0].components["panel"].parameters[3] = 0.1
    assert check_panel_scales(algorithm.experiments) == [3]
