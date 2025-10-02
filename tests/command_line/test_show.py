from __future__ import annotations

import os
import shutil
import subprocess

from dxtbx.serialize import load

from dials.command_line.show import model_connectivity, run


def test_dials_show(dials_data):
    path = str(dials_data("experiment_test_data", pathlib=True) / "experiment_1.json")
    result = subprocess.run(
        [shutil.which("dials.show"), path],
        env={"DIALS_NOBANNER": "1", **os.environ},
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    output = result.stdout.decode("latin-1")
    output = [_f for _f in (s.rstrip() for s in output.split("\n")) if _f]
    assert (
        "\n".join(output[6:])
        == """
Detector:
Panel:
  name: Panel
  type: SENSOR_PAD
  identifier:
  pixel_size:{0.172,0.172}
  image_size: {2463,2527}
  trusted_range: {-1,495976}
  thickness: 0
  material:
  mu: 0
  gain: 1
  pedestal: 0
  fast_axis: {1,0,0}
  slow_axis: {0,-1,0}
  origin: {-212.478,220.002,-190.18}
  distance: 190.18
  pixel to millimeter strategy: SimplePxMmStrategy
Max resolution (at corners): 1.008178
Max resolution (inscribed):  1.204283
Beam:
    probe: x-ray
    wavelength: 0.9795
    sample to source direction : {0,0,1}
    divergence: 0
    sigma divergence: 0
    polarization normal: {0,1,0}
    polarization fraction: 0.999
    flux: 0
    transmission: 1
    sample to source distance: 0
Beam centre:
    mm: (212.48,220.00)
    px: (1235.34,1279.08)
Scan:
    number of images:   9
    image range:   {1,9}
    epoch:    1.36032e+09
    exposure time:    0.2
    oscillation:   {0,0.2}
Goniometer:
    Rotation axis:   {1,0,0}
    Fixed rotation:  {1,0,0,0,1,0,0,0,1}
    Setting rotation:{1,0,0,0,1,0,0,0,1}
Crystal:
    Unit cell: 42.272, 42.272, 39.670, 90.000, 89.999, 90.000
    Space group: P 4 2 2
    U matrix:  {{ 0.833603, -0.535987, -0.133506},
                {-0.179804, -0.034770, -0.983088},
                { 0.522280,  0.843510, -0.125357}}
    B matrix:  {{ 0.023656,  0.000000,  0.000000},
                {-0.000000,  0.023656,  0.000000},
                {-0.000000,  0.000000,  0.025208}}
    A = UB:    {{ 0.019720, -0.012679, -0.003365},
                {-0.004253, -0.000823, -0.024781},
                { 0.012355,  0.019954, -0.003160}}
    Mosaicity:  0.157000
""".strip()
    )


def test_dials_show_i04_weak_data(dials_data):
    data_dir = dials_data("i04_weak_data", pathlib=True)
    path = data_dir / "experiments_import.json"
    result = subprocess.run(
        [shutil.which("dials.show"), path],
        env={"DIALS_NOBANNER": "1", **os.environ},
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    output = result.stdout.decode("latin-1")
    output = [_f for _f in (s.rstrip() for s in output.split("\n")) if _f]
    assert (
        "\n".join(output[6:])
        == """
Detector:
Panel:
  name: Panel
  type: SENSOR_PAD
  identifier:
  pixel_size:{0.172,0.172}
  image_size: {2463,2527}
  trusted_range: {-1,161977}
  thickness: 0
  material:
  mu: 0
  gain: 1
  pedestal: 0
  fast_axis: {1,0,0}
  slow_axis: {0,-1,0}
  origin: {-210.76,205.277,-265.27}
  distance: 265.27
  pixel to millimeter strategy: SimplePxMmStrategy
Max resolution (at corners): 1.161261
Max resolution (inscribed):  1.509475
Beam:
    probe: x-ray
    wavelength: 0.97625
    sample to source direction : {0,0,1}
    divergence: 0
    sigma divergence: 0
    polarization normal: {0,1,0}
    polarization fraction: 0.999
    flux: 0
    transmission: 1
    sample to source distance: 0
Beam centre:
    mm: (210.76,205.28)
    px: (1225.35,1193.47)
Scan:
    number of images:   540
    image range:   {1,540}
    epoch:    1.37284e+09
    exposure time:    0.067
    oscillation:   {82,0.15}
Goniometer:
    Rotation axis:   {1,0,0}
    Fixed rotation:  {1,0,0,0,1,0,0,0,1}
    Setting rotation:{1,0,0,0,1,0,0,0,1}
""".strip()
    )


def test_dials_show_centroid_test_data(dials_data):
    result = subprocess.run(
        [shutil.which("dials.show")]
        + sorted(dials_data("centroid_test_data", pathlib=True).glob("centroid_*.cbf")),
        env={"DIALS_NOBANNER": "1", **os.environ},
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    output = result.stdout.decode("latin-1")
    output = [_f for _f in (s.rstrip() for s in output.split("\n")) if _f]
    assert (
        "\n".join(output[7:])
        == """
Detector:
Panel:
  name: Panel
  type: SENSOR_PAD
  identifier:
  pixel_size:{0.172,0.172}
  image_size: {2463,2527}
  trusted_range: {0,495976}
  thickness: 0.32
  material: Si
  mu: 3.96039
  gain: 1
  pedestal: 0
  fast_axis: {1,0,0}
  slow_axis: {0,-1,0}
  origin: {-212.478,220.002,-190.18}
  distance: 190.18
  pixel to millimeter strategy: ParallaxCorrectedPxMmStrategy
    mu: 3.96039
    t0: 0.32
Max resolution (at corners): 1.008375
Max resolution (inscribed):  1.204621
Beam:
    probe: x-ray
    wavelength: 0.9795
    sample to source direction : {0,0,1}
    divergence: 0
    sigma divergence: 0
    polarization normal: {0,1,0}
    polarization fraction: 0.999
    flux: 0
    transmission: 1
    sample to source distance: 0
Beam centre:
    mm: (212.48,220.00)
    px: (1235.34,1279.08)
Scan:
    number of images:   9
    image range:   {1,9}
    epoch:    1.36032e+09
    exposure time:    0.2
    oscillation:   {0,0.2}
Goniometer:
    Rotation axis:   {1,0,0}
    Fixed rotation:  {1,0,0,0,1,0,0,0,1}
    Setting rotation:{1,0,0,0,1,0,0,0,1}
""".strip()
    )


def test_dials_show_multi_panel_i23(dials_data):
    data_dir = dials_data("image_examples", pathlib=True)
    path = data_dir / "DLS_I23_germ_13KeV_0001.cbf"

    result = subprocess.run(
        [shutil.which("dials.show"), path],
        env={"DIALS_NOBANNER": "1", **os.environ},
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    output = result.stdout.decode("latin-1")
    output = [_f for _f in (s.rstrip() for s in output.split("\n")) if _f]
    assert (
        "\n".join(output[7:27])
        == """
Detector:
Panel:
  name: row-00
  type: SENSOR_PAD
  identifier:
  pixel_size:{0.172,0.172}
  image_size: {2463,195}
  trusted_range: {0,1e+06}
  thickness: 0.32
  material: Si
  mu: 3.663
  gain: 1
  pedestal: 0
  fast_axis: {-1,0,0}
  slow_axis: {0,-0.143467,-0.989655}
  origin: {191.952,-243.628,52.4929}
  distance: 248.638
  pixel to millimeter strategy: ParallaxCorrectedPxMmStrategy
    mu: 3.663
    t0: 0.32
""".strip()
    )

    assert (
        "\n".join(output[-47:])
        == """
Panel:
  name: row-23
  type: SENSOR_PAD
  identifier:
  pixel_size:{0.172,0.172}
  image_size: {2463,195}
  trusted_range: {0,1e+06}
  thickness: 0.32
  material: Si
  mu: 3.663
  gain: 1
  pedestal: 0
  fast_axis: {-1,0,0}
  slow_axis: {-0,-0.0638966,0.997957}
  origin: {191.952,251.939,-0.791525}
  distance: 251.373
  pixel to millimeter strategy: ParallaxCorrectedPxMmStrategy
    mu: 3.663
    t0: 0.32
Max resolution (at corners): 0.624307
Max resolution (inscribed):  0.829324
Beam:
    probe: x-ray
    wavelength: 0.95373
    sample to source direction : {0,0,1}
    divergence: 0
    sigma divergence: 0
    polarization normal: {0,1,0}
    polarization fraction: 0.999
    flux: 0
    transmission: 1
    sample to source distance: 0
Beam centre:
    mm: panel 12, (191.95,7.22)
    px: panel 12, (1116.00,41.96)
    mm, raw image: (191.95,444.63)
    px, raw image: (1116.00,2585.96)
Scan:
    number of images:   1
    image range:   {1,1}
    epoch:    1.41893e+09
    exposure time:    0.2
    oscillation:   {0,0.1}
Goniometer:
    Rotation axis:   {-1,0,0}
    Fixed rotation:  {1,0,0,0,1,0,0,0,1}
    Setting rotation:{1,0,0,0,1,0,0,0,1}
""".strip()
    )


def test_dials_show_reflection_table(dials_data):
    """Test the output of dials.show on a reflection_table pickle file"""
    result = subprocess.run(
        [
            shutil.which("dials.show"),
            dials_data("centroid_test_data", pathlib=True) / "integrated.pickle",
        ],
        env={"DIALS_NOBANNER": "1", **os.environ},
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    output = result.stdout.decode("latin-1")
    output = [_f for _f in (s.rstrip() for s in output.split("\n")) if _f]

    assert output[4] == "Reflection list contains 2269 reflections"
    headers = ["Column", "min", "max", "mean"]
    for header in headers:
        assert header in output[6]
    row_names = [
        "background.mean",
        "background.sum.value",
        "background.sum.variance",
        "d",
        "dqe",
        "flags",
        "id",
        "imageset_id",
        "intensity.prf.value",
        "intensity.prf.variance",
        "intensity.sum.value",
        "intensity.sum.variance",
        "lp",
        "miller_index",
        "num_pixels.background",
        "num_pixels.background_used",
        "num_pixels.foreground",
        "num_pixels.valid",
        "panel",
        "partial_id",
        "partiality",
        "profile.correlation",
        "profile.rmsd",
        "rlp",
        "s1",
        "shoebox",
        "summed I",
        "N pix",
        "N valid foreground pix",
        "xyzcal.mm",
        "xyzcal.px",
        "xyzobs.mm.value",
        "xyzobs.mm.variance",
        "xyzobs.px.value",
        "xyzobs.px.variance",
        "zeta",
    ]
    for name, out in zip(row_names, output[8:-1]):
        assert name in out


def test_dials_show_image_statistics(dials_data):
    # Run on one multi-panel image
    data_dir = dials_data("image_examples", pathlib=True)
    path = data_dir / "DLS_I23_germ_13KeV_0001.cbf"
    result = subprocess.run(
        [shutil.which("dials.show"), "image_statistics.show_raw=true", path],
        env={"DIALS_NOBANNER": "1", **os.environ},
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    output = result.stdout.decode("latin-1")
    output = [_f for _f in (s.rstrip() for s in output.split("\n")) if _f]
    assert (
        output[-1]
        == "DLS_I23_germ_13KeV_0001.cbf: Min: -2.0 Q1: 9.0 Med: 12.0 Q3: 16.0 Max: 1070079.0"
    )


def test_dials_show_on_scaled_data(dials_data):
    """Test that dials.show works on scaled data."""
    location = dials_data("l_cysteine_4_sweeps_scaled", pathlib=True)
    refl = location / "scaled_30.refl"
    expt = location / "scaled_30.expt"

    result = subprocess.run(
        [shutil.which("dials.show"), refl, expt], capture_output=True
    )
    assert not result.returncode and not result.stderr


def test_model_connectivity(dials_data):
    """Test that dials.show experiments_has_model option."""
    location = dials_data("l_cysteine_dials_output", pathlib=True)
    expts = load.experiment_list(location / "indexed.expt", check_format=False)
    assert (
        model_connectivity(expts)
        == """\
Experiment / Models

Detector:
              0  1
Experiment 0  x  .
Experiment 1  x  .
Experiment 2  x  .
Experiment 3  .  x

Crystal:
              0
Experiment 0  x
Experiment 1  x
Experiment 2  x
Experiment 3  x

Beam:
              0
Experiment 0  x
Experiment 1  x
Experiment 2  x
Experiment 3  x"""
    )


def test_dials_show_shared_models(dials_data, capsys):
    """Test that dials.show experiments_has_model option."""
    location = dials_data("l_cysteine_dials_output", pathlib=True)
    run([str(location / "indexed.expt"), "show_shared_models=True"])
    stdout, stderr = capsys.readouterr()
    assert not stderr
    assert "Experiment / Models" in stdout


def test_dials_show_centroid_test_data_image_zero(dials_data, tmp_path):
    """Integration test: import image 0; show import / show works"""

    im1 = dials_data("centroid_test_data", pathlib=True) / "centroid_0001.cbf"
    im0 = tmp_path / "centroid_0000.cbf"

    shutil.copyfile(im1, im0)

    result = subprocess.run(
        [shutil.which("dials.import"), im0], cwd=tmp_path, capture_output=True
    )
    assert not result.returncode and not result.stderr

    result = subprocess.run(
        [shutil.which("dials.show"), "imported.expt"], cwd=tmp_path, capture_output=True
    )
    assert not result.returncode and not result.stderr
