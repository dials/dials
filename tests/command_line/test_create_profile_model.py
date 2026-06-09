from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest


def test_basic(dials_data, tmp_path):
    data_dir = dials_data("insulin_processed")
    # Call dials.create_profile_model
    result = subprocess.run(
        [
            shutil.which("dials.create_profile_model"),
            data_dir / "indexed.expt",
            data_dir / "indexed.refl",
            "sigma_m_algorithm=basic",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert os.path.exists(tmp_path / "models_with_profiles.expt")

    from dxtbx.model.experiment_list import ExperimentListFactory

    experiments = ExperimentListFactory.from_json_file(
        tmp_path / "models_with_profiles.expt", check_format=False
    )
    sigma_b = experiments[0].profile.sigma_b(deg=True)
    sigma_m = experiments[0].profile.sigma_m(deg=True)
    assert sigma_b == pytest.approx(0.0539, abs=1e-3)
    assert sigma_m == pytest.approx(0.1873, abs=1e-3)


def test_extended(dials_data: Path, tmp_path):
    data_dir = dials_data("insulin_processed")
    # Call dials.create_profile_model
    result = subprocess.run(
        [
            shutil.which("dials.create_profile_model"),
            data_dir / "indexed.expt",
            data_dir / "indexed.refl",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr
    assert os.path.exists(tmp_path / "models_with_profiles.expt")

    from dxtbx.model.experiment_list import ExperimentListFactory

    experiments = ExperimentListFactory.from_json_file(
        tmp_path / "models_with_profiles.expt", check_format=False
    )
    sigma_b = experiments[0].profile.sigma_b(deg=True)
    sigma_m = experiments[0].profile.sigma_m(deg=True)
    assert sigma_b == pytest.approx(0.0539, abs=1e-3)
    assert sigma_m == pytest.approx(0.1584, abs=1e-3)
