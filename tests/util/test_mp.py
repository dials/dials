from __future__ import annotations

import os

import dials.util.mp


def test_available_cores():
    # We can't make any assumptions about the testing environment,
    # but we know there will be at least one available core, and
    # the function must return a positive integer in any case.
    assert dials.util.mp.available_cores() >= 1


def test_available_cores_nslots(monkeypatch):
    # For now at least the NSLOTS environment variable should
    # override whatever available_cores() returns otherwise.
    true_cores = dials.util.mp.available_cores()
    monkeypatch.setenv("NSLOTS", str(true_cores + 1))
    assert dials.util.mp.available_cores() == true_cores + 1


def test_available_cores_condor_job_ad(monkeypatch, tmp_path):
    true_cores = dials.util.mp.available_cores()
    job_ad = tmp_path / ".job.ad"
    monkeypatch.setenv("_CONDOR_JOB_AD", os.fspath(job_ad))

    # Test that we handle gracefully absence of file
    assert dials.util.mp.available_cores() == true_cores

    # Now write something sensible to the file
    job_ad.write_text(
        f"""\
CpusProvisioned = {true_cores + 1}
"""
    )
    assert dials.util.mp.available_cores() == true_cores + 1

    # Now check it is robust against parsing failures
    job_ad.write_text(
        """\
CpusProvisioned =
"""
    )
    assert dials.util.mp.available_cores() == true_cores
