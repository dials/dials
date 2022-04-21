from __future__ import annotations

from dials.util import parse_htcondor_job_classad


def test_parse_htcondor_job_classad(tmp_path):
    machine_ad = tmp_path / "condor_job_ad"
    machine_ad.write_text(
        """\
NumShadowStarts = 1
MemoryProvisioned = 4096
CpusProvisioned = 8
TransferIn = false
"""
    )
    class_ad = parse_htcondor_job_classad(machine_ad)
    assert class_ad.cpus_provisioned == 8
    assert class_ad.memory_provisioned == 4096
