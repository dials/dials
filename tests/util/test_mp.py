import dials.util.mp


def test_core_detection():
    # We can't make any assumptions about the testing environment,
    # but we know there will be at least one available core, and
    # the function must return a positive integer in any case.
    assert dials.util.mp.available_cores() >= 1


def test_core_detection_override(monkeypatch):
    # For now at least the NSLOTS environment variable should
    # override whatever available_cores() returns otherwise.
    true_cores = dials.util.mp.available_cores()
    monkeypatch.setenv("NSLOTS", str(true_cores + 1))
    assert dials.util.mp.available_cores() == true_cores + 1
