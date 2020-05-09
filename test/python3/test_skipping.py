import six


def test_that_this_test_runs_only_on_python3(python3):
    assert six.PY3
