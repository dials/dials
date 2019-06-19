"""
Test for dials.assign_experiment_identifiers
"""
from __future__ import absolute_import, division, print_function

import os

import procrunner
from dials.array_family import flex
from dxtbx.serialize import load


def run_assign_identifiers(pickle_path_list, sweep_path_list, extra_args):
    command = (
        ["dials.assign_experiment_identifiers"]
        + pickle_path_list
        + sweep_path_list
        + extra_args
    )
    print(command)
    procrunner.run(command).check_returncode()
    assert os.path.exists("assigned.expt")
    assert os.path.exists("assigned.refl")


def test_assign_identifiers(dials_regression, run_in_tmpdir):
    """Test for dials.assign_experiment_identifiers"""
    pickle_path_list = []
    sweep_path_list = []
    data_dir = os.path.join(dials_regression, "xia2-28")
    for i in [20, 25]:
        pickle_path_list.append(os.path.join(data_dir, str(i) + "_integrated.pickle"))
        sweep_path_list.append(
            os.path.join(data_dir, str(i) + "_integrated_experiments.json")
        )

    run_assign_identifiers(pickle_path_list, sweep_path_list, extra_args=[])

    r = flex.reflection_table.from_pickle("assigned.refl")
    e = load.experiment_list("assigned.expt", check_format=False)
    r.assert_experiment_identifiers_are_consistent(e)
    assert list(r.experiment_identifiers().values()) == ["0", "1"]
    assert list(r.experiment_identifiers().keys()) == [0, 1]
    assert list(e.identifiers()) == ["0", "1"]

    # now run again, with already assigned data
    pickle_path_list = ["assigned.refl"]
    sweep_path_list = ["assigned.expt"]
    run_assign_identifiers(pickle_path_list, sweep_path_list, extra_args=[])

    r = flex.reflection_table.from_pickle("assigned.refl")
    e = load.experiment_list("assigned.expt", check_format=False)
    r.assert_experiment_identifiers_are_consistent(e)
    assert list(r.experiment_identifiers().values()) == ["0", "1"]
    assert list(r.experiment_identifiers().keys()) == [0, 1]
    assert list(e.identifiers()) == ["0", "1"]

    # now run again, with adding more data
    pickle_path_list = ["assigned.refl"]
    sweep_path_list = ["assigned.expt"]
    for i in [30, 35]:
        pickle_path_list.append(os.path.join(data_dir, str(i) + "_integrated.pickle"))
        sweep_path_list.append(
            os.path.join(data_dir, str(i) + "_integrated_experiments.json")
        )

    run_assign_identifiers(
        pickle_path_list, sweep_path_list, extra_args=["identifiers=0 5 10 15"]
    )

    r = flex.reflection_table.from_pickle("assigned.refl")
    e = load.experiment_list("assigned.expt", check_format=False)
    r.assert_experiment_identifiers_are_consistent(e)
    assert list(r.experiment_identifiers().values()) == ["0", "5", "10", "15"]
    assert list(r.experiment_identifiers().keys()) == [0, 1, 2, 3]
    assert list(e.identifiers()) == ["0", "5", "10", "15"]
