"""Test the dials.space_group command line script."""
import os
import procrunner
from dxtbx.serialize import load


def test_systematic_absences_script(dials_data, run_in_tmpdir):
    """Test the command line script with real data. Proteinase K in P41"""
    location = dials_data("vmxi_proteinase_k_sweeps")

    command = ["dials.space_group", "output.experiments=symmetrized.expt"]
    command.append(location.join("experiments_0.json").strpath)
    command.append(location.join("reflections_0.pickle").strpath)

    result = procrunner.run(command)
    assert not result.returncode and not result.stderr
    assert os.path.exists("absences.html")
    assert os.path.exists("symmetrized.expt")
    exps = load.experiment_list("symmetrized.expt", check_format=False)
    assert str(exps[0].crystal.get_space_group().info()) == "P 41"

    # Now try with a d_min
    command += ["d_min=4.0"]
    result = procrunner.run(command)
    assert not result.returncode and not result.stderr
    exps = load.experiment_list("symmetrized.expt", check_format=False)
    assert str(exps[0].crystal.get_space_group().info()) == "P 41"
