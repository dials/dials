from __future__ import absolute_import, division, print_function

import json
import math
import os
import pickle
import shutil

from dials.array_family import flex
import procrunner


def test2(dials_data, tmpdir):
    # Call dials.integrate
    result = procrunner.run(
        [
            "dials.integrate",
            dials_data("centroid_test_data").join("experiments.json"),
            "profile.fitting=False",
            "integration.integrator=3d",
            "prediction.padding=0",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    with tmpdir.join("integrated.refl").open("rb") as fh:
        table = pickle.load(fh)
    mask = table.get_flags(table.flags.integrated, all=False)
    assert len(table) == 1996
    assert mask.count(True) == 1666
    assert "id" in table
    for row in table:
        assert row["id"] == 0

    originaltable = table

    tmpdir.join("integrated.refl").remove()

    for i in range(1, 10):
        source = dials_data("centroid_test_data").join("centroid_000%d.cbf" % i)
        destination = source.new(
            dirname=tmpdir.strpath, basename="centroid_001%d.cbf" % i
        )
        source.copy(destination)

    with dials_data("centroid_test_data").join("experiments.json").open("r") as fh:
        j = json.load(fh)
    assert j["scan"][0]["image_range"] == [1, 9]
    j["scan"][0]["image_range"] = [11, 19]
    assert j["scan"][0]["oscillation"] == [0.0, 0.2]
    j["scan"][0]["oscillation"] = [360.0, 0.2]
    with tmpdir.join("models.expt").open("w") as fh:
        json.dump(j, fh)

    # Call dials.integrate
    result = procrunner.run(
        [
            "dials.integrate",
            "models.expt",
            "profile.fitting=False",
            "integration.integrator=3d",
            "prediction.padding=0",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    with tmpdir.join("integrated.refl").open("rb") as fh:
        table = pickle.load(fh)
    mask1 = table.get_flags(table.flags.integrated, all=False)
    assert len(table) == 1996
    assert mask1.count(True) == 1666
    mask2 = originaltable.get_flags(table.flags.integrated, all=False)
    assert mask1.all_eq(mask2)
    t1 = table.select(mask1)
    t2 = originaltable.select(mask1)
    Cal_P1 = t1["xyzcal.mm"].parts()[2]
    Cal_Z1 = t1["xyzcal.px"].parts()[2]
    Obs_Z1 = t1["xyzobs.px.value"].parts()[2]
    # Obs_P1 = t1['xyzobs.mm.value'].parts()[2]
    Cal_Z2 = t2["xyzcal.px"].parts()[2]
    Cal_P2 = t2["xyzcal.mm"].parts()[2]
    Obs_Z2 = t2["xyzobs.px.value"].parts()[2]
    # Obs_P2 = t2['xyzobs.mm.value'].parts()[2]
    diff_I = t1["intensity.sum.value"] - t2["intensity.sum.value"]
    diff_Cal_Z = Cal_Z1 - (Cal_Z2 + 10)
    diff_Obs_Z = Obs_Z1 - (Obs_Z2 + 10)
    diff_Cal_P = Cal_P1 - (Cal_P2 + 2 * math.pi)
    # diff_Obs_P = Obs_P1 - (Obs_P2 + 2*math.pi)
    assert flex.abs(diff_I).all_lt(1e-7)
    assert flex.abs(diff_Cal_Z).all_lt(1e-7)
    assert flex.abs(diff_Cal_P).all_lt(1e-7)
    assert flex.abs(diff_Obs_Z).all_lt(1e-7)
    # assert(flex.abs(diff_Obs_P).all_lt(1e-7))


def test_integration_with_sampling(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.integrate",
            dials_data("centroid_test_data").join("experiments.json"),
            "profile.fitting=False",
            "sampling.integrate_all_reflections=False",
            "prediction.padding=0",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    with tmpdir.join("integrated.refl").open("rb") as fh:
        table = pickle.load(fh)
    assert len(table) == 1000


def test_integration_with_sample_size(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.integrate",
            dials_data("centroid_test_data").join("experiments.json"),
            "profile.fitting=False",
            "sampling.integrate_all_reflections=False",
            "sampling.minimum_sample_size=500",
            "prediction.padding=0",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr

    with tmpdir.join("integrated.refl").open("rb") as fh:
        table = pickle.load(fh)
    assert len(table) == 500


def test_multi_sweep(dials_regression, run_in_tmpdir):
    result = procrunner.run(
        [
            "dials.integrate",
            os.path.join(
                dials_regression,
                "integration_test_data",
                "multi_sweep",
                "experiments.json",
            ),
            os.path.join(
                dials_regression,
                "integration_test_data",
                "multi_sweep",
                "indexed.pickle",
            ),
            "prediction.padding=0",
        ]
    )
    assert not result.returncode and not result.stderr
    assert os.path.exists("integrated.refl")

    with open("integrated.refl", "rb") as fh:
        table = pickle.load(fh)
    assert len(table) == 4020

    # Check the results
    T1 = table[:2010]
    T2 = table[2010:]
    ID1 = list(set(T1["id"]))
    ID2 = list(set(T2["id"]))
    assert len(ID1) == 1
    assert len(ID2) == 1
    assert ID1[0] == 0
    assert ID2[0] == 1
    I1 = T1["intensity.prf.value"]
    I2 = T2["intensity.prf.value"]
    F1 = T1.get_flags(T1.flags.integrated_prf)
    F2 = T2.get_flags(T2.flags.integrated_prf)
    assert F1 == F2
    I1 = I1.select(F1)
    I2 = I2.select(F2)
    assert flex.abs(I1 - I2) < 1e-6


def test_multi_lattice(dials_regression, tmpdir):
    result = procrunner.run(
        [
            "dials.integrate",
            os.path.join(
                dials_regression,
                "integration_test_data",
                "multi_lattice",
                "experiments.json",
            ),
            os.path.join(
                dials_regression,
                "integration_test_data",
                "multi_lattice",
                "indexed.pickle",
            ),
            "prediction.padding=0",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("integrated.refl").check()

    table = flex.reflection_table.from_file(tmpdir.join("integrated.refl"))
    assert len(table) == 5605

    # Check output contains from two lattices
    exp_id = list(set(table["id"]))
    assert len(exp_id) == 2

    # Check both lattices have integrated reflections
    mask = table.get_flags(table.flags.integrated_prf)
    table = table.select(mask)
    exp_id = list(set(table["id"]))
    assert len(exp_id) == 2


def test_output_rubbish(dials_data, tmpdir):
    result = procrunner.run(
        [
            "dials.index",
            dials_data("centroid_test_data").join("datablock.json"),
            dials_data("centroid_test_data").join("strong.pickle"),
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("indexed.expt").check(file=1)
    assert tmpdir.join("indexed.refl").check(file=1)

    # Call dials.integrate
    result = procrunner.run(
        [
            "dials.integrate",
            "indexed.expt",
            "indexed.refl",
            "profile.fitting=False",
            "prediction.padding=0",
        ],
        working_directory=tmpdir,
    )
    assert not result.returncode and not result.stderr
    assert tmpdir.join("integrated.refl").check(file=1)

    with tmpdir.join("integrated.refl").open("rb") as fh:
        table = pickle.load(fh)

    assert table.get_flags(table.flags.bad_reference) > 0
    assert "id" in table
    for row in table:
        assert row["id"] == 0


def test_integrate_with_kapton(dials_regression, tmpdir):
    tmpdir.chdir()
    loc = tmpdir.strpath

    pickle_name = "idx-20161021225550223_indexed.pickle"
    json_name = "idx-20161021225550223_refined_experiments.json"
    image_name = "20161021225550223.pickle"
    pickle_path = os.path.join(
        dials_regression, "integration_test_data", "stills_PSII", pickle_name
    )
    json_path = os.path.join(
        dials_regression, "integration_test_data", "stills_PSII", json_name
    )
    image_path = os.path.join(
        dials_regression, "integration_test_data", "stills_PSII", image_name
    )

    assert os.path.exists(pickle_path)
    assert os.path.exists(json_path)
    shutil.copy(pickle_path, loc)
    shutil.copy(image_path, loc)

    with open(json_name, "wb") as w, open(json_path, "rb") as r:
        w.write(r.read() % loc.replace("\\", "\\\\"))

    templ_phil = """
      output {
        experiments = 'idx-20161021225550223_integrated_experiments_%s.expt'
        reflections = 'idx-20161021225550223_integrated_%s.refl'
      }
      integration {
        lookup.mask = '%s'
        integrator = stills
        profile.fitting = False
        background.algorithm = simple
        debug {
          output = True
          separate_files = False
          split_experiments = False
        }
      }
      profile {
        gaussian_rs.min_spots.overall = 0
      }
      absorption_correction {
        apply = %s
        algorithm = fuller_kapton
        fuller_kapton {
          smart_sigmas = True
        }
      }
"""
    without_kapton_phil = templ_phil % (
        "nokapton",
        "nokapton",
        os.path.join(
            dials_regression, "integration_test_data", "stills_PSII", "mask.pickle"
        ).replace("\\", "\\\\"),
        "False",
    )
    with_kapton_phil = templ_phil % (
        "kapton",
        "kapton",
        os.path.join(
            dials_regression, "integration_test_data", "stills_PSII", "mask.pickle"
        ).replace("\\", "\\\\"),
        "True",
    )

    with open("integrate_without_kapton.phil", "wb") as f:
        f.write(without_kapton_phil)

    with open("integrate_with_kapton.phil", "wb") as f:
        f.write(with_kapton_phil)

    # Call dials.integrate with and without kapton correction
    for phil in "integrate_without_kapton.phil", "integrate_with_kapton.phil":
        result = procrunner.run(["dials.integrate", pickle_name, json_name, phil])
        assert not result.returncode and not result.stderr

    results = []
    for mode in "kapton", "nokapton":
        result = os.path.join(loc, "idx-20161021225550223_integrated_%s.refl" % mode)
        with open(result, "rb") as f:
            table = pickle.load(f)
        millers = table["miller_index"]
        test_indices = {"zero": (-5, 2, -6), "low": (-2, -20, 7), "high": (-1, -10, 4)}
        test_rows = {k: millers.first_index(v) for k, v in test_indices.items()}
        test_I_sigsqI = {
            k: (table[v]["intensity.sum.value"], table[v]["intensity.sum.variance"])
            for k, v in test_rows.items()
        }
        results.append(test_I_sigsqI)
    assert results[0]["zero"][0] == results[1]["zero"][0]
    assert results[0]["zero"][1] - results[1]["zero"][1] < 0.0001
    assert False not in [results[0]["low"][i] > results[1]["low"][i] for i in (0, 1)]
    assert False not in [results[0]["high"][i] > results[1]["high"][i] for i in (0, 1)]
