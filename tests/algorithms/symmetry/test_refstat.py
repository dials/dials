from __future__ import annotations

import os
import shutil
import subprocess
import time

import pytest

from cctbx import crystal
from iotbx import reflection_file_reader, reflection_file_utils
from smtbx.regression.test_data import fnames

from dials.algorithms.symmetry import refstat


def test_registry():
    xr = refstat.registry()
    description = xr.describe()
    assert (
        description
        == """b--
	-x,y+1/2,z
c--
	-x,y,z+1/2
n--
	-x,y+1/2,z+1/2
d--
	-x,y+1/4,z+1/4
-a-
	x+1/2,-y,z
-c-
	x,-y,z+1/2
-n-
	x+1/2,-y,z+1/2
-d-
	x+1/4,-y,z+1/4
--a
	x+1/2,y,-z
--b
	x,y+1/2,-z
--n
	x+1/2,y+1/2,-z
--d
	x+1/4,y+1/4,-z
21--
	x+1/2,-y,-z
  Shadowed by: -a- -n- --a --n
-21-
	-x,y+1/2,-z
  Shadowed by: b-- n-- --b --n
--21
	-x,-y,z+1/2
  Shadowed by: c-- n-- -c- -n-
31
	-y,x-y,z+1/3
	-x+y,-x,z+2/3
41
	-y,x,z+1/4
	-x,-y,z+1/2
	y,-x,z+3/4
  Shadowed by: c-- n-- -c- -n- --21 d--
42
	-y,x,z+1/2
	-x,-y,z
	y,-x,z+1/2
  Shadowed by: c-- n-- -c- -n- --21
61
	x-y,x,z+1/6
	-y,x-y,z+1/3
	-x,-y,z+1/2
	-x+y,-x,z+2/3
	y,-x+y,z+5/6"""
    )

    assert xr.elements["-21-"].is_shadowed_by([xr.elements["--n"]])
    assert not xr.elements["-21-"].is_shadowed_by([xr.elements["-n-"]])

    sgs = ["I 41/a m d", "P 1 21/c 1", "C 1 2/c 1", "P n a 21", "P 43 3 2"]
    expected_strs = [
        "Extinction elements for I 41/a m d: n-- -n- 21-- -21- --21 41 --n",
        "Extinction elements for P 1 21/c 1: -c- -21-",
        "Extinction elements for C 1 2/c 1: -n- -21-",
        "Extinction elements for P n a 21: n-- -a- --21",
        "Extinction elements for P 43 3 2: 21-- -21- --21 41",
    ]
    for sgn, expected in zip(sgs, expected_strs):
        assert xr.show_extinctions_for(sgn) == expected


def test_reflections():
    # Currently limited to just THPP data from Olex2 sample data.
    file_base = fnames.thpp_ins[:-4]

    ins_file = file_base + ".ins"
    if not os.path.exists(ins_file):
        ins_file = file_base + ".res"
    assert os.path.exists(ins_file)
    hkl_file = file_base + ".hkl"
    assert os.path.exists(hkl_file)
    cell = None
    for l in open(ins_file).readlines():
        if l.startswith("CELL"):
            cell = [float(x) for x in l.split()[2:]]
            break
    assert cell is not None and len(cell) == 6

    cs = crystal.symmetry(cell, "P1")
    reflections_server = reflection_file_utils.reflection_file_server(
        crystal_symmetry=cs,
        reflection_files=[
            reflection_file_reader.any_reflection_file(
                "hklf4=%s" % (hkl_file), strict=False
            )
        ],
    )
    miller_array = reflections_server.get_miller_arrays(None)[0]
    print("Read in %s reflections" % (len(miller_array.indices())))
    miller_array = miller_array.merge_equivalents(algorithm="gaussian").array()
    data = miller_array.data()
    sigmas = miller_array.sigmas()
    print("Uniq in P1: %s" % (len(data)))
    xr = refstat.registry()
    timex = 10
    t = time.time()
    for r in range(timex):
        xr.process(miller_array.indices(), data, sigmas)
    print("CPP processing time: %.3f" % (time.time() - t))
    t = time.time()
    try:
        for r in range(timex):
            xr.process_omp(miller_array.indices(), data, sigmas, -1)
        print("CPP_omp processing time: %.3f" % (time.time() - t))
    except RuntimeError as e:
        if "Not implemented" in str(e):
            print("CPP_omp processing not available.")
            pass

    for x in xr:
        x.reset()

    sa = refstat.extinctions(miller_array)
    sa.analyse(scale_I_to=10000)
    stats = sa.show_stats()
    assert (
        stats
        == """b--  (  439):            94.23(  0.81) -
c--  (  446):           111.75(  0.82) -
n--  (  437):            64.44(  0.34) -
d--  (  658):            87.07(  0.53) -
-a-  (  212):            48.70(  0.41) -
-c-  (  212):            48.64(  0.41) -
-n-  (  204):            -0.03(  0.03) +
-d-  (  307):            91.55(  2.41) -
--a  (  307):            44.25(  0.19) -
--b  (  302):           129.15(  1.65) -
--n  (  305):           145.22(  1.64) -
--d  (  457):           165.07(  1.97) -
21-- (   10):             0.99(  0.13) +-
-21- (   20):             0.05(  0.09) +
--21 (   14):             0.08(  0.12) +-
31   (   18):            60.29(  1.52) -
41   (   20):            11.39(  0.25) -
42   (   14):             0.08(  0.12) +-
61   (   22):            49.35(  1.24) -"""
    )

    assert pytest.approx(sa.meanI, 0.001) == 53.880
    assert pytest.approx(sa.mean_sig, 0.001) == 10.450
    assert sa.ref_count == 11372

    expected_output = """Inconsistent eq: 30, r_int: 4.264, w: -0.025(0.47)/224 -0.054, s: 190.118(49.62)/231 3.832
P 1 21/n 1: 40
Inconsistent eq: 775, r_int: 48.289, w: 96.537(19.86)/975 4.861, s: 142.917(41.86)/1011 3.414
P c n n: 40
Inconsistent eq: 775, r_int: 48.289, w: 81.558(21.44)/543 3.804, s: 168.369(46.87)/544 3.593
P m n n: 60
Inconsistent eq: 775, r_int: 48.289, w: 55.842(11.51)/984 4.852, s: 183.298(44.96)/1002 4.077
P b n a: 60
Inconsistent eq: 775, r_int: 48.289, w: 90.476(19.73)/982 4.586, s: 149.169(42.03)/1004 3.549
P c n b: 60
Inconsistent eq: 775, r_int: 48.289, w: 70.926(21.21)/550 3.345, s: 139.000(36.46)/1436 3.813
P m n b: 80
Inconsistent eq: 775, r_int: 48.289, w: 60.224(13.59)/687 4.431, s: 151.839(39.52)/1299 3.842
P b n m: 80"""
    matches = sa.get_all_matching_space_groups()
    lines = []
    for sg, mp in matches:
        t = refstat.merge_test(miller_array.indices(), data, sigmas)
        weak_stats = t.sysabs_test(sg, sa.scale)
        wI = weak_stats.weak_I_sum / weak_stats.weak_count
        wIs = (weak_stats.weak_sig_sq_sum / weak_stats.weak_count) ** 0.5
        if wI > 5 * wIs:
            continue
        merge_stats = t.merge_test(sg)
        sI = weak_stats.strong_I_sum / weak_stats.strong_count
        sIs = (weak_stats.strong_sig_sq_sum / weak_stats.strong_count) ** 0.5
        lines.append(
            "Inconsistent eq: %s, r_int: %.3f, w: %.3f(%.2f)/%s %.3f, s: %.3f(%.2f)/%s %.3f"
            % (
                merge_stats.inconsistent_count,
                merge_stats.r_int * 100,
                wI,
                wIs,
                weak_stats.weak_count,
                wI / wIs,
                sI,
                sIs,
                weak_stats.strong_count,
                sI / sIs,
            )
        )
        lines.append("%s: %s" % (sg.name, int(mp * 100)))
    assert "\n".join(lines) == expected_output
    filtered_matches = sa.get_filtered_matching_space_groups(matches=matches)
    assert len(filtered_matches) == 1
    assert filtered_matches[0].name == "P 1 21/n 1"


def test_refstat_symmetry_analysis_check_file(tmp_path):
    result = subprocess.run(
        (
            shutil.which("dev.dials.refstat_symmetry_analysis"),
            f"check_file={fnames.thpp_ins}",
        ),
        cwd=tmp_path,
        capture_output=True,
        text=True,  # Convert bytes to strings and normalizes line endings
    )
    assert not result.check_returncode()
    assert result.stdout.endswith(
        "P21/n: 40% matches\nAll matches: P21/n\nTop acentric matches: \nTop centric matches: P21/n\n"
    )


def test_refstat_symmetry_analysis_dials_input(dials_data, tmp_path):
    quartz = dials_data("quartz_processed")
    result = subprocess.run(
        (
            shutil.which("dev.dials.refstat_symmetry_analysis"),
            quartz / "integrated.expt",
            quartz / "integrated.refl",
        ),
        cwd=tmp_path,
        capture_output=True,
        text=True,  # Convert bytes to strings and normalizes line endings
    )
    assert not result.check_returncode()
    assert result.stdout.endswith(
        "P3221: 100% matches\nAll matches: P31, P32, P3121, P3221\nTop acentric matches: P31, P32, P3121\nTop centric matches: \n"
    )
