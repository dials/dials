from __future__ import annotations

import os
import time
from importlib import reload

from cctbx import crystal
from iotbx import reflection_file_reader, reflection_file_utils
from smtbx.regression.test_data import fnames

from dials.algorithms.symmetry import refstat

reload(refstat)

xr = refstat.registry()
xr.describe()

assert xr.elements["-21-"].is_shadowed_by([xr.elements["--n"]])
assert not xr.elements["-21-"].is_shadowed_by([xr.elements["-n-"]])

sgs = ["I 41/a m d", "P 1 21/c 1", "C 1 2/c 1", "P n a 21", "P 43 3 2"]
for sgn in sgs:
    xr.show_extinctions_for(sgn)


def test_reflections():
    # Currently limited to just THPP data from Olex2 sample data. Other datasets
    # can be added as a filebase fixture
    file_base = fnames.thpp_ins[:-4]

    ins_file = file_base + ".ins"
    if not os.path.exists(ins_file):
        ins_file = file_base + ".res"
    assert os.path.exists(ins_file)
    hkl_file = file_base + ".hkl"
    assert os.path.exists(hkl_file)
    # this seems not to work??
    # xs = cctbx.xray.structure.from_shelx(ins_file, strictly_shelxl=False)
    # test_reflections_(xs.unit_cell(), hkl_file)
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
    timex = 10
    t = time.time()
    for r in range(timex):
        xr.process(miller_array.indices(), data, sigmas)
    print("CPP processing time: %.3f" % (time.time() - t))
    t = time.time()
    for r in range(timex):
        xr.process_omp(miller_array.indices(), data, sigmas, -1)
    print("CPP_omp processing time: %.3f" % (time.time() - t))

    for x in xr:
        x.reset()

    sa = refstat.extinctions(miller_array)
    sa.analyse(scale_I_to=10000)
    sa.print_stats()
    print("Mean I(sig): %.3f(%.2f)/%s" % (sa.meanI, sa.mean_sig, sa.ref_count))
    matches = sa.get_all_matching_space_groups()

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
        print(
            "Inonsistent eq: %s, r_int: %.3f, w: %.3f(%.2f)/%s %.3f, s: %.3f(%.2f)/%s %.3f"
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
        print("%s: %s" % (sg.name, int(mp * 100)))

    print(
        "Matches: %s"
        % (
            ", ".join(
                [
                    sg.name
                    for sg in sa.get_filtered_matching_space_groups(matches=matches)
                ]
            )
        )
    )
