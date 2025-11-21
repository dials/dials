from __future__ import annotations

import os
import time

import iotbx
from cctbx import crystal
from iotbx import reflection_file_reader, reflection_file_utils

# import cctbx.sgtbx.refstat as refstat
from dials.algorithms.symmetry import refstat

# this is needed only if reloading module from within Olex2
# reload(refstat)

xr = refstat.registry()


def basics():
    xr.describe()

    assert xr.elements["-21-"].is_shadowed_by([xr.elements["--n"]])
    assert not xr.elements["-21-"].is_shadowed_by([xr.elements["-n-"]])

    sgs = ["I 41/a m d", "P 1 21/c 1", "C 1 2/c 1", "P n a 21", "P 43 3 2"]
    for sgn in sgs:
        xr.show_extinctions_for(sgn)


def get_cs_hkl(file_base):
    ins_file = file_base + ".ins"
    if not os.path.exists(ins_file):
        ins_file = file_base + ".res"
    if not os.path.exists(ins_file):
        return None, None
    hkl_file = file_base + ".hkl"
    if not os.path.exists(hkl_file):
        return None, None
    cs = iotbx.shelx.crystal_symmetry_from_ins.extract_from(file_name=ins_file)
    return cs, hkl_file


def test_reflections(file_base):
    cs, hkl_file = get_cs_hkl(file_base)
    assert cs != None
    print(
        "Original space group: %s"
        % (cs.space_group().match_tabulated_settings().hermann_mauguin())
    )
    test_reflections_(cs.unit_cell(), hkl_file)


def get_miller_array(cell, hkl_file):
    cs = crystal.symmetry(cell, "P1")
    reflections_server = reflection_file_utils.reflection_file_server(
        crystal_symmetry=cs,
        reflection_files=[
            reflection_file_reader.any_reflection_file(
                "hklf4=%s" % (hkl_file), strict=False
            )
        ],
    )
    return reflections_server.get_miller_arrays(None)[0]


def test_reflections_(cell, hkl_file):
    miller_array = get_miller_array(cell, hkl_file)
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

    xr.reset()

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


def test_olx():
    """This can be run from within Olex2 on currently loaded structure"""
    try:
        import olx

        cell = [float(x) for x in olx.xf.au.GetCell().split(",")]
        test_reflections(cell, olx.HKLSrc())
    except Exception as e:
        print("Failed to run olx test: %s" % str(e))


def test_samples(samples_dir):
    test_list = [
        "C:/Program Files/Olex2-1.5-alpha/sample_data/THPP/thpp",
        "C:/Program Files/Olex2-1.5-alpha/sample_data/ZP2/ZP2",
    ]
    for file_base in test_list:
        sample_base = os.path.join(samples_dir, file_base)
        try:
            print("Testing: %s" % sample_base)
            test_reflections(sample_base)
        except Exception as e:
            import traceback

            print(traceback.format_exc())
            print("Failed to test %s: %s " % (sample_base, str(e)))


def test_dir(root_):
    def get_matches(cs, hkl_file, centering):
        miller_array = get_miller_array(cs.unit_cell(), hkl_file)
        # print("Read in %s reflections" % (len(miller_array.indices())))
        miller_array = miller_array.merge_equivalents(algorithm="gaussian").array()
        xr.reset()
        sa = refstat.extinctions(miller_array)
        sa.analyse(scale_I_to=10000)
        # sa.print_stats()
        # print("Mean I(sig): %.3f(%.2f)/%s" % (sa.meanI, sa.mean_sig, sa.ref_count))
        matches = sa.get_all_matching_space_groups(centering=centering)
        if not matches:
            return (None, None)
        return matches, sa.get_filtered_matching_space_groups(matches=matches)

    stats = {
        "0": 0,
        "P1": 0,
        "+": 0,
        "100": 0,
    }
    for root, dirs, files in os.walk(root_):
        for f in files:
            name, ext = os.path.splitext(f)
            if ext.lower() not in [".res", ".ins"]:
                continue
            file_full = os.path.join(root, f)
            file_base = os.path.splitext(file_full)[0]
            print(os.path.join(root, f))
            try:
                cs, hkl_path = get_cs_hkl(file_base)
                if cs is None:
                    continue
                original_sg_name = (
                    cs.space_group().match_tabulated_settings().hermann_mauguin()
                )
                if not original_sg_name:
                    continue
                print("Original space group: %s" % (original_sg_name))
                matches, filtred_matches = get_matches(
                    cs, hkl_path, original_sg_name[0]
                )
                if matches is None:
                    if original_sg_name not in ("P 1", "P -1"):
                        stats["0"] += 1
                    else:
                        stats["P1"] += 1
                    continue

                filtred_matches_names = [sg.name for sg in filtred_matches]
                if original_sg_name in filtred_matches_names:
                    stats["100"] += 1
                else:
                    stats["+"] += 1
                print("Matches: %s" % (", ".join(filtred_matches_names)))

            except Exception as e:
                import traceback

                print(traceback.format_exc())
                print("Failed on: %s, %s" % (file_full, str(e)))
    return stats


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser(usage="unleash_olex2.py [options]")
    parser.add_option(
        "--sample_dir",
        dest="sample_dir",
        default="C:/Program Files/Olex2-1.5-alpha/sample_data/",
        help="Path to the directory contains date for test_samles",
    )
    parser.add_option(
        "--test_dir",
        dest="test_dir",
        default="d:/devel/data/",
        help="Path to the directory contains date for test_dir."
        " Runs stats on all structure in the directory.",
    )
    parser.add_option(
        "--test_file",
        dest="test_file",
        default="",
        help="Path to a file to test for chasing particular examples",
    )
    options = parser.parse_args()[0]

    if options.test_file:
        if os.path.exists(options.test_file):
            test_base = os.path.splitext(options.test_file)[0]
            print("Testing: %s" % test_base)
            test_reflections(test_base)
        else:
            print("Specified test_file doe snot exist: " % options.test_file)
        exit(0)

    basics()

    if os.path.exists(options.sample_dir):
        test_samples(options.sample_dir)
    else:
        print("Skipping test_samples -  dir does not exist")

    if os.path.exists(options.test_dir):
        stats = test_dir(options.test_dir)
        print(stats)
    else:
        print("Skipping test_dir -  dir does not exist")
