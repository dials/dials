# LIBTBX_SET_DISPATCHER_NAME dev.dials.refstat_symmetry_analysis
from __future__ import annotations

import logging
import os
import time
from pathlib import Path

import iotbx
import libtbx.phil
from cctbx import crystal
from iotbx import reflection_file_reader, reflection_file_utils

import dials.util.log
from dials.algorithms.symmetry import refstat
from dials.util.options import ArgumentParser
from dials.util.version import dials_version

xr = refstat.registry()

logger = logging.getLogger("dials.command_line.refstat_symmetry_analysis")

phil_scope = libtbx.phil.parse(
    """
    sample_dir = None
        .type = path
        .help = "Path to the directory containing Olex2 sample data"
                "(e.g. Olex2/sample_data/)"

    check_dir = None
        .type = path
        .help = "Path to a parent directory from which statistics will be"
                "computed containing for all .res and .ins file pairs found"
                "in this and subdirectories."

    check_file = None
        .type = path
        .help="Path to an .ins or .res file to calculate statistics for a"
              "particular example"

    output {
        log = dials.refstat_symmetry_analysis.log
            .type = str
            .help = "The log filename"
    }
"""
)

help_message = """
Run refstat symmetry analysis on reflection data from SHELX files.

Examples::

  dev.dials.refstat_symmetry_analysis sample_dir=/path/to/Olex2-1.5/sample_data/

  dev.dials.refstat_symmetry_analysis check_dir=/path/to/check_dir/

  dev.dials.refstat_symmetry_analysis check_file=/path/to/dials.ins
"""


def basics():
    val = xr.describe()

    assert xr.elements["-21-"].is_shadowed_by([xr.elements["--n"]])
    assert not xr.elements["-21-"].is_shadowed_by([xr.elements["-n-"]])

    sgs = ["I 41/a m d", "P 1 21/c 1", "C 1 2/c 1", "P n a 21", "P 43 3 2"]
    extinctions = [xr.show_extinctions_for(sgn) for sgn in sgs]
    return val + "\n".join(extinctions)


def get_cs_hkl(file_base):
    ins_file = str(file_base) + ".ins"
    if not os.path.exists(ins_file):
        ins_file = str(file_base) + ".res"
    if not os.path.exists(ins_file):
        return None, None
    hkl_file = str(file_base) + ".hkl"
    if not os.path.exists(hkl_file):
        return None, None
    cs = iotbx.shelx.crystal_symmetry_from_ins.extract_from(file_name=ins_file)
    return cs, hkl_file


def check_reflections(file_base):
    cs, hkl_file = get_cs_hkl(str(file_base))
    assert cs is not None
    logger.info(
        "Original space group: %s"
        % (cs.space_group().match_tabulated_settings().hermann_mauguin())
    )
    check_reflections_(cs.unit_cell(), hkl_file)


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


def check_reflections_(cell, hkl_file):
    miller_array = get_miller_array(cell, hkl_file)
    logger.info("Read in %s reflections" % (len(miller_array.indices())))
    miller_array = miller_array.merge_equivalents(algorithm="gaussian").array()
    data = miller_array.data()
    sigmas = miller_array.sigmas()
    logger.info("Uniq in P1: %s" % (len(data)))
    timex = 10
    t = time.time()
    for r in range(timex):
        xr.process(miller_array.indices(), data, sigmas)
    logger.info("CPP processing time: %.3f" % (time.time() - t))
    t = time.time()
    try:
        for r in range(timex):
            xr.process_omp(miller_array.indices(), data, sigmas, -1)
        logger.info("CPP_omp processing time: %.3f" % (time.time() - t))
    except RuntimeError as e:
        if "Not implemented" in str(e):
            logger.info("CPP_omp processing not available.")
            pass

    xr.reset()

    sa = refstat.extinctions(miller_array)
    sa.analyse(scale_I_to=10000)
    logger.info(sa.show_stats())
    logger.info("Mean I(sig): %.3f(%.2f)/%s" % (sa.meanI, sa.mean_sig, sa.ref_count))
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
        logger.info(
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
        logger.info("%s: %s" % (sg.name, int(mp * 100)))

    logger.info(
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


def check_samples(samples_dir):
    samples_dir = Path(samples_dir)
    test_list = [
        samples_dir / "THPP" / "thpp",
        samples_dir / "ZP2" / "ZP2",
    ]
    for sample_base in test_list:
        try:
            logger.info("Testing: %s" % sample_base)
            check_reflections(sample_base)
        except Exception as e:
            import traceback

            logger.info(traceback.format_exc())
            logger.info("Failed to test %s: %s " % (sample_base, str(e)))


def check_dir(root_):
    def get_matches(cs, hkl_file, centering):
        miller_array = get_miller_array(cs.unit_cell(), hkl_file)
        miller_array = miller_array.merge_equivalents(algorithm="gaussian").array()
        xr.reset()
        sa = refstat.extinctions(miller_array)
        sa.analyse(scale_I_to=10000)
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
            logger.info(os.path.join(root, f))
            try:
                cs, hkl_path = get_cs_hkl(file_base)
                if cs is None:
                    continue
                original_sg_name = (
                    cs.space_group().match_tabulated_settings().hermann_mauguin()
                )
                if not original_sg_name:
                    continue
                logger.info("Original space group: %s" % (original_sg_name))
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
                logger.info("Matches: %s" % (", ".join(filtred_matches_names)))

            except Exception as e:
                import traceback

                logger.info(traceback.format_exc())
                logger.info("Failed on: %s, %s" % (file_full, str(e)))
    return stats


@dials.util.show_mail_on_error()
def run(args: list[str] = None, phil: libtbx.phil.scope = phil_scope) -> None:
    usage = "dev.dials.refstat_symmetry_analysis [options]"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=False,
        read_experiments=False,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging.
    dials.util.log.config(options.verbose, logfile=params.output.log)

    # Log the dials version
    logger.info(dials_version())

    # Log the difference between the PHIL scope definition and the active PHIL scope,
    # which will include the parsed user inputs.
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    if [params.sample_dir, params.check_dir, params.check_file].count(None) == 3:
        logger.info("No test paths provided. Only performing a basic test.")
        logger.info(basics())

    if params.check_file and os.path.exists(params.check_file):
        check_base = os.path.splitext(params.check_file)[0]
        logger.info("Testing: %s" % check_base)
        check_reflections(check_base)

    if params.sample_dir and os.path.exists(params.sample_dir):
        check_samples(params.sample_dir)

    if params.check_dir and os.path.exists(params.check_dir):
        stats = check_dir(params.check_dir)
        logger.info(stats)


if __name__ == "__main__":
    run()
