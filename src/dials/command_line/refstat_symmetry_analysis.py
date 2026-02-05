# LIBTBX_SET_DISPATCHER_NAME dev.dials.refstat_symmetry_analysis
from __future__ import annotations

import logging
import os
import sys
import time
from pathlib import Path

import cctbx.miller
import iotbx
import libtbx.phil
from cctbx import crystal
from iotbx import reflection_file_reader, reflection_file_utils

import dials.util.log
from dials.algorithms.symmetry import (
    get_subset_for_symmetry,
    refstat,
    resolution_filter_from_reflections_experiments,
)
from dials.array_family import flex
from dials.command_line.symmetry import phil_scope as symmetry_phil_scope
from dials.util.filter_reflections import filter_reflection_table
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
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

    sigma_level = 5.0
        .type = float
        .help = "Sigma level to use to identify systematic absences"

    output {
        log = dev.dials.refstat_symmetry_analysis.log
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
    return val + "\n" + "\n".join(extinctions)


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


def load_miller_array_and_centering_from_hkl(file_base):
    cs, hkl_file = get_cs_hkl(str(file_base))
    assert cs is not None
    original_sg_name = cs.space_group().match_tabulated_settings().hermann_mauguin()
    centering = original_sg_name[0]
    logger.info(
        "Original space group: %s"
        % (cs.space_group().match_tabulated_settings().hermann_mauguin())
    )
    miller_array = get_miller_array(cs.unit_cell(), hkl_file)
    logger.info("Read in %s reflections" % (len(miller_array.indices())))
    return miller_array, centering


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


def format_sg_name(name):
    toks = name.split()
    if len(toks) == 4:
        if toks[1] == "1" and toks[3] == "1":
            return "%s%s" % (toks[0], toks[2])
    return "".join(toks)


def check_reflections(miller_array, centering="P", sigma_level=5.0):
    miller_array = miller_array.merge_equivalents(algorithm="gaussian").array()
    data = miller_array.data()
    sigmas = miller_array.sigmas()
    logger.info("Uniq in P1: %s" % (len(data)))
    timex = 1
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

    sa = refstat.extinctions(miller_array, sigma_level=sigma_level)
    sa.analyse(scale_I_to=1)
    logger.info(sa.show_stats())
    logger.info("Mean I(sig): %.3f(%.2f)/%s" % (sa.meanI, sa.mean_sig, sa.ref_count))
    matches = sa.get_all_matching_space_groups(centering=centering)
    filtered_matches = sa.get_filtered_matching_space_groups(matches=matches)
    # merge_test object
    t = refstat.merge_test(miller_array.indices(), data, sigmas)
    for sg, mp in matches:
        # limit to the filteres selection, maybe for non-verbose only?
        if sg not in filtered_matches:
            continue
        weak_stats = t.sysabs_test(sg, sa.scale)
        wI = weak_stats.weak_I_sum / weak_stats.weak_count
        wIs = (weak_stats.weak_sig_sq_sum / weak_stats.weak_count) ** 0.5
        if wI > 5 * wIs:
            continue
        merge_stats = t.merge_test(sg)
        sI = weak_stats.strong_I_sum / weak_stats.strong_count
        sIs = (weak_stats.strong_sig_sq_sum / weak_stats.strong_count) ** 0.5
        logger.info(
            "Inconsistent equivalents: %s, r_int: %.3f, weak: %.3f(%.2f)/%s %.3f, strong: %.3f(%.2f)/%s %.3f"
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
        logger.info("%s: %s%% matches" % (format_sg_name(sg.name), int(mp * 100)))

    logger.info(
        "All matches: %s"
        % (", ".join([format_sg_name(sg.name) for sg in filtered_matches]))
    )
    centric, acentric = [], []
    max_top = 3
    for sg in filtered_matches:
        if sg.is_centric():
            if len(centric) < max_top:
                centric.append(sg)
        else:
            if len(acentric) < max_top:
                acentric.append(sg)
    logger.info(
        "Top acentric matches: %s"
        % (", ".join([format_sg_name(sg.name) for sg in acentric]))
    )
    logger.info(
        "Top centric matches: %s"
        % (", ".join([format_sg_name(sg.name) for sg in centric]))
    )


def check_samples(samples_dir, sigma_level=5.0):
    samples_dir = Path(samples_dir)
    test_list = [
        samples_dir / "THPP" / "thpp",
        samples_dir / "ZP2" / "ZP2",
    ]
    for sample_base in test_list:
        try:
            logger.info("Testing: %s" % sample_base)
            ma, centering = load_miller_array_and_centering_from_hkl(sample_base)
            check_reflections(ma, centering, sigma_level=sigma_level)
        except Exception as e:
            import traceback

            logger.info(traceback.format_exc())
            logger.info("Failed to test %s: %s " % (sample_base, str(e)))


def check_dir(root_, sigma_level=5.0):
    def get_matches(cs, hkl_file, centering):
        miller_array = get_miller_array(cs.unit_cell(), hkl_file)
        miller_array = miller_array.merge_equivalents(algorithm="gaussian").array()
        xr.reset()
        sa = refstat.extinctions(miller_array, sigma_level=sigma_level)
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


def check_dials_input(experiments, reflections):
    # Perform the same steps that dials.symmetry does to prepare the data
    params = symmetry_phil_scope.extract()
    refls_for_sym = get_subset_for_symmetry(
        experiments, reflections, params.exclude_images
    )
    d_min = resolution_filter_from_reflections_experiments(
        refls_for_sym,
        experiments,
        params.min_i_mean_over_sigma_mean,
        params.min_cc_half,
    )
    d_max = None
    if len(reflections) > 1:
        reflection_table = flex.reflection_table()
        for table in refls_for_sym:
            reflection_table.extend(table)
    else:
        reflection_table = refls_for_sym[0]

    # Filter reflections and make an intensity choice
    partiality_threshold = 0.99
    if (
        "inverse_scale_factor" in reflection_table
        and "intensity.scale.value" in reflection_table
    ):
        logger.info("Performing systematic absence checks on scaled data")
        reflections = filter_reflection_table(
            reflection_table,
            intensity_choice=["scale"],
            d_min=d_min,
            d_max=d_max,
            partiality_threshold=partiality_threshold,
        )
        reflections["intensity"] = reflections["intensity.scale.value"]
        reflections["variance"] = reflections["intensity.scale.variance"]
    elif "intensity.prf.value" in reflection_table:
        logger.info(
            "Performing systematic absence checks on unscaled profile-integrated data"
        )
        reflections = filter_reflection_table(
            reflection_table,
            intensity_choice=["profile"],
            d_min=d_min,
            d_max=d_max,
            partiality_threshold=partiality_threshold,
        )
        reflections["intensity"] = reflections["intensity.prf.value"]
        reflections["variance"] = reflections["intensity.prf.variance"]
    else:
        logger.info(
            "Performing systematic absence checks on unscaled summation-integrated data"
        )
        reflections = filter_reflection_table(
            reflection_table,
            intensity_choice=["sum"],
            d_min=d_min,
            d_max=d_max,
            partiality_threshold=partiality_threshold,
        )
        reflections["intensity"] = reflections["intensity.sum.value"]
        reflections["variance"] = reflections["intensity.sum.variance"]

    cs = crystal.symmetry(experiments[0].crystal.get_unit_cell(), "P1")
    centering = (
        experiments[0]
        .crystal.get_space_group()
        .match_tabulated_settings()
        .hermann_mauguin()[0]
    )
    miller_set = cctbx.miller.set(
        crystal_symmetry=cs,
        indices=reflections["miller_index"],
        anomalous_flag=False,
    )
    i_obs = cctbx.miller.array(miller_set, data=reflections["intensity"])
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas(cctbx.array_family.flex.sqrt(reflections["variance"]))
    i_obs.set_info(
        cctbx.miller.array_info(source="DIALS", source_type="reflection_tables")
    )

    check_reflections(i_obs, centering)


@dials.util.show_mail_on_error()
def run(args: list[str] = None, phil: libtbx.phil.scope = phil_scope) -> None:
    usage = "dev.dials.refstat_symmetry_analysis [options]"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
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

    # Run analysis in the various modes
    if params.input.experiments and params.input.reflections:
        reflections, experiments = reflections_and_experiments_from_files(
            params.input.reflections, params.input.experiments
        )
        reflections = parse_multiple_datasets(reflections)
        try:
            experiments, reflections = assign_unique_identifiers(
                experiments, reflections
            )
            check_dials_input(experiments, reflections)
        except ValueError:
            pass

    elif params.check_file and os.path.exists(params.check_file):
        check_base = os.path.splitext(params.check_file)[0]
        logger.info("Testing: %s" % check_base)
        ma, centering = load_miller_array_and_centering_from_hkl(check_base)
        check_reflections(ma, centering, sigma_level=params.sigma_level)
        sys.exit(0)

    elif params.sample_dir and os.path.exists(params.sample_dir):
        check_samples(params.sample_dir, sigma_level=params.sigma_level)
        sys.exit(0)

    elif params.check_dir and os.path.exists(params.check_dir):
        stats = check_dir(params.check_dir, sigma_level=params.sigma_level)
        logger.info(stats)
        sys.exit(0)

    else:
        parser.print_help()
        logger.info("No valid input files. Only performing a basic test.")
        logger.info(basics())


if __name__ == "__main__":
    run()
