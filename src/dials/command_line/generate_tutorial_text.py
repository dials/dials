# LIBTBX_SET_DISPATCHER_NAME dev.dials.generate_tutorial_text

from __future__ import annotations

import argparse
import functools
import json
import os
import pathlib
import shutil
import sys
import tempfile
import time

import procrunner

import dials.util
import dials_data.download


def _command_runner(
    command, output_directory=None, store_command=None, store_output=None, **kwargs
):
    """Run a command and write its output to a defined location"""
    print(" ".join(str(e) for e in command))
    if output_directory and store_command is None:
        store_command = output_directory / f"{command[0]}.cmd"
    if output_directory and store_output is None:
        store_output = output_directory / f"{command[0]}.log"
    if store_command:
        store_command.parent.mkdir(parents=True, exist_ok=True)
        store_command.write_text(" ".join(str(e) for e in command))
    start = time.perf_counter()
    result = procrunner.run(
        command, environment_override={"DIALS_NOBANNER": "1"}, **kwargs
    )
    print(f"running command took {time.perf_counter() - start:.1f} seconds\n")
    assert not result.returncode, "Command execution failed"
    if store_output:
        store_output.parent.mkdir(parents=True, exist_ok=True)
        store_output.write_bytes(result.stdout)


def generate_processing_detail_text_thaumatin(options):
    print("Generating thaumatin processing tutorial output")

    tmpdir = pathlib.Path("tmp-thaumatin")
    tmpdir.mkdir(exist_ok=True)
    outdir = options.output / "thaumatin"
    outdir.mkdir(parents=True, exist_ok=True)
    runcmd = functools.partial(
        _command_runner, output_directory=outdir, working_directory=tmpdir
    )

    df = dials_data.download.DataFetcher()
    runcmd(["dials.import", df("thaumatin_i04", pathlib=True) / "th_8_2_0*cbf"])
    runcmd(["dials.find_spots", "imported.expt", "nproc=4"])
    runcmd(["dials.index", "imported.expt", "strong.refl"])
    runcmd(["dials.refine_bravais_settings", "indexed.expt", "indexed.refl"])
    runcmd(["dials.reindex", "indexed.refl", "change_of_basis_op=a,b,c"])
    runcmd(
        [
            "dials.refine",
            "bravais_setting_9.expt",
            "reindexed.refl",
            "scan_varying=false",
        ]
    )
    runcmd(
        ["dials.refine", "refined.expt", "refined.refl", "scan_varying=true"],
        store_command=outdir / "dials.sv_refine.cmd",
        store_output=outdir / "dials.sv_refine.log",
    )
    runcmd(["dials.integrate", "refined.expt", "refined.refl", "nproc=4"])
    runcmd(["dials.report", "integrated.expt", "integrated.refl"])
    shutil.copy(tmpdir / "dials.report.html", outdir)
    runcmd(["dials.export", "integrated.refl", "integrated.expt"])

    print(f"Updated result files written to {outdir}")
    if not options.keep:
        shutil.rmtree(tmpdir)


def generate_processing_detail_text_mpro_x0692(options):
    print("Generating mpro-x0692 / PDB ID 5REL processing tutorial output")
    # See: https://www.ebi.ac.uk/pdbe/entry/pdb/5REL
    #      https://doi.org/10.5281/zenodo.3730940

    tmpdir = pathlib.Path("tmp-mpro_x0692")
    tmpdir.mkdir(exist_ok=True)
    outdir = options.output / "mpro_x0692"
    outdir.mkdir(parents=True, exist_ok=True)
    runcmd = functools.partial(
        _command_runner, output_directory=outdir, working_directory=tmpdir
    )

    # Find/validate the data input - until we've decided to integrate this
    # into the main release, have a DLS default or otherwise let it be
    # specified via an environment variable.
    df = dials_data.download.DataFetcher()
    runcmd(["dials.import", df("mpro_x0692", pathlib=True) / "Mpro-x0692_1_0*.cbf"])
    runcmd(["dials.find_spots", "imported.expt", "nproc=4"])
    runcmd(["dials.index", "imported.expt", "strong.refl"])
    # Because it's hard to robustly extract the *final* unindexed count
    # using sphinx's built-in literalinclude, we extract it specially here
    extract_last_indexed_spot_count(
        outdir / "dials.index.log", outdir / "dials.index.log.extract_unindexed"
    )
    runcmd(
        [
            "dials.refine_bravais_settings",
            "indexed.expt",
            "indexed.refl",
            "best_monoclinic_beta=False",
        ]
    )
    cb_op = extract_rbs_cb_op(
        tmpdir / "bravais_summary.json",
        outdir / "dials.refine_bravais_settings.log.cb_op",
        2,
    )
    runcmd(["dials.reindex", "indexed.refl", f"change_of_basis_op={cb_op}"])
    runcmd(
        [
            "dials.refine",
            "bravais_setting_2.expt",
            "reindexed.refl",
            "scan_varying=false",
        ]
    )
    runcmd(
        ["dials.refine", "refined.expt", "refined.refl", "scan_varying=true"],
        store_command=outdir / "dials.sv_refine.cmd",
        store_output=outdir / "dials.sv_refine.log",
    )
    runcmd(["dials.integrate", "refined.expt", "refined.refl", "nproc=4"])
    runcmd(
        [
            "dials.symmetry",
            "integrated.expt",
            "integrated.refl",
            "best_monoclinic_beta=False",
        ]
    )
    runcmd(["dials.scale", "symmetrized.expt", "symmetrized.refl"])
    runcmd(["dials.estimate_resolution", "scaled.expt", "scaled.refl"])
    d_min = extract_resolution(outdir / "dials.estimate_resolution.log", "cc_half")
    runcmd(
        ["dials.scale", "scaled.expt", "scaled.refl", f"d_min={d_min:.2f}"],
        store_command=outdir / "dials.scale_cut.cmd",
        store_output=outdir / "dials.scale_cut.log",
    )
    runcmd(["dials.report", "scaled.expt", "scaled.refl"])
    shutil.copy(tmpdir / "dials.report.html", outdir)

    print(f"Updated result files written to {outdir}")
    if not options.keep:
        shutil.rmtree(tmpdir)


def generate_processing_detail_text_betalactamase(options):
    print("Generating Beta-lactamase processing tutorial output")

    tmpdir = pathlib.Path("tmp-betalactamase")
    tmpdir.mkdir(exist_ok=True)
    outdir = options.output / "betalactamase"
    outdir.mkdir(parents=True, exist_ok=True)
    runcmd = functools.partial(
        _command_runner, output_directory=outdir, working_directory=tmpdir
    )

    # Find/validate the data input - until we've decided to integrate this
    # into the main release, have a DLS default or otherwise let it be
    # specified via an environment variable.
    datadir = os.getenv("BETALACTAMASE_TUTORIAL_DATA")
    if datadir:
        datadir = pathlib.Path(datadir)
    else:
        datadir = pathlib.Path(
            os.getenv(
                "CCP4_TUTORIAL_DATA",
                "/dls/i03/data/2017/mx19576-1/tutorial_data/summed/summed/C2sum_1*.cbf.gz",
            )
        ).parent
    if not datadir.is_dir() or not datadir.glob("C2sum_1*.cbf.gz"):
        sys.exit(
            """Error:  Could not find Betalactamase data: skipping text generation.
        Please download C2sum_1 from https://zenodo.org/record/1014387 and extract,
        then set environment variable BETALACTAMASE_TUTORIAL_DATA to the folder with C2sum_1_*.cbf.gz"""
        )
    print(f"Using data: {datadir}")

    runcmd(["dials.import", datadir / "C2sum_1*.cbf.gz"])
    runcmd(["dials.find_spots", "imported.expt", "nproc=4"])
    runcmd(["dials.index", "imported.expt", "strong.refl"])
    # Because it's hard to robustly extract the *final* unindexed count
    # using sphinx's built-in literalinclude, we extract it specially here
    extract_last_indexed_spot_count(
        outdir / "dials.index.log", outdir / "dials.index.log.extract_unindexed"
    )
    runcmd(["dials.refine_bravais_settings", "indexed.expt", "indexed.refl"])
    runcmd(["dials.reindex", "indexed.refl", "change_of_basis_op=a+b,-a+b,c"])
    runcmd(
        [
            "dials.refine",
            "bravais_setting_2.expt",
            "reindexed.refl",
        ]
    )
    runcmd(["dials.integrate", "refined.expt", "refined.refl", "nproc=4"])
    runcmd(["dials.symmetry", "integrated.expt", "integrated.refl"])
    runcmd(["dials.scale", "symmetrized.expt", "symmetrized.refl"])
    runcmd(
        ["dials.scale", "scaled.expt", "scaled.refl", "d_min=1.4"],
        store_command=outdir / "dials.scale_cut.cmd",
        store_output=outdir / "dials.scale_cut.log",
    )
    runcmd(["dials.report", "scaled.expt", "scaled.refl"])
    shutil.copy(tmpdir / "dials.report.html", outdir)

    print(f"Updated result files written to {outdir}")
    if not options.keep:
        shutil.rmtree(tmpdir)


def generate_multi_crystal_symmetry_and_scaling(options):
    print("Generating multi-crystal symmetry analysis and scaling output")

    tmpdir = pathlib.Path(tempfile.mkdtemp("_multi_crystal", dir="."))
    tmpdir.mkdir(exist_ok=True)
    outdir = options.output / "multi_crystal"
    outdir.mkdir(parents=True, exist_ok=True)
    runcmd = functools.partial(
        _command_runner, output_directory=outdir, working_directory=tmpdir
    )

    df = dials_data.download.DataFetcher()
    experiment_files = sorted(
        df("vmxi_proteinase_k_sweeps", pathlib=True).glob("experiments_*.expt")
    )
    reflection_files = sorted(
        df("vmxi_proteinase_k_sweeps", pathlib=True).glob("reflections_*.refl")
    )
    input_files = []
    for src in experiment_files + reflection_files:
        dst = tmpdir / src.name
        os.symlink(src, dst)
        input_files.append(dst.name)
    runcmd(["xia2.multiplex"] + input_files)
    shutil.copy(tmpdir / "xia2.multiplex.html", outdir)
    runcmd(["dials.cosym"] + input_files)
    shutil.copy(tmpdir / "dials.cosym.html", outdir)
    runcmd(["dials.scale", "symmetrized.expt", "symmetrized.refl"])
    runcmd(["dials.estimate_resolution", "scaled.expt", "scaled.refl"])
    d_min = extract_resolution(outdir / "dials.estimate_resolution.log", "cc_half")
    runcmd(
        ["dials.scale", "scaled.expt", "scaled.refl", f"d_min={d_min:.2f}"],
        store_command=outdir / "dials.scale_cut.cmd",
        store_output=outdir / "dials.scale_cut.log",
    )
    runcmd(["dials.compute_delta_cchalf", "scaled.refl", "scaled.expt"])
    runcmd(
        ["dials.scale", "scaled.expt", "scaled.refl", f"d_min={d_min:.2f}"],
        store_command=outdir / "dials.scale_exclude.cmd",
        store_output=outdir / "dials.scale_exclude.log",
    )
    shutil.copy(tmpdir / "dials.scale.html", outdir)
    runcmd(["dials.symmetry", "scaled.expt", "scaled.refl", "laue_group=None"])
    shutil.copy(tmpdir / "dials.symmetry.html", outdir)
    runcmd(["dials.merge", "symmetrized.expt", "symmetrized.refl"])

    print(f"Updated result files written to {outdir}")
    if not options.keep:
        shutil.rmtree(tmpdir)


def find_in_line(string, lines, start=0):
    """Find the next line index containing a given string"""
    for n, line in enumerate(lines[start:], start):
        if string in line:
            return n
    raise RuntimeError(f"Could not find line '{string}' in lines")


def write_extract(destination, start, end, lines):
    """Write lines to a file, in the correct line position, with markers.

    This can be used to provide sphinx with an easily extractable literalinclude
    block, that preserves the correct line numbers from the original file.
    """
    assert start > 0
    out_lines = []
    for n, line in enumerate(lines):
        if n == start - 1:
            out_lines.append("[START_EXTRACT]")
        elif n >= start and n <= end:
            out_lines.append(line.strip())
        elif n == end + 1:
            out_lines.append("[END_EXTRACT]")
        else:
            out_lines.append("")

    destination.write_text("\n".join(out_lines), "latin-1")


def extract_last_indexed_spot_count(source, destination):
    """Extract the last 'unindexed' entry from an index log."""
    lines = source.read_text("latin-1").split("\n")

    # Find the last entry for '# unindexed'
    next_ui = len(lines) - find_in_line("# unindexed", list(reversed(lines))) - 1

    # We now have the last entry for '# unindexed'. Find the end of the table
    end_ui = find_in_line("---", lines, next_ui + 2)

    # Range to extract = (next_ui, end_ui)
    write_extract(destination, next_ui - 1, end_ui, lines)


def extract_resolution(source, method):
    """Extract the resolution from a dials.estimate_resolution log."""
    lines = source.read_text("latin-1").split("\n")

    # Find the Resolution line
    resolution_line = lines[find_in_line(f"Resolution {method}", lines)]
    # Parse and return the suggested resolution value
    return float(resolution_line.split(":")[-1].strip())


def extract_rbs_cb_op(source, destination, solution):
    cb_op = json.loads(source.read_text())[f"{solution}"]["cb_op"]
    write_extract(destination, 1, 3, [cb_op])
    return cb_op


@dials.util.show_mail_handle_errors()
def run(args=None):
    parser = argparse.ArgumentParser(
        description="Generate tutorial logs for DIALS documentation website"
    )
    parser.add_argument("-?", action="help", help=argparse.SUPPRESS)
    parser.add_argument(
        "--beta",
        action="store_true",
        help="Generate betalactamase tutorial logs",
    )
    parser.add_argument(
        "--mpro_x0692",
        action="store_true",
        help="Generate Mpro x0692 tutorial logs",
    )
    parser.add_argument(
        "--thaum",
        action="store_true",
        help="Generate thaumatin tutorial logs",
    )
    parser.add_argument(
        "--multi_crystal",
        action="store_true",
        help="Generate multi-crystal tutorial logs",
    )
    parser.add_argument(
        "--keep",
        action="store_true",
        help="Keep temporary directories on successful generation",
    )
    parser.add_argument(
        "--output",
        action="store",
        type=pathlib.Path,
        default=".",
        help="Write output to this location",
    )
    options = parser.parse_args(args)

    targets = []
    if options.beta:
        targets.append(generate_processing_detail_text_betalactamase)
    if options.mpro_x0692:
        targets.append(generate_processing_detail_text_mpro_x0692)
    if options.thaum:
        targets.append(generate_processing_detail_text_thaumatin)
    if options.multi_crystal:
        targets.append(generate_multi_crystal_symmetry_and_scaling)

    if not targets:
        parser.error(
            "You need to specify at least one dataset\n"
            "to generate logs for (eg. --beta, --thaum)."
        )

    for target in targets:
        target(options)


if __name__ == "__main__":
    run()
