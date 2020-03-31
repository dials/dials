# LIBTBX_SET_DISPATCHER_NAME dev.dials.generate_tutorial_text

from __future__ import absolute_import, division, print_function

import functools
import os
import sys
import tempfile
from optparse import SUPPRESS_HELP, OptionParser

import dials_data.download
import procrunner
import py


def run(
    command, output_directory=None, store_command=None, store_output=None, **kwargs
):
    """Run a command and write its output to a defined location"""
    print(" ".join(str(e) for e in command))
    if output_directory and store_command is None:
        store_command = output_directory.join(command[0] + ".cmd")
    if output_directory and store_output is None:
        store_output = output_directory.join(command[0] + ".log")
    if store_command:
        store_command.write(" ".join(str(e) for e in command), ensure=True)
    result = procrunner.run(
        command, environment_override={"DIALS_NOBANNER": "1"}, **kwargs
    )
    print("running command took {:.1f} seconds\n".format(result["runtime"]))
    assert not result.returncode, "Command execution failed"
    if store_output:
        store_output.write_binary(result.stdout, ensure=True)


def generate_processing_detail_text_thaumatin(options):
    print("Generating thaumatin processing tutorial output")

    tmpdir = py.path.local("./tmp-thaumatin")
    tmpdir.ensure(dir=1)
    outdir = py.path.local(options.output).join("thaumatin")
    runcmd = functools.partial(run, output_directory=outdir, working_directory=tmpdir)

    df = dials_data.download.DataFetcher()
    runcmd(["dials.import", df("thaumatin_i04").join("th_8_2_0*cbf")])
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
    tmpdir.join("dials-report.html").copy(outdir.join("dials-report.html"))
    runcmd(["dials.export", "integrated.refl", "integrated.expt"])

    print("Updated result files written to {}".format(outdir.strpath))
    if not options.keep:
        tmpdir.remove(rec=1)


def generate_processing_detail_text_betalactamase(options):
    print("Generating Beta-lactamase processing tutorial output")

    tmpdir = py.path.local("./tmp-betalactamase")
    tmpdir.ensure(dir=1)
    outdir = py.path.local(options.output).join("betalactamase")
    runcmd = functools.partial(run, output_directory=outdir, working_directory=tmpdir)

    # Find/validate the data input - until we've decided to integrate this
    # into the main release, have a DLS default or otherwise let it be
    # specified via an environment variable.
    datadir = os.getenv("BETALACTAMASE_TUTORIAL_DATA")
    if datadir:
        datadir = py.path.local(datadir)
    else:
        datadir = py.path.local(
            os.getenv(
                "CCP4_TUTORIAL_DATA",
                "/dls/i03/data/2017/mx19576-1/tutorial_data/summed/summed/C2sum_1*.cbf.gz",
            )
        ).dirpath()
    if not datadir.check(dir=1) or not datadir.listdir("C2sum_1*.cbf.gz"):
        sys.exit(
            """Error:  Could not find Betalactamase data: skipping text generation.
        Please download C2sum_1 from https://zenodo.org/record/1014387 and extract,
        then set environment variable BETALACTAMASE_TUTORIAL_DATA to the folder with C2sum_1_*.cbf.gz"""
        )
    print("Using data: {}".format(datadir.strpath))

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
            "scan_varying=false",
        ]
    )
    runcmd(
        ["dials.refine", "refined.expt", "refined.refl", "scan_varying=true"],
        store_command=outdir / "dials.sv_refine.cmd",
        store_output=outdir / "dials.sv_refine.log",
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
    tmpdir.join("dials-report.html").copy(outdir.join("dials-report.html"))

    print("Updated result files written to {}".format(outdir.strpath))
    if not options.keep:
        tmpdir.remove(rec=1)


def generate_multi_crystal_symmetry_and_scaling(options):
    print("Generating multi-crystal symmetry analysis and scaling output")

    tmpdir = py.path.local(tempfile.mkdtemp("_multi_crystal", dir="."))
    tmpdir.ensure(dir=1)
    outdir = py.path.local(options.output).join("multi_crystal")
    runcmd = functools.partial(run, output_directory=outdir, working_directory=tmpdir)

    df = dials_data.download.DataFetcher()
    experiment_files = sorted(
        df("vmxi_proteinase_k_sweeps").visit("experiments_*.expt")
    )
    reflection_files = sorted(
        df("vmxi_proteinase_k_sweeps").visit("reflections_*.refl")
    )
    input_files = []
    for src in experiment_files + reflection_files:
        dst = tmpdir.join(src.basename)
        os.symlink(src.strpath, dst.strpath)
        input_files.append(dst.basename)
    runcmd(["xia2.multiplex"] + input_files)
    tmpdir.join("xia2.multiplex.html").copy(outdir.join("xia2.multiplex.html"))
    runcmd(["dials.cosym"] + input_files)
    tmpdir.join("dials.cosym.html").copy(outdir.join("dials.cosym.html"))
    runcmd(["dials.scale", "symmetrized.expt", "symmetrized.refl"])
    runcmd(["dials.resolutionizer", "scaled.expt", "scaled.refl"])
    d_min = extract_resolution(outdir / "dials.resolutionizer.log", "cc_half")
    runcmd(
        ["dials.scale", "scaled.expt", "scaled.refl", "d_min=%.2f" % d_min],
        store_command=outdir / "dials.scale_cut.cmd",
        store_output=outdir / "dials.scale_cut.log",
    )
    runcmd(["dials.compute_delta_cchalf", "scaled.refl", "scaled.expt"])
    runcmd(
        ["dials.scale", "scaled.expt", "scaled.refl", "d_min=%.2f" % d_min],
        store_command=outdir / "dials.scale_exclude.cmd",
        store_output=outdir / "dials.scale_exclude.log",
    )
    tmpdir.join("dials.scale.html").copy(outdir.join("dials.scale.html"))
    runcmd(["dials.symmetry", "scaled.expt", "scaled.refl", "laue_group=None"])
    tmpdir.join("dials.symmetry.html").copy(outdir.join("dials.symmetry.html"))
    runcmd(["dials.merge", "symmetrized.expt", "symmetrized.refl"])

    print("Updated result files written to {}".format(outdir.strpath))
    if not options.keep:
        tmpdir.remove(rec=1)


def find_in_line(string, lines, start=0):
    """Find the next line index containing a given string"""
    for n, line in enumerate(lines[start:], start):
        if string in line:
            return n
    raise RuntimeError("Could not find line '{}' in lines".format(string))


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

    destination.write_text(u"\n".join(out_lines), "latin-1")


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
    """Extract the resolution from a dials.resolutionizer log."""
    lines = source.read_text("latin-1").split("\n")

    # Find the Resolution line
    resolution_line = lines[find_in_line("Resolution %s" % method, lines)]
    # Parse and return the suggested resolution value
    return float(resolution_line.split(":")[-1].strip())


if __name__ == "__main__":
    parser = OptionParser(
        description="Generate tutorial logs for DIALS documentation website"
    )
    parser.add_option("-?", action="help", help=SUPPRESS_HELP)
    parser.add_option(
        "--beta",
        dest="beta",
        action="store_true",
        default=False,
        help="Generate betalactamase tutorial logs",
    )
    parser.add_option(
        "--thaum",
        dest="thaum",
        action="store_true",
        default=False,
        help="Generate thaumatin tutorial logs",
    )
    parser.add_option(
        "--multi_crystal",
        dest="multi_crystal",
        action="store_true",
        default=False,
        help="Generate multi-crystal tutorial logs",
    )
    parser.add_option(
        "--keep",
        dest="keep",
        action="store_true",
        default=False,
        help="Keep temporary directories on successful generation",
    )
    parser.add_option(
        "--output",
        dest="output",
        action="store",
        type="string",
        default=".",
        help="Write output to this location",
    )
    options, _ = parser.parse_args()

    targets = []
    if options.beta:
        targets.append(generate_processing_detail_text_betalactamase)
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
