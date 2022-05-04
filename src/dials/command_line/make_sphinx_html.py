# LIBTBX_SET_DISPATCHER_NAME dev.dials.make_sphinx_html

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import procrunner

import dials.util


@dials.util.show_mail_handle_errors()
def run(args=None):
    # Navigate to the DIALS root
    dials_dir = Path(dials.__file__).parent.parent.parent
    sphinx_dir = dials_dir / "doc" / "sphinx"
    if not sphinx_dir.is_dir():
        sys.exit(
            "Error: Could not find DIALS sphinx folder at ../../doc/sphinx - are you running from a checkout?"
        )
    tutorial_doc_dir = sphinx_dir / "documentation" / "tutorials"

    parser = argparse.ArgumentParser(
        description="Generate documentation website for DIALS",
        prog="dev.dials.make_sphinx_html",
    )
    parser.add_argument("-?", action="help", help=argparse.SUPPRESS)
    parser.add_argument(
        "-s",
        "--strict",
        action="store_true",
        default=True,
        help="Run in strict mode and stop on encountering any errors or warnings (default)",
    )
    parser.add_argument(
        "-i",
        "--ignore",
        dest="strict",
        action="store_false",
        help="Ignore any errors or warnings",
    )
    parser.add_argument(
        "-l",
        "--logs",
        dest="logs",
        type=Path,
        help="Use generated dials output logs from this location.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=sphinx_dir / "build" / "html",
        help="Generate output in this location. Defaults to dials doc/sphinx dir.",
    )
    parser.add_argument(
        "--legacy-version",
        metavar="VERSION",
        help="Add a warning message to every page saying that this documentation "
        "relates to the out-of-date version VERSION",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Empty the output directory before building the documentation",
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Build documentation in parallel",
    )
    options = parser.parse_args(args)
    print(options)

    output_dir = Path(options.output).resolve()
    if options.clean:
        print(f"Cleaning out {output_dir}")
        shutil.rmtree(output_dir)
    print(f"Generating documentation into {output_dir}")

    command = ["libtbx.sphinx-build", "-b", "html", ".", output_dir]
    if options.parallel:
        command.extend(["-j", "auto"])
    if options.strict:
        command.append("-W")
    if options.legacy_version:
        command.extend(["-A", "legacy_version=" + options.legacy_version])
    if not options.logs:
        # Try to automatically determine this
        branch = "master"
        if options.legacy_version:
            branch = "dials-" + options.legacy_version
        log_dir = dials_dir / ".." / "dials.github.io" / "tutorial_data" / branch
        if log_dir.is_dir():
            options.logs = Path(log_dir)
            print(f"Using DIALS logs from {options.logs}")

    if options.logs:
        command.extend(["-D", f"dials_logs={options.logs.resolve()}"])
        for report in ("betalactamase", "thaumatin"):
            shutil.copy(
                Path(options.logs) / report / "dials.report.html",
                tutorial_doc_dir / f"{report}-report.html",
            )
    else:
        sys.exit(
            "You must specify the location of the tutorial data output with the '-l' option"
        )

    env = {}
    result = procrunner.run(
        command, environment_override=env, working_directory=sphinx_dir
    )
    if result.returncode:
        sys.exit(f"Sphinx build failed with exit code {result.returncode}")


if __name__ == "__main__":
    run()
