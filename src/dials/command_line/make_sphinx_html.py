# LIBTBX_SET_DISPATCHER_NAME dev.dials.make_sphinx_html

from __future__ import annotations

import os
import sys
from optparse import SUPPRESS_HELP, OptionParser

import procrunner
import py

import dials.util


@dials.util.show_mail_handle_errors()
def run(args=None):
    dials_dir = py.path.local(dials.__file__).dirpath()
    sphinx_dir = dials_dir / "doc" / "sphinx"
    tutorial_doc_dir = sphinx_dir / "documentation" / "tutorials"

    parser = OptionParser(description="Generate documentation website for DIALS")
    parser.add_option("-?", action="help", help=SUPPRESS_HELP)
    parser.add_option(
        "-s",
        "--strict",
        dest="strict",
        action="store_true",
        default=True,
        help="Run in strict mode and stop on encountering any errors or warnings (default)",
    )
    parser.add_option(
        "-i",
        "--ignore",
        dest="strict",
        action="store_false",
        help="Ignore any errors or warnings",
    )
    parser.add_option(
        "-l",
        "--logs",
        dest="logs",
        action="store",
        type="string",
        default=None,
        help="Use generated dials output logs from this location",
    )
    parser.add_option(
        "-o",
        "--output",
        dest="output",
        action="store",
        type="string",
        default=(sphinx_dir / "build" / "html").strpath,
        help="Generate output in this location",
    )
    parser.add_option(
        "--legacy-version",
        dest="legacy_version",
        action="store",
        type="string",
        default=None,
        metavar="VERSION",
        help="Add a warning message to every page saying that this documentation "
        "relates to the out-of-date version VERSION",
    )
    parser.add_option(
        "--clean",
        dest="clean",
        action="store_true",
        default=False,
        help="Empty the output directory before building the documentation",
    )
    parser.add_option(
        "--parallel",
        dest="parallel",
        action="store_true",
        default=False,
        help="Build documentation in parallel",
    )
    options, _ = parser.parse_args(args)

    output_dir = py.path.local(options.output)
    if options.clean:
        print("Cleaning out", output_dir.strpath)
        output_dir.remove()
    print("Generating documentation into", output_dir.strpath)

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
        if log_dir.check(dir=1):
            options.logs = log_dir.strpath
            print("Using DIALS logs from", options.logs)

    if options.logs:
        command.extend(["-D", "dials_logs=" + os.path.abspath(options.logs)])
        for report in ("betalactamase", "thaumatin"):
            py.path.local(options.logs).join(report).join("dials.report.html").copy(
                tutorial_doc_dir.join(report + "-report.html")
            )
    else:
        sys.exit(
            "You must specify the location of the tutorial data output with the '-l' option"
        )

    env = {}
    # Disable all HTTPS verification. This is to work around an issue
    # in biopython, possibly biopython relying on unreliable servers.
    # env["PYTHONHTTPSVERIFY"] = "0"

    result = procrunner.run(
        command, environment_override=env, working_directory=sphinx_dir
    )
    if result.returncode:
        sys.exit("Sphinx build failed with exit code %d" % result.returncode)


if __name__ == "__main__":
    run()
