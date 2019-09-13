# LIBTBX_SET_DISPATCHER_NAME dev.dials.make_sphinx_html

from __future__ import absolute_import, division, print_function

import json
import os
import re
import sys
from datetime import datetime
from optparse import SUPPRESS_HELP, OptionParser

import dials
import procrunner
import py

# Disable all HTTPS verification. This is to work around an issue
# in biopython, possibly biopython relying on unreliable servers.
os.environ["PYTHONHTTPSVERIFY"] = "0"


try:  # Python 3
    from urllib.request import urlopen, Request
except ImportError:  # Python 2
    from urllib2 import urlopen, Request


def update_dials_download_links():
    print("Checking DIALS release status")
    url_request = Request("https://api.github.com/repos/dials/dials/releases/latest")
    conn = urlopen(url_request)
    try:
        release_info = json.load(conn)
    finally:
        conn.close()

    dials_dir = py.path.local(dials.__file__).dirpath()
    release_file = dials_dir / "doc" / "sphinx" / "installation.stable_release"

    with release_file.open("w") as release:
        caption = "Stable Release"
        if "name" in release_info:
            caption = caption + ": " + release_info["name"]
            print("Most recent major DIALS release is:", release_info["name"])
        else:
            print("Could not determine most recent major DIALS release")
        release.write(caption + "\n" + "=" * len(caption) + "\n\n")

        release.write(
            "The current stable release can be downloaded from `Github <https://github.com/dials/dials/releases/latest>`_,\n"
        )
        release.write(
            "where you can also find further `release notes <https://github.com/dials/dials/releases/latest>`_.\n\n"
        )

        def download_button(text, version, link):
            print("  %s %s -> %s" % (version, text, link))
            return ".. button::\n   :text: DIALS %s %s\n   :link: %s\n\n" % (
                version,
                text,
                link,
            )

        assets = {}
        for a in release_info.get("assets", []):
            tag = re.search("dials-v([^-]+)-([^-]+)-([^-]+)-(.+)", a["name"])
            if tag:
                shortname = tag.group(4)
                version = ".".join(tag.group(1, 2, 3))
                last_update = datetime.strptime(
                    a["updated_at"], "%Y-%m-%dT%H:%M:%SZ"
                )  # - datetime(1970,1,1)).total_seconds()
                if shortname not in assets or assets[shortname][0] < last_update:
                    assets[shortname] = (
                        last_update,
                        version,
                        a.get("browser_download_url"),
                    )

        long_names = {
            "macosx.pkg": "Mac installer",
            "macosx.tar.gz": "Mac tar archive",
            "macosx-10.6.pkg": "Mac installer (OS X 10.6)",
            "macosx-10.6.tar.gz": "Mac tar archive (OS X 10.6)",
            "linux-x86_64.tar.xz": "Linux installer",
            "source.tar.xz": "Source installer",
        }

        buttons = [
            download_button(long_names.get(asset, asset), _version, link)
            for asset, (_, _version, link) in assets.items()
        ]

        release.write("".join(sorted(buttons)))


if __name__ == "__main__":
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
        default=None,
        help="Copy generated output to this location",
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
        "--no-clean",
        dest="clean",
        action="store_false",
        default=True,
        help="Don't run 'make clean' before building the documentation",
    )
    parser.add_option(
        "--parallel",
        dest="parallel",
        action="store_true",
        default=False,
        help="Build documentation in parallel",
    )
    options, _ = parser.parse_args()

    try:
        update_dials_download_links()
    except Exception as e:
        if options.strict:
            raise
        print("Ignoring error:", e)
    dials_dir = py.path.local(dials.__file__).dirpath()
    sphinx_dir = dials_dir / "doc" / "sphinx"
    tutorial_doc_dir = sphinx_dir / "documentation" / "tutorials"

    sphinx_options = ""
    if options.parallel:
        print("DIALS 1.14 variant does not support parallel building")
    if options.strict:
        sphinx_options += " -W"
    if options.legacy_version:
        if " " in options.legacy_version:
            exit("legacy version must not contain spaces")
        sphinx_options += " -A legacy_version=" + options.legacy_version
    if options.logs:
        sphinx_options += " -D dials_logs=" + options.logs
        for report in ("betalactamase", "thaumatin"):
            py.path.local(options.logs).join(report).join("dials-report.html").copy(
                tutorial_doc_dir.join(report + "-report.html")
            )
    else:
        sys.exit(
            "You must specify the location of the tutorial data output with the '-l' option"
        )
    env = {"SPHINXOPTS": sphinx_options}

    # Disable all HTTPS verification. This is to work around an issue
    # in biopython, possibly biopython relying on unreliable servers.
    # env["PYTHONHTTPSVERIFY"] = "0"

    if options.clean:
        result = procrunner.run(
            ["make", "clean"], environment_override=env, working_directory=sphinx_dir.strpath
        )
        assert not result["exitcode"], (
            "make clean failed with exit code %d" % result["exitcode"]
        )
    result = procrunner.run(
        ["make", "html"], environment_override=env, working_directory=sphinx_dir.strpath
    )
    assert not result["exitcode"], (
        "make html failed with exit code %d" % result["exitcode"]
    )
    if options.output:
        print("Copying HTML pages to", options.output)
        dest_dir = py.path.local(options.output)
        sphinx_dir.join("build").join("html").copy(dest_dir)
