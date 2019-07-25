# LIBTBX_SET_DISPATCHER_NAME dev.dials.make_sphinx_html

from __future__ import absolute_import, division, print_function

import os
import shutil
from optparse import SUPPRESS_HELP, OptionParser

import libtbx.load_env
import procrunner

# Disable all HTTPS verification. This is to work around an issue
# in biopython, possibly biopython relying on unreliable servers.
os.environ["PYTHONHTTPSVERIFY"] = "0"


def recursive_overwrite(src, dest, ignore=None):
    if os.path.isdir(src):
        if not os.path.isdir(dest):
            os.makedirs(dest)
        files = os.listdir(src)
        if ignore is not None:
            ignored = ignore(src, files)
        else:
            ignored = set()
        for f in files:
            if f not in ignored:
                recursive_overwrite(os.path.join(src, f), os.path.join(dest, f), ignore)
    else:
        shutil.copyfile(src, dest)


def update_dials_download_links():
    dials_dir = libtbx.env.find_in_repositories("dials")
    release_file = os.path.join(
        dials_dir, "doc", "sphinx", "installation.stable_release"
    )
    release_json = os.path.join(
        dials_dir, "doc", "sphinx", "installation.stable_release.json"
    )

    release_info = None
    from libtbx.auto_build.bootstrap import Toolbox
    from datetime import datetime
    import json
    import re

    print("Checking DIALS release status: ", end="")
    if Toolbox().download_to_file(
        "https://api.github.com/repos/dials/dials/releases/latest",
        release_json,
        cache=False,
    ):
        with open(release_json, "r") as json_data:
            release_info = json.load(json_data)
    try:
        os.remove(release_json)
    except OSError:
        pass

    if not release_info:
        release_info = {}

    with open(release_file, "w") as release:
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
            download_button(long_names.get(asset, asset), version, link)
            for asset, (_, version, link) in assets.items()
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
        "--no-clean",
        dest="clean",
        action="store_false",
        default=True,
        help="Don't run 'make clean' before building the documentation",
    )
    options, args = parser.parse_args()

    try:
        update_dials_download_links()
    except Exception as e:
        if options.strict:
            raise
        print("Ignoring error:", e)
    dials_dir = libtbx.env.find_in_repositories("dials")
    dials_github_io = libtbx.env.find_in_repositories("dials.github.io")
    assert (
        dials_github_io is not None
    ), "Repository dials.github.io needs to be checked out and configured in modules directory"
    dest_dir = dials_github_io
    os.chdir(os.path.join(dials_dir, "doc", "sphinx"))

    env = {}
    if options.strict:
        env["SPHINXOPTS"] = "-W"

    if options.clean:
        result = procrunner.run(["make", "clean"], environment_override=env)
        assert not result.returncode, (
            "make clean failed with exit code %d" % result.returncode
        )
    result = procrunner.run(["make", "html"], environment_override=env)
    assert not result.returncode, (
        "make html failed with exit code %d" % result.returncode
    )
    print("Copying HTML pages to", dest_dir)
    recursive_overwrite("build/html", dest_dir)
