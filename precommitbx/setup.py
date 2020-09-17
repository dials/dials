from __future__ import absolute_import, division, print_function

import re

from setuptools import find_packages, setup

with open("_precommitbx/__init__.py") as fh:
    init_file = fh.read()
    version_match = re.search(r"^__version__ = \"([0-9.]+)\"", init_file, re.M)
    pc_version_match = re.search(
        r"^__precommit_min_version__ = \"([0-9.]+)\"", init_file, re.M
    )
    if not version_match or not pc_version_match:
        raise RuntimeError("Unable to find version string")

version = version_match.group(1)
pre_commit_min_version = pc_version_match.group(1)

setup(
    author="DIALS",
    author_email="dials-support@lists.sourceforge.net",
    description="A cctbx-adapter for pre-commit",
    install_requires=["pre-commit>=" + pre_commit_min_version],
    python_requires=">=3.5",
    license="BSD license",
    include_package_data=True,
    name="precommitbx",
    packages=find_packages(include=["_precommitbx"]),
    version=version,
    zip_safe=False,
)
