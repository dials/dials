from __future__ import absolute_import, division, print_function

import re
from setuptools import find_packages, setup

with open("_precommitbx/__init__.py") as fh:
    version_match = re.search(r"^__version__ = \"([0-9.]+)\"", fh.read(), re.M)
    if not version_match:
        raise RuntimeError("Unable to find version string")
version = version_match.group(1)

setup(
    author="DIALS",
    author_email="dials-support@lists.sourceforge.net",
    description="A cctbx-adapter for pre-commit",
    install_requires=["pre-commit"],
    python_requires=">=3.5",
    license="BSD license",
    include_package_data=True,
    name="precommitbx",
    packages=find_packages(include=["_precommitbx"]),
    version=version,
    zip_safe=False,
)
