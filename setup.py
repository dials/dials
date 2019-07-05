#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from setuptools import setup, find_packages

requirements = ["dials-data", "Jinja2", "procrunner", "six"]

setup_requirements = []
needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
if needs_pytest:
    setup_requirements.append("pytest-runner")

test_requirements = ["mock", "pytest"]

setup(
    author="Diamond Light Source",
    author_email="scientificsoftware@diamond.ac.uk",
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: Implementation :: CPython",
    ],
    description="Diffraction Integration for Advanced Light Sources",
    install_requires=requirements,
    license="BSD license",
    include_package_data=True,
    keywords="dials",
    name="dials",
    packages=find_packages(),
    package_dir={"dials": "../dials"},
    data_files=[
        ("dials", ["conftest.py", "__init__.py", "libtbx_refresh.py", "run_tests.py"])
    ],
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/dials/dials",
    version="0.0.1",
    zip_safe=False,
)
