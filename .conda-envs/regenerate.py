#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import shlex
import subprocess
import sys

root = pathlib.Path(__file__).parent.parent

for platform, conda_platform in [
    ("macos", "osx"),
    ("linux", "linux"),
    ("windows", "win"),
]:
    command = [
        sys.executable,
        str(root / "util" / "parse_dependency_selectors.py"),
        "--no-channel",
        "--header",
        "--platform",
        conda_platform,
        str(root / "dependencies.yaml"),
        str(root / ".conda-envs" / "cctbx-dependencies.yaml"),
        str(root / ".conda-envs" / "fallback-dxtbx-dependencies.yaml"),
    ]
    try:
        data = subprocess.run(
            command,
            capture_output=True,
            check=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Error running:\n    {shlex.join(command)}:\n" + e.stdout)
        raise
    (root / ".conda-envs" / f"{platform}.txt").write_text(data.stdout)
