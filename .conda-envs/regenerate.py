#!/usr/bin/env python3
from __future__ import annotations

import pathlib
import subprocess
import sys

root = pathlib.Path(__file__).parent.parent

for platform in ["macos", "linux", "windows"]:
    command = [
        sys.executable,
        str(root / "util" / "parse_dependency_selectors.py"),
        "--no-channel",
        "--header",
        str(root / "dependencies.yaml"),
        str(root / ".conda-envs" / "cctbx-dependencies.yaml"),
    ]
    data = subprocess.run(
        command,
        capture_output=True,
        check=True,
        text=True,
    )
    (root / ".conda-envs" / f"{platform}.txt").write_text(data.stdout)
