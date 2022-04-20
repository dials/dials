"""
Handle dynamic aspects of setup.py and building.

This is separate because the non-dynamic stuff can generally be moved
out of a setup.py, but mainly because at the moment it's how poetry
offloads the unresolved build phases.
"""
from __future__ import annotations

import ast
import itertools
import re
import sys
from pathlib import Path
from typing import Any


def get_entry_point(filename: Path, prefix: str, import_path: str) -> list[str]:
    """Returns the entry point string for a given path.

    This looks for LIBTBX_SET_DISPATCHER_NAME, and a root function
    named 'run'. It can return multiple results for each file, if more
    than one dispatcher name is bound.

    Args:
        filename:
            The python file to parse. Will look for a run() function
            and any number of LIBTBX_SET_DISPATCHER_NAME.
        prefix: The prefix to output the entry point console script with
        import_path: The import path to get to the package the file is in

    Returns:
        A list of entry_point specifications
    """
    contents = filename.read_text(encoding="utf-8")
    tree = ast.parse(contents)
    # Find root functions named "run"
    has_run = any(
        x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == "run"
    )
    if not has_run:
        return []
    # Find if we need an alternate name via LIBTBX_SET_DISPATCHER_NAME
    alternate_names = re.findall(
        r"^#\s*LIBTBX_SET_DISPATCHER_NAME\s+(.*)$", contents, re.M
    )
    if alternate_names:
        return [f"{name}={import_path}.{filename.stem}:run" for name in alternate_names]

    return [f"{prefix}.{filename.stem}={import_path}.{filename.stem}:run"]


def build(setup_kwargs: dict[str, Any]) -> None:
    """Called by setup.py to inject any dynamic configuration"""
    package_path = Path(__file__).parent / "src" / "dials"
    entry_points = setup_kwargs.setdefault("entry_points", {})
    console_scripts = entry_points.setdefault("console_scripts", [])
    # Work out what dispatchers to add
    all_dispatchers = sorted(
        itertools.chain.from_iterable(
            get_entry_point(f, "dials", "dials.command_line")
            for f in (package_path / "command_line").glob("*.py")
        )
    )
    console_scripts.extend(x for x in all_dispatchers if x not in console_scripts)
    libtbx_dispatchers = entry_points.setdefault("libtbx.dispatcher.script", [])
    libtbx_dispatchers.extend(
        "{name}={name}".format(name=x.split("=")[0]) for x in console_scripts
    )

    print(f"Found {len(entry_points['console_scripts'])} dials dispatchers")


if __name__ == "__main__":
    sys.exit("Cannot call build.py directly, please use setup.py instead")
