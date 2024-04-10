#!/usr/bin/env python

from __future__ import annotations

import os
import re
import sys

re_selector = re.compile(r"# *\[([^#]+)]$")
re_pin = re.compile(r"""{{ *pin_compatible *\( *['"]([^'"]+)['"]""")


def _split_dependency_line(line):
    # type: (str) -> tuple[str | None, str | None, str]

    # Lines that are templated get ignored here
    if "{" in line:
        return (None, None, line)
    pending = line
    # Strip off the comment/selector
    if "#" in line:
        pending = pending[: pending.index("#")].strip()
    # If we have a version spec and no space, this is an error
    if " " not in pending and (set(pending) & set("><=!")):
        raise RuntimeError(
            "Error: Versioned requirement '%s' has no space" % (pending,)
        )
    vers = None
    if " " in pending:
        pending, vers = pending.split(" ", maxsplit=1)
    return (pending, vers, line)


class DependencySelectorParser(object):
    """
    Parse simple conda-build selectors syntax, with optional variables.

    Supported:
    - Variables linux, osx, win, in addition to anything passed into __init__
    - Variable inversion e.g. "not osx"
    - Basic "And" combinations e.g. "bootstrap and not osx"
    """

    def __init__(self, **kwargs):
        self._vars = dict(kwargs)
        self._vars.update(
            {
                "osx": sys.platform == "darwin",
                "linux": sys.platform.startswith("linux"),
                "win": os.name == "nt",
            }
        )

    def _parse_fragment(self, fragment, pos=0):
        # type: (str, int) -> bool
        """Recursively parse a single expression or fragment of an expression."""
        # print(f"Parsing fragment: {fragment}")
        # if "not bootstrap" in fragment:
        #     breakpoint()
        if fragment in self._vars:
            return self._vars[fragment]
        if " and " in fragment:
            left, right = fragment.split(" and ", maxsplit=1)
            return self._parse_fragment(left, pos) and self._parse_fragment(
                right, pos + fragment.index(" and ")
            )
        if fragment.startswith("not "):
            return not self._parse_fragment(fragment[4:].strip(), pos + 4)
        raise ValueError("Could not parse selector fragment '" + fragment + "'")

    def preprocess(self, data):
        # type: (str) -> str
        """Apply preprocessing selectors to file data"""
        output_lines = []
        for line in data.splitlines():
            match = re_selector.search(line)

            if match:
                if self._parse_fragment(match.group(1)):
                    # print(f"... Passed: {line}")
                    output_lines.append(line)
            elif re_pin.search(line):
                # Ignore pin_compatible dependencies
                continue
            else:
                output_lines.append(line)
        return "\n".join(output_lines)

    def parse_file(self, filename):
        # type: (str) -> dict[str, list[tuple[str, str|None, str]]]
        """
        Parse a dependency file into a structured dictionary.

        The dictionary has structure:
        {
            "section": [
                ("dependency_name", "dependency_version", "raw_line"),
                ...
            ]
        }
        """
        with open(filename, "rt") as f:
            data = self.preprocess(f.read())
        output = {}
        current_section = None
        for n, line in enumerate(data.splitlines()):
            # print(f"Examining line {n}: {line} (current: {current_section})")
            if "#" in line:
                line = line[: line.index("#")]
            line = line.strip()
            if line.endswith(":"):
                current_section = line[:-1].strip()
                output[current_section] = []
            elif line.startswith("-"):
                if not current_section:
                    raise RuntimeError(
                        "Error parsing "
                        + filename
                        + ":"
                        + str(n + 1)
                        + "; No current section on line '"
                        + line
                        + "'"
                    )
                output[current_section].append(_split_dependency_line(line[1:].strip()))
            else:
                if line:
                    raise RuntimeError(
                        "Error parsing "
                        + filename
                        + ":"
                        + str(n + 1)
                        + "; Uncategorised line '"
                        + line
                        + "'"
                    )
        return output


def test_parser():
    parser = DependencySelectorParser(bootstrap=True, prebuilt_cctbx=False)
    assert parser._parse_fragment("osx")
    assert parser._parse_fragment("bootstrap")
    assert parser._parse_fragment("osx and bootstrap")
    assert not parser._parse_fragment("linux and bootstrap")
    assert not parser._parse_fragment("prebuilt_cctbx and osx and not bootstrap")


if __name__ == "__main__":
    from argparse import ArgumentParser
    from pprint import pprint

    parser = ArgumentParser()
    parser.add_argument(
        "-k",
        "--kind",
        choices=["bootstrap", "conda-build"],
        help="Choose the target for handling dependency lists. Default: %(default)s",
        default="bootstrap",
    )
    # parser.add_argument("--conda-build", action="store_const", const="conda-build", dest="kind", help="Run as though constructing a conda-build recipe")
    parser.add_argument(
        "--prebuilt-cctbx", help="Mark as using prebuilt cctbx. Implied by conda-build."
    )
    parser.add_argument("source", nargs="+", help="Dependency files to merge")
    args = parser.parse_args()
    # if not args.kind:
    #     print("Must provide source files")
    #     parser.print_usage()
    #     sys.exit(1)

    # print("Processing requirements for target:", args.kind)
    deps = DependencySelectorParser(
        bootstrap=args.kind == "bootstrap", prebuilt_cctbx=args.prebuilt_cctbx
    )

    def _merge_deps(source, merge_into):
        # type: (list[tuple[str|None,str|None, str]], list[tuple[str|None,str|None, str]]) -> None
        indices = {x[0]: i for i, x in enumerate(merge_into)}
        for pkg, ver, line in source:
            if pkg is None:
                # Lines that don't define a package always get added
                merge_into.append(pkg, ver, line)
            elif pkg in indices:
                # This already exists in the target. Should we replace it?
                other_ver = merge_into[indices[pkg]][1]
                if not other_ver and ver:
                    print(f"Merging '{line}' over {merge_into[indices[pkg]]}")
                    merge_into[indices[pkg]] = (pkg, ver, line)
                elif other_ver and ver and ver != other_ver:
                    raise RuntimeError(
                        "Cannot merge requirements for %s: '%s' and '%s'"
                        % (pkg, ver, other_ver)
                    )
            else:
                merge_into.append((pkg, ver, line))
                indices[pkg] = len(merge_into) - 1

    # Map of section to list of (dependency_name, dependency_version, raw_line)
    reqs = {}  # type: dict[str, list[tuple[str, str|None, str]]]
    for source in args.source:
        source_reqs = deps.parse_file(source)
        # Now, merge this into the existing requirements
        for section, items in source_reqs.items():
            _merge_deps(items, reqs.setdefault(section, []))
    # breakpoint()
    # If bootstrap, then we further merge everything down
    if args.kind == "bootstrap":
        merged_req = []
        for items in reqs.values():
            _merge_deps(items, merged_req)
        for pkg, ver, _ in sorted(merged_req, key=lambda x: x[0]):
            if pkg == "python":
                # Bootstrap handles this dependency implicitly
                continue
            print(f"conda-forge::{pkg}" + (f"{ver}" if ver else ""))

    else:
        pprint(reqs)
        # for item in items:
        #     _merge_deps(item)
        # for pkg, ver, line in items:
        #     if
    #             if " " not in item and (set(item) & set("")):
    #                 raise RuntimeError("Error: Versioned requirement has no space")
    #             package_ver =
    #             if " " in item:
    #                 # We have a version.

    # pprint(reqs)
# from pprint import pprint

# if len(sys.argv) > 1:
#     for arg in sys.argv[1:]:
#         print("Parsing " + arg)
#         pprint(parser.parse_file(arg))
