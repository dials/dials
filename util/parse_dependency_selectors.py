#!/usr/bin/env python

from __future__ import annotations

import os
import re
import sys

re_selector = re.compile(r"# *\[([^#]+)]$")
re_pin = re.compile(r"""{{ *pin_compatible *\( *['"]([^'"]+)['"]""")


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
                    print(f"... Passed: {line}")
                    output_lines.append(line)
            elif re_pin.search(line):
                # Ignore pin_compatible dependencies
                continue
            else:
                output_lines.append(line)
        return "\n".join(output_lines)

    def parse_file(self, filename):
        # type: (str) -> dict[str, list[str]]
        """Parse a dependency file into a structured dictionary"""
        with open(filename, "rt") as f:
            data = self.preprocess(f.read())
        output = {}
        current_section = None
        for n, line in enumerate(data.splitlines()):
            print(f"Examining line {n}: {line} (current: {current_section})")
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
                output[current_section].append(line[1:].strip())
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
        help="Choose the target for handling dependency lists.",
        default="bootstrap",
    )
    # parser.add_argument("--conda-build", action="store_const", const="conda-build", dest="kind", help="Run as though constructing a conda-build recipe")
    parser.add_argument(
        "--prebuilt-cctbx", help="Mark as using prebuilt cctbx. Implied by conda-build."
    )
    parser.add_argument("source", nargs="+", help="Dependency files to merge")
    args = parser.parse_args()
    if not args.kind:
        print("Must provide source files")
        parser.print_usage()
        sys.exit(1)

    deps = DependencySelectorParser(
        bootstrap=args.kind == "bootstrap", prebuilt_cctbx=args.prebuilt_cctbx
    )
    pprint(deps.parse_file(args.source[0]))
# from pprint import pprint

# if len(sys.argv) > 1:
#     for arg in sys.argv[1:]:
#         print("Parsing " + arg)
#         pprint(parser.parse_file(arg))
