# LIBTBX_SET_DISPATCHER_NAME dials.remove_columns

from __future__ import annotations

import fnmatch
import logging
import os
import sys

import iotbx.phil

import dials.util
from dials.util import log
from dials.util.options import ArgumentParser
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.remove_columns")

help_message = """

Remove named columns from a reflection table, for example to reduce the size
of a file on disk by discarding data which are not needed for the analysis to
be performed.

Columns are named with the remove= option. Several columns may be given at
once, separated by commas, or the option may be repeated. Shell-style
wildcards are matched against the columns which are actually present.

Examples::

  dials.remove_columns indexed.refl remove=xyzobs.mm.value,xyzobs.mm.variance \\
    output=small.refl

  dials.remove_columns strong.refl remove=shoebox

  dials.remove_columns integrated.refl remove=xyzobs.mm.*
"""

phil_scope = iotbx.phil.parse(
    """
  remove = None
    .type = strings
    .multiple = True
    .help = "Names of the columns to remove. Several columns may be given in"
            "one go, separated by commas, e.g."
            "remove=xyzobs.mm.value,xyzobs.mm.variance"
            "remove='xyzobs.mm.value xyzobs.mm.variance'"
            "or the option may be repeated. Shell-style wildcards are matched"
            "against the columns present, e.g. remove=xyzobs.mm.*"

  output = stripped.refl
    .type = path
    .help = "The output reflection file"

  log = None
    .type = path
    .help = "The log filename. If unset, output is written to the terminal only."
"""
)

# Columns which most other DIALS programs depend upon. Removing these is
# allowed, but is unlikely to be what was intended, so warn about it.
LOAD_BEARING_COLUMNS = ("id", "panel")


def format_file_size(nbytes):
    """Format a size in bytes for display, e.g. 1234567 -> '1.2 MB'."""

    size = float(nbytes)
    for unit in ("B", "kB", "MB"):
        if size < 1000:
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1000
    return f"{size:.1f} GB"


def expand_column_patterns(commands, keys):
    """Expand the remove= command line options against the columns present.

    :param commands: list of lists of strings, as returned by a multiple phil
                     strings parameter, e.g. [["xyzobs.mm.value,shoebox"]]
    :param keys: the column names present in the reflection table
    :return: (set of matching column names, list of patterns matching nothing)
    """

    patterns = []
    # Note, when extracted rather than parsed, an unset multiple phil strings
    # parameter becomes [None].
    for command in filter(lambda x: x is not None, commands):
        for value in command:
            patterns.extend(p for p in value.split(",") if p)

    matched = set()
    unmatched = []

    for pattern in patterns:
        # fnmatchcase rather than fnmatch, as the latter is case insensitive
        # on platforms with case insensitive filesystems
        found = [k for k in keys if fnmatch.fnmatchcase(k, pattern)]
        if found:
            matched.update(found)
        else:
            unmatched.append(pattern)

    return matched, unmatched


@dials.util.show_mail_handle_errors()
def run(args=None):
    usage = "dials.remove_columns [options] indexed.refl remove=column,column"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        epilog=help_message,
    )

    params, _ = parser.parse_args(args, show_diff_phil=True)

    if len(params.input.reflections) != 1:
        parser.print_help()
        return

    log.config(logfile=params.log)
    logger.info(dials_version())

    filename = params.input.reflections[0].filename
    reflections = params.input.reflections[0].data

    keys = sorted(reflections.keys())

    if all(command is None for command in params.remove):
        parser.print_help()
        sys.exit("\nNo columns to remove were given, e.g. remove=shoebox")

    remove, unmatched = expand_column_patterns(params.remove, keys)

    for pattern in unmatched:
        logger.warning(f"No column matching {pattern} in {filename}")

    if not remove:
        sys.exit("None of the requested columns are present")

    if len(remove) == len(keys):
        sys.exit("Refusing to remove every column from the reflection table")

    load_bearing = [column for column in LOAD_BEARING_COLUMNS if column in remove]
    if load_bearing:
        logger.warning(
            f"Removing {', '.join(load_bearing)}: the output is unlikely to be "
            "usable by other DIALS programs"
        )

    for column in sorted(remove):
        del reflections[column]

    logger.info(f"Removed columns: {', '.join(sorted(remove))}")
    logger.info(f"Retained columns: {', '.join(sorted(reflections.keys()))}")

    reflections.as_file(params.output)

    logger.info(f"Saved {reflections.size()} reflections to {params.output}")
    logger.info(
        f"File size: {format_file_size(os.path.getsize(filename))} -> "
        f"{format_file_size(os.path.getsize(params.output))}"
    )


if __name__ == "__main__":
    run()
