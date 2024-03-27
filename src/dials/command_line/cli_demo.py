from __future__ import annotations

import time

from libtbx.phil import parse

from dials.command_line.show import show_experiments
from dials.util.command_line import OptionParser

master_scope = parse(
    """
verbose = False
  .type = bool
"""
)


def main():
    op = OptionParser(phil=master_scope)

    print(op)

    if op.params().verbose:
        for input_file in sorted(op.input_experiments()):
            print(input_file)

    # only at this point do we actually read any data - and only first time,
    # after this it is cached in the object not a global scope

    t0 = time.time()

    experiments = op.experiments
    print(show_experiments(experiments))

    t1 = time.time()

    experiments = op.experiments
    print(show_experiments(experiments))

    t2 = time.time()

    print("%.4fs %.4fs" % (t1 - t0, t2 - t1))


if __name__ == "__main__":
    main()
