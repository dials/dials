from __future__ import absolute_import, division, print_function

from dials.util.command_line import OptionParser
from dials.command_line.show import show_experiments
from libtbx.phil import parse

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

    # only at this point do we actually read any data
    experiments = op.input_experiments_as_data()

    print(show_experiments(experiments))


if __name__ == "__main__":
    main()
