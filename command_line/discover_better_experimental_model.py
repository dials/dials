from __future__ import absolute_import, division, print_function

from dials.command_line.search_beam_position import *

print("*" * 74)
print("*      dials.discover_better_experimental_model has been renamed to      *")
print("*      dials.search_beam_position. Please use that command instead.      *")
print("*" * 74)
print()

if __name__ == '__main__':
  import sys
  try:
    run(sys.argv[1:])
  finally:
    print()
    print("*" * 74)
    print("*      dials.discover_better_experimental_model has been renamed to      *")
    print("*      dials.search_beam_position. Please use that command instead.      *")
    print("*" * 74)
