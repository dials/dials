from __future__ import division
from toytbx import make_list

def main(args):
    assert(args)
    n = int(args[0])
    print make_list(n)

if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
