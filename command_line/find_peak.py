from __future__ import division
from dxtbx.imageset import ImageSetFactory
from dials.algorithms.peak_finding.spot_finder_lui import SpotFinderLui
#from dials.scratch.luiso_s.to_dials_reg.func_fnd_pk import *

import numpy
def fnd_pk():
    import sys
    from dxtbx.format.Registry import Registry

    arg_lst = sys.argv[1:]
    algrm = 'none'
    filenames = []
    times = 5
    shift = 10
    n_blocks_x = 5
    n_blocks_y = 12
    for arg in arg_lst:
        if '=' in arg:
            leng_01 = arg.find('=')
            lft_str = arg[:leng_01].lower()
            print lft_str
            if lft_str == 'algr' or lft_str == 'al':
                if arg[leng_01 + 1:] == 'lui':
                    algrm = 'lui'
                else:
                    algrm = 'xds'
            elif lft_str == 'times' or lft_str == 'tm':
                times = int(arg[leng_01 + 1:])
            elif lft_str == 'shf' or lft_str == 'shift':
                shift = int(arg[leng_01 + 1:])
            elif lft_str == 'nbx' or lft_str == 'numblockx':
                n_blocks_x = int(arg[leng_01 + 1:])
            elif lft_str == 'nby' or lft_str == 'numblocky':
                n_blocks_y = int(arg[leng_01 + 1:])

        else:
            filenames.append(arg)
    print len(filenames), "images given"
    print 'following', algrm, 'algorithm'

    print 'algrm =', algrm

    print'___________________________'
    print 'filenames =', filenames

    if algrm == 'none':
        print 'no algorithm spesifyed'
    else:
        if algrm == 'lui':
            sweep = ImageSetFactory.new(filenames)
            assert(len(sweep) == 1)
            sweep = sweep[0]
            find_spots = SpotFinderLui()
            reflection_list = find_spots(sweep, times, shift, n_blocks_x, n_blocks_y)

            import pickle
            output_file = 'lui_reflections.pkl'
            pickle.dump(reflection_list, open(output_file, 'wb'))
        elif algrm == 'xds':
            print 'xds algorithm needs to be connected'

if __name__ == '__main__':
    fnd_pk()
