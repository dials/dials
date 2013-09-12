from __future__ import division

def print_profile(r):
    s = r.shoebox
    b = r.bounding_box
    _i, _j, _k = s.all()
    for i in range(_i):
        for j in range(_j):
            for k in range(_k):
                print '%5d' % int(s[i,j,k]),
            print ''
        print ''
        print ''

def show_profiles(integrated_pickle, isig_limit = None):
    from dials.model.data import ReflectionList # implicit import
    import cPickle as pickle
    import math

    integrated_data = pickle.load(open(integrated_pickle, 'r'))

    for j, r in enumerate(integrated_data):
        if not isig_limit is None:
            if r.intensity <= 0:
                continue
            if r.intensity / math.sqrt(r.intensity_variance) < isig_limit:
                continue

            print_profile(r)

if __name__ == '__main__':
    import sys
    show_profiles(sys.argv[1], isig_limit = float(sys.argv[2]))
