from __future__ import division
from dials.algorithms.centroid.filtered_centroid import FilteredCentroid

def tst_filtered_centroid():
    import libtbx.load_env
    try:
        dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
        print 'FAIL: dials_regression not configured'
        return

    import os

    frames = [os.path.join(dials_regression, 'centroid_test_data',
                           'centroid_%04d.cbf' % j) for j in range(1, 10)]

    from dxtbx.imageset import ImageSetFactory

    sweep = ImageSetFactory.new(frames)
    assert(len(sweep) == 1)
    sweep = sweep[0]

    from dials.model.data import Reflection, ReflectionList
    ref = Reflection()
    ref.bounding_box = (1070, 1083, 1022, 1036, 0, 8)

    ref.shoebox = sweep.to_array(
        (ref.bounding_box[4], ref.bounding_box[5],
         ref.bounding_box[2], ref.bounding_box[3],
         ref.bounding_box[0], ref.bounding_box[1])).as_double()

    reflections = ReflectionList()
    reflections.append(ref)

    tc = FilteredCentroid(reflections, n_sigma=3)

    centroids = tc.get_reflections()

    print 'OK'
    #for r in reflections:
    #    print r

if __name__ == '__main__':
    tst_filtered_centroid()
