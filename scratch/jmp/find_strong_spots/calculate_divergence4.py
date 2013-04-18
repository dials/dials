

if __name__ == '__main__':


    import os
    import libtbx.load_env
    from dials.util.nexus import NexusFile
    from glob import glob
    from dxtbx.sweep import SweepFactory
    from math import pi
    from scitbx import matrix

    try:
        dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
        print 'FAIL: dials_regression not configured'
        raise

    # The XDS values
    xds_sigma_d = 0.060
    xds_sigma_m = 0.154

    # Read the reflections file
    filename = os.path.join(dials_regression,
        'centroid_test_data', 'reflections.h5')
    #filename = os.path.join('/home/upc86896/Data/X1_strong_M1S1_1_', 'reflections.h5')
    handle = NexusFile(filename, 'r')

    # Get the reflection list
    print 'Reading reflections.'
    predicted = handle.get_reflections()
    print 'Read {0} reflections.'.format(len(predicted))

    # Read images
    template = os.path.join(dials_regression,
        'centroid_test_data', 'centroid_*.cbf')
    #template = os.path.join('/home/upc86896/Data/X1_strong_M1S1_1_', 'X1_strong_M1S1_1_*.cbf')
    filenames = glob(template)

    # Load the sweep
    print 'Loading sweep'
    sweep = SweepFactory.sweep(filenames)
    print 'Loaded sweep of {0} images.'.format(len(sweep))

    from dials.algorithms.peak_finding.spot_finder import SpotFinder
    sweep.reader().set_max_cache(1)
    sweep.get_detector().set_trusted_range((0, 20000))
    find_spots = SpotFinder()

    observed = find_spots(sweep)

    from divergence_and_mosaicity import BeamDivergenceAndMosaicity

    calculate_params = BeamDivergenceAndMosaicity(sweep)
    calculate_params(observed, predicted)
