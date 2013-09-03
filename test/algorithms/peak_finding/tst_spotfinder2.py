from __future__ import division

class Test:

    def __init__(self):
        from dials.algorithms.peak_finding import XDSThresholdStrategy
        from dials.algorithms.peak_finding import FindSpots
        from dials.model.serialize import load
        import libtbx.load_env
        import os

        try:
            dials_regression = libtbx.env.dist_path( 'dials_regression' )
        except KeyError, e:
            print 'FAIL: dials_regression not configured'
            return

        filename = os.path.join(dials_regression,
            'centroid_test_data', 'sweep.json')

        # Get the sweep
        self.sweep = load.sweep(filename)

        # Create the threshold strategy
        threshold_image = XDSThresholdStrategy(
                kernel_size=(3, 3), gain=None,
                n_sigma_b=6, n_sigma_s=3)

        # Create the spot finder
        self.find_spots = FindSpots(threshold_image=threshold_image)

    def run(self):

        # Get all the shoeboxes
        shoeboxes = self.find_spots(self.sweep)

        # Expected number
        assert(len(shoeboxes) == 1756)

        # Test passed
        print 'OK'


if __name__ == '__main__':
    test = Test()
    test.run()
