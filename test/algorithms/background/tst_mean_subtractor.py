from __future__ import division

class Test(object):

    EPS = 1e-7

    def __init__(self):
        from dials.algorithms.background import MeanSubtractor

        # Create the subtraction algorithm
        self.subtract = MeanSubtractor()

    def run(self):

        from scitbx.array_family import flex
        from dials.model.data import Reflection

        # Create shoebox and mask data
        shoebox_d = flex.random_double(5 * 5 * 5) * 100
        shoebox = flex.int(flex.grid(5, 5, 5))
        for idx in range(len(shoebox_d)):
            shoebox[idx] = int(shoebox_d[idx])
        mask = flex.int(flex.grid(5, 5, 5), (1 << 0))
        for k in range(1, 4):
            for j in range(1, 4):
                for i in range(1, 4):
                    mask[k,j,i] = (1 << 1)

        # Create the reflection data
        reflection = Reflection()
        reflection.shoebox = shoebox.deep_copy()
        reflection.shoebox_mask = mask.deep_copy()


        # Subtract background
        self.subtract(reflection)

        # Test to see if we're correct
        pixels = flex.double(flex.select(shoebox, flags=(mask == (1 << 0))))
        background = flex.mean(pixels)
        shoebox = shoebox - int(background)
        for p1, p2 in zip(shoebox, reflection.shoebox):
            assert(abs(p1 - p2) < self.EPS)

        # Test passed
        print 'OK'

if __name__ == '__main__':
    test = Test()
    test.run()
