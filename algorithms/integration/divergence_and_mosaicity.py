#
# beam_divergence_and_mosaicity.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class ComputeEsdBeamDivergence(object):
    '''Calculate the E.s.d of the beam divergence.'''

    def __init__(self, detector):
        '''Initialise the algorithm.

        Params:
            detector The detector class

        '''
        self._detector = detector

    def __call__(self, reflections):
        '''Calculate the ESD of the beam divergence.

        Params:
            reflections The list of reflections

        Returns:
            E.s.d of the beam divergence

        '''
        # Calculate the beam direction variances
        variance = self._beam_direction_variance_list(reflections)

        # Calculate and return the e.s.d of the beam divergence
        return self._compute_esd(variance)

    def _beam_direction_variance_list(self, reflections):
        '''Calculate the variance in beam direction for each spot.

        Params:
            reflections The list of reflections

        Returns:
            The list of variances

        '''
        from scitbx.array_family import flex
        from scitbx import matrix
        import numpy

        # The detector
        detector = self._detector

        # Loop through all the reflections
        variance = []
        for r in reflections:

            # Find the active pixels
            zi, yi, xi = numpy.where(r.shoebox_mask.as_numpy_array() != 0)
            index = zip(map(int, zi), map(int, yi), map(int, xi))

            # Extract the pixel values
            values = flex.double([r.shoebox[k, j, i] for k, j, i in index])

            # Get the pixel coordinates centres
            xs, xf, ys, yf, zs, zf = r.bounding_box
            xp = xi + xs + 0.5
            yp = yi + ys + 0.5
            zp = zi + zs + 0.5

            # Calculate the beam directions to each pixel
            s1 = [detector.get_pixel_lab_coord((x, y)) for x, y in zip(xp, yp)]

            # Calculate the beam vector at the centroid
            xc, yc, zc = r.centroid_position
            s1_centroid = detector.get_pixel_lab_coord((xc, yc))

            # Calculate the variance in beam vector directions
            var = self._beam_direction_variance(s1_centroid, s1, values)
            variance.append(var)

        # Return a list of variances
        return variance

    def _beam_direction_variance(self, s1_centroid, s1, values):
        '''Calculate the angles between the s1 centroid and s1 vectors for
        each pixel and then calculate the variance of the angles.

        Params:
            s1_centroid The centroid beam direction
            s1 The list of beam directions
            values The list of pixel values

        Returns:
            The variance in the beam direction

        '''
        from scitbx.array_family import flex
        from scitbx import matrix

        # Calculate angles between vectors
        s1_centroid = matrix.col(s1_centroid)
        angles = flex.double([s1_centroid.angle(matrix.col(s)) for s in s1])

        # Calculate the variance of the angles
        return flex.sum(values * (angles**2)) / flex.sum(values)

    def _compute_esd(self, variance):
        '''Calculate the beam divergence as the sum of centroid variance of the
        intensity weighted diffracted beam directions.

        Params:
            variance The variance of the beam directions

        Returns:
            The e.s.d of the beam divergence

        '''
        from math import sqrt

        # Return the beam divergence as the sum / num reflections
        return sqrt(sum(variance) / len(variance))


class FractionOfObservedIntensity(object):
    '''Calculate the fraction of observed intensity for different sigma_m.'''

    def __init__(self, reflections, sweep):
        '''Initialise the algorithm. Calculate the list of tau and zetas.

        Params:
            reflections The list of reflections
            sweep The sweep object

        '''
        self.dphi = sweep.get_scan().get_oscillation(deg=False)[1]
        self.tau, self.zeta = self._calculate_tau_and_zeta(reflections, sweep)

    def _calculate_tau_and_zeta(self, reflections, sweep):
        '''Calculate the list of tau and zeta needed for the calculation.

        Params:
            reflections The list of reflections
            sweep The sweep object.

        Returns:
            (list of tau, list of zeta)

        '''
        from scitbx.array_family import flex
        from dials.algorithms.integration import zeta_factor

        # Calculate the list of frames and z coords
        bbox = [r.bounding_box for r in reflections]
        frames = [range(b[4], b[5]) for b in bbox]
        phi = [r.rotation_angle for r in reflections]

        # Calculate the zeta list
        s0 = sweep.get_beam().get_s0()
        m2 = sweep.get_goniometer().get_rotation_axis()
        zeta = [zeta_factor(s0, r.beam_vector, m2) for r in reflections]

        # Calculate the list of tau values
        tau = []
        zeta2 = []
        scan = sweep.get_scan()
        for rf, p, z in zip(frames, phi, zeta):
            for f in rf:
                phi0 = scan.get_angle_from_array_index(int(f), deg=False)
                phi1 = scan.get_angle_from_array_index(int(f)+1, deg=False)
                tau.append((phi1 + phi0) / 2.0 - p)
                zeta2.append(z)

        # Return the list of tau and zeta
        assert(len(zeta2) == len(tau))
        return flex.double(tau), flex.double(zeta2)

    def __call__(self, sigma_m):
        '''Calculate the fraction of observed intensity for each observation.

        Params:
            sigma_m The mosaicity

        Returns:
            A list of fractions of length n

        '''
        from math import sqrt, erf, exp
        from scitbx.array_family import flex
        import numpy

        # Tiny value
        TINY = 1e-10

        # Oscillation range / 2
        dphi2 = self.dphi / 2

        # Calculate the denominator to the fraction
        den =  sqrt(2) * exp(sigma_m) / flex.abs(self.zeta)

        # Calculate the two components to the fraction
        a = flex.double([erf(x) for x in (self.tau + dphi2) / den])
        b = flex.double([erf(x) for x in (self.tau - dphi2) / den])

        # Calculate the fraction of observed reflection intensity
        R = (a - b) / 2.0

        # Set any points <= 0 to 1e-10 (otherwise will get a floating
        # point error in log calculation below).
        bad_index = numpy.where(R.as_numpy_array() < TINY)[0]
        for i in bad_index:
            R[int(i)] = TINY

        # Return the logarithm of r
        return flex.log(R)


class MaximumLikelihoodEstimator(object):
    '''Estimate E.s.d reflecting range by maximum likelihood estimation.'''

    def __init__(self, reflections, sweep):
        '''Initialise the optmization.'''
        from scitbx import simplex
        from scitbx.array_family import flex
        from math import pi
        import random

        # Initialise the function used in likelihood estimation.
        self._R = FractionOfObservedIntensity(reflections, sweep)

        # Set the starting values to try
        start = random.random() * pi / 180
        stop = random.random() * pi / 180
        starting_simplex = [flex.double([start]), flex.double([stop])]

        # Initialise the optimizer
        self._optimizer = simplex.simplex_opt(1, matrix=starting_simplex,
            evaluator=self, tolerance=1e-7)

    def target(self, sigma):
        '''Return the target function evaluated at this valid.

        Params:
            sigma The parameters

        Returns:
            The value of the target function at the given value.

        '''
        return -self.likelihood(sigma[0])

    def likelihood(self, sigma_m):
        '''Return the likelihood of the current sigma

        Params:
            sigma_m the estimated mosaicity

        Returns:
            The likelihood of this value.

        '''
        from scitbx.array_family import flex
        return flex.sum(self._R(sigma_m))

    def __call__(self):
        '''Perform maximum likelihood estimation of sigma_m

        Returns:
            The value of sigma_m

        '''
        from math import exp
        return exp(self._optimizer.get_solution()[0])


class ComputeEsdReflectingRange(object):
    '''Calculate the E.s.d of the reflecting range (mosaicity).'''

    def __init__(self, sweep):
        '''Initialise the algorithm with the scan.

        Params:
            scan The scan object

        '''
        self._sweep = sweep

    def __call__(self, reflections):
        '''Calculate the value for the reflecting range (mosaicity).

        Params:
            reflections The list of reflections

        Returns:
            The calculated value of sigma_m

        '''
        maximum_likelihood = MaximumLikelihoodEstimator(
            reflections, self._sweep)
        return maximum_likelihood()



class FractionOfObservedIntensityPerFrame(object):
    '''Calculate the fraction of observed intensity for different sigma_m.'''

    def __init__(self, reflections, sweep, index):
        '''Initialise the algorithm. Calculate the list of tau and zetas.

        Params:
            reflections The list of reflections
            sweep The sweep object

        '''
        self.dphi = sweep.get_scan().get_oscillation(deg=False)[1]
        self.tau, self.zeta = self._calculate_tau_and_zeta(reflections, sweep, index)

    def _calculate_tau_and_zeta(self, reflections, sweep, index):
        '''Calculate the list of tau and zeta needed for the calculation.

        Params:
            reflections The list of reflections
            sweep The sweep object.

        Returns:
            (list of tau, list of zeta)

        '''
        from scitbx.array_family import flex
        from dials.algorithms.integration import zeta_factor

        # Calculate the list of frames and z coords
        bbox = [r.bounding_box for r in reflections]
        frames = [range(b[4], b[5]) for b in bbox]
        phi = [r.rotation_angle for r in reflections]

        # Calculate the zeta list
        s0 = sweep.get_beam().get_s0()
        m2 = sweep.get_goniometer().get_rotation_axis()
        zeta = [zeta_factor(s0, r.beam_vector, m2) for r in reflections]

        # Calculate the list of tau values
        tau = []
        zeta2 = []
        scan = sweep.get_scan()
        for rf, p, z in zip(frames, phi, zeta):
            for f in rf:
                if f == index:#index - 1 <= f <= index + 1:
                    phi0 = scan.get_angle_from_array_index(int(f), deg=False)
                    phi1 = scan.get_angle_from_array_index(int(f)+1, deg=False)
                    tau.append((phi1 + phi0) / 2.0 - p)
                    zeta2.append(z)

#        print len(zeta2)

        # Return the list of tau and zeta
        assert(len(zeta2) == len(tau))
        return flex.double(tau), flex.double(zeta2)

    def __call__(self, sigma_m):
        '''Calculate the fraction of observed intensity for each observation.

        Params:
            sigma_m The mosaicity

        Returns:
            A list of fractions of length n

        '''
        from math import sqrt, erf, exp
        from scitbx.array_family import flex
        import numpy

        # Tiny value
        TINY = 1e-10

        # Ensure value for sigma_m is valid
        #if sigma_m < TINY:
        #    raise ValueError('sigma_m must be > 0')

        # Oscillation range / 2
        dphi2 = self.dphi / 2

        # Calculate the denominator to the fraction
        den =  sqrt(2) * exp(sigma_m) / flex.abs(self.zeta)

        # Calculate the two components to the fraction
        a = flex.double([erf(x) for x in (self.tau + dphi2) / den])
        b = flex.double([erf(x) for x in (self.tau - dphi2) / den])

        # Calculate the fraction of observed reflection intensity
        R = (a - b) / 2.0

        # Set any points <= 0 to 1e-10 (otherwise will get a floating
        # point error in log calculation below).
        bad_index = numpy.where(R.as_numpy_array() < TINY)[0]
        for i in bad_index:
            R[int(i)] = TINY

        # Return the logarithm of r
        return flex.log(R)


class MaximumLikelihoodEstimatorPerFrame(object):
    '''Estimate E.s.d reflecting range by maximum likelihood estimation.'''

    def __init__(self, reflections, sweep, index):
        '''Initialise the optmization.'''
        from scitbx import simplex
        from scitbx.array_family import flex
        from math import pi
        import random

        # Initialise the function used in likelihood estimation.
        self._R = FractionOfObservedIntensityPerFrame(reflections, sweep, index)

        # Set the starting values to try
        start = 0.1 * random.random() * pi / 180
        stop = 0.1 * random.random() * pi / 180
        starting_simplex = [flex.double([start]), flex.double([stop])]

        # Initialise the optimizer
        self._optimizer = simplex.simplex_opt(1, matrix=starting_simplex,
            evaluator=self, tolerance=1e-7)

    def target(self, sigma):
        '''Return the target function evaluated at this valid.

        Params:
            sigma The parameters

        Returns:
            The value of the target function at the given value.

        '''
        return -self.likelihood(sigma[0])

    def likelihood(self, sigma_m):
        '''Return the likelihood of the current sigma

        Params:
            sigma_m the estimated mosaicity

        Returns:
            The likelihood of this value.

        '''
        from scitbx.array_family import flex
        return flex.sum(self._R(sigma_m))

    def __call__(self):
        '''Perform maximum likelihood estimation of sigma_m

        Returns:
            The value of sigma_m

        '''
        return self._optimizer.get_solution()[0]


class ComputeEsdReflectingRangePerFrame(object):
    '''Calculate the E.s.d of the reflecting range (mosaicity).'''

    def __init__(self, sweep):
        '''Initialise the algorithm with the scan.

        Params:
            scan The scan object

        '''
        self._sweep = sweep

    def __call__(self, reflections):
        '''Calculate the value for the reflecting range (mosaicity).

        Params:
            reflections The list of reflections

        Returns:
            The calculated value of sigma_m

        '''
        sigma_m = list()
        for index in self._sweep.indices():
            try:
                maximum_likelihood = MaximumLikelihoodEstimatorPerFrame(
                    reflections, self._sweep, index)
                sigma_m.append(maximum_likelihood())
            except(ValueError):
                sigma_m.append(0)
            print sigma_m[-1]

        return sigma_m


class BeamDivergenceAndMosaicity(object):
    '''An algorithm to calculate the beam divergence and mosaicity params.'''

    def __init__(self, sweep):
        '''Initialise the algorithm.

        Params:
            sweep The sweep object
            max_separation Max pixel dist between predicted and observed spot

        '''
        # Setup the algorithm objects
        self._compute_sigma_d = ComputeEsdBeamDivergence(sweep.get_detector())
        self._compute_sigma_m = ComputeEsdReflectingRange(sweep)

        # Set the internal reflection list to None
        self._data = None

    def reflections(self):
        '''Return the updated reflections.'''
        return self._data

    def __call__(self, reflections):
        '''Calculate the divegence/mosaicity parameters.

        First match the observed reflections to the predicted reflections
        by finding the nearest neighbour based on the observed centroid
        and predicted bragg maximum. Then calculate the standard deviation
        of the beam divergence followed by the standard deviation of the
        mosaicity.

        The updated list of reflections can be accessed via the
        self.reflections member function.

        Params:
            reflections The list of reflections

        Returns:
            sigma_d, sigma_m

        '''
        from math import pi
        from dials.util.command_line import Command

        # Setup the output
        Command.indent = 4
        print '\nCalculating e.s.d of beam divergence and mosaicity...'

        # Save the reflections
        self._data = reflections

        # Calculate the standard deviation of the beam divergence
        Command.start('Calculating e.s.d of the beam divergence')
        sigma_d = self._calculate_esd_beam_divergence(self._data)
        Command.end('Calculated e.s.d of the beam divergence = {0:.3f} deg' \
            .format(sigma_d * 180 / pi))

        # Calculate the standard deviation of the reflecting range (mosaicity)
        Command.start('Calculating the mosaicity')
        sigma_m = self._calculate_esd_reflecting_range(self._data)
        #from matplotlib import pylab
        #import numpy
        #pylab.plot(numpy.array(sigma_m) * 180 / pi)
        #pylab.show()
        #print 1/0
        Command.end('Calculated mosaicity = {0:.3f} deg'\
            .format(sigma_m * 180 / pi))

        # Return the parameters
        return sigma_d, sigma_m

    def _calculate_esd_beam_divergence(self, reflections):
        return self._compute_sigma_d(reflections)

    def _calculate_esd_reflecting_range(self, reflections):
        return self._compute_sigma_m(reflections)
