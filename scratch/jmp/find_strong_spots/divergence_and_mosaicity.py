
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

        # Calculate the sum of s^2
        sum_variance = reduce(lambda x, y: x + y, variance)

        # Return the beam divergence as the sum / num reflections
        return sqrt(sum_variance / len(variance))


class ComputeEsdReflectingRange(object):

    def __init__(self):
        pass

    def __call__(self, reflections):
        pass


class BeamDivergenceAndMosaicity(object):

    def __init__(self, detector, max_separation=2):
        '''Initialise the algorithm.

        Params:
            max_separation Max pixel dist between predicted and observed spot

        '''
        self._compute_sigma_d = ComputeEsdBeamDivergence(detector)

        # Set the parameters
        self._max_separation = max_separation

        # Set the internal reflection list to None
        self._data = None

    def reflections(self):
        '''Return the updated reflections.'''
        return self._reflections

    def __call__(self, observed, predicted):
        '''Calculate the divegence/mosaicity parameters.

        First match the observed reflections to the predicted reflections
        by finding the nearest neighbour based on the observed centroid
        and predicted bragg maximum. Then calculate the standard deviation
        of the beam divergence followed by the standard deviation of the
        mosaicity.

        The updated list of reflections can be accessed via the
        self.reflections member function.

        Params:
            observed The observed list of reflections
            predicted The predicted list of reflections

        Returns:
            sigma_d, sigma_m

        '''
        # Map observed to predicted reflections
        print 'Matching observed and predicted reflections.'
        self._data = self._match_observed_w_predicted(observed, predicted)
        print '{0} reflections remaining.'.format(len(self._data))

        # Calculate the standard deviation of the beam divergence
        sigma_d = self._calculate_esd_beam_divergence(self._data)
        from math import pi
        print sigma_d * 180 / pi

        # Calculate the standard deviation of the reflecting range (mosaicity)
        #sigma_m = self._calculate_esd_reflecting_range()

        # Return the parameters
        #return sigma_d, sigma_m

    def _match_observed_w_predicted(self, observed, predicted):
        '''Match the observed reflections with the predicted.

        Params:
            observed The list of observed reflections.
            predicted The list of predicted reflections.

        Returns:
            The list of matched reflections

        '''
        from dials.model.data import ReflectionList

        # Find the nearest neighbours and distances
        nn, dist = self._find_nearest_neighbours(observed, predicted)

        # Filter the matches by distance
        index = self._filter_by_distance(nn, dist)

        # Copy all of the reflection data for the matched reflections
        reflections = ReflectionList()
        for i in index:
            o = observed[i]
            p = predicted[nn[i]]
            o.miller_index = p.miller_index
            o.rotation_angle = p.rotation_angle
            o.beam_vector = p.beam_vector
            o.image_coord_px = p.image_coord_px
            o.image_coord_mm = p.image_coord_mm
            o.panel_number = p.panel_number
            o.frame_number = p.frame_number
            reflections.append(o)

        # Return the list of reflections
        return reflections


    def _calculate_esd_beam_divergence(self, reflections):
        return self._compute_sigma_d(reflections)

    def _calculate_esd_reflecting_range(self):
        pass

    def _find_nearest_neighbours(self, observed, predicted):
        '''Find the nearest predicted spot to the observed spot.

        Params:
            observed The observed reflections
            predicted The predicted reflections

        Returns:
            (nearest neighbours, distance)

        '''
        from annlib_ext import AnnAdaptor
        from scitbx.array_family import flex
        from math import sqrt

        # Get the predicted coordinates
        predicted_xyz = []
        for r in predicted:
            x, y = r.image_coord_px
            z = r.frame_number
            predicted_xyz.append((x, y, z))

        # Create the KD Tree
        ann = AnnAdaptor(flex.double(predicted_xyz).as_1d(), 3)

        # Get the observed coordinates
        observed_xyz = [r.centroid_position for r in observed]

        # Query to find all the nearest neighbours
        ann.query(flex.double(observed_xyz).as_1d())

        # Return the nearest neighbours and distances
        return ann.nn, flex.sqrt(ann.distances)

    def _filter_by_distance(self, nn, dist):
        '''Filter the matches by distance.

        Params:
            nn The nearest neighbour list
            dist The distances

        Returns:
            A reduced list of nearest neighbours

        '''
        from scitbx.array_family import flex
        index = range(len(nn))
        return flex.int([i for i in index if dist[i] <= self._max_separation])
