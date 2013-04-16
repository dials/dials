from __future__ import division

class SpotFinder(object):
    '''A class to perform spot finding operations on a sweep of images.'''

    def __init__(self, min_spot_size=6):
        '''Initialise the algorithm with some parameters.'''
        self._min_spot_size = min_spot_size

    def __call__(self, sweep):
        '''The main function of the spot finder. Select the pixels from
        the sweep and then group the pixels into spots. Return the data
        in the form of a reflection list.

        Params:
            sweep The sweep object

        Returns:
            The reflection list

        '''
        # Extract the image pixels from the sweep
        print 'Extracting pixels from sweep'
        coords, intensity = self._extract_pixels(sweep)
        print 'Extracted {0} strong pixels.'.format(len(coords))

        # Label the pixels and group into spots
        print 'Labelling the pixels'
        labels = self._label_pixels(coords, sweep)
        print 'Formed {0} spots'.format(max(labels) + 1)

        # Filter spots that are too small
        print 'Filtering spots by size'
        spots = self._filter_spots(labels)
        print 'Filtered {0} spots'.format(len(spots))

        # Calculate the bounding box for each spot
        print 'Calculating spot bounding box.'
        bbox = self._calculate_bbox(coords, spots, sweep)
        print 'Calculated bounding box for {0} spots.'.format(len(bbox))

        # Calculate the spot centroids
        print 'Calculating centroids.'
        cpos, cvar = self._calculate_centroids(coords, intensity, spots)

        # Create a reflection list and return
        print 'Creating reflection list.'
        return self._create_reflection_list(coords, intensity, spots, bbox,
            cpos, cvar)

    def _extract_pixels(self, sweep):
        '''Extract the pixels from the sweep

        Params:
            sweep The sweep object

        Returns:
            The list of selected pixels

        '''
        from dials.util.command_line import ProgressBar
        from scitbx.array_family import flex
        from dials.algorithms.peak_finding import flex_vec3_int

        # Initialise the pixel arrays
        coords = flex_vec3_int()
        intensity = flex.int()

        # Get the start index and trusted range from the sweep
        start = sweep.get_array_range()[0]
        trusted_range = sweep.get_detector().get_trusted_range()

        # Loop through all the images in the sweep and extract the pixels
        # from each of the images
        progress = ProgressBar()
        for frame, image in enumerate(sweep):
            c, i = self._extract_image_pixels(image, frame + start,
                trusted_range)
            coords.extend(c)
            intensity.extend(i)
            progress.update(100.0 * float(frame + 1) / len(sweep))

        progress.finished()

        # Reuturn the pixel values
        return coords, intensity

    def _extract_image_pixels(self, flex_image, frame, trusted_range):
        '''Extract pixels from a single image

        Params:
            image The image from which to extract the pixels
            frame The frame number of the current image
            trusted_range The trusted range of pixel values

        Returns:
            The list of selected pixels

        '''
        from scitbx.array_family import flex
        from dials.algorithms.peak_finding import flex_vec3_int
        import numpy

        image = flex_image.as_numpy_array()
        height, width = image.shape

        # Calculate the threshold
        threshold = self._calculate_threshold(image, trusted_range)

        # Extract the pixel indices
        index = numpy.where(image >= threshold)[0]

        # Create the array of pixel coordinates
        z = [frame] * len(index)
        y = map(int, index // width)
        x = map(int, index  % width)
        coords = flex_vec3_int(zip(x, y, z))

        # Get the array of pixel intensities
        intensity = flex.int(image[index])

        # Return the pixel values
        return coords, intensity

    def _calculate_threshold(self, image, trusted_range):
        '''Calculate the threshold for this image.

        Params:
            image The image to process
            trusted_range The trusted range of pixel values

        Returns:
            The threshold value

        '''
        from scipy.ndimage.measurements import histogram
        from thresholding import maximum_deviation
        import numpy

        # Cap pixels to within trusted range
        image.shape = -1
        ind = numpy.where(image < trusted_range[0])
        image[ind] = trusted_range[0]
        ind = numpy.where(image > trusted_range[1])
        image[ind] = trusted_range[1]

        # Histogram the pixels
        histo = histogram(image, trusted_range[0],
            trusted_range[1], trusted_range[1])
        histo = histo / numpy.sum(histo)

        # Calculate the threshold and add to list
        return maximum_deviation(histo)

    def _label_pixels(self, pixels, sweep):
        '''Do a connected component labelling of the pixels to get
        groups of spots.

        Params:
            pixels The pixel coordinates
            sweep The sweep object

        Returns:
            The pixel-spot mapping.

        '''
        from dials.algorithms.peak_finding import LabelPixels

        # Get the grid size needed by the labelling algorithm
        image_size = sweep.get_image_size()
        scan_size = sweep.get_array_range()[1]
        grid_size = (image_size[0], image_size[1], scan_size)

        # Label the pixels
        label_pixels = LabelPixels(grid_size)
        labels = label_pixels(pixels)

        # Return the labels
        return labels

    def _filter_spots(self, labels):
        '''Filter spots that are too small.

        Params:
            spots The input spots

        Returns:
            The filtered spots

        '''
        # Create a list of lists of pixels representing the spots
        spots = [list() for i in range(max(labels) + 1)]
        for i, l in enumerate(labels):
            spots[l].append(i)

        # Filter by spot size
        return [s for s in spots if len(s) >= self._min_spot_size]

    def _calculate_bbox(self, coords, spots, sweep):
        '''Calculate the bounding boxes for each spot.

        Params:
            coords The pixel coordinates
            spots The pixel-spot mapping
            sweep The sweep object

        Returns:
            The bounding boxes for each spot

        '''
        # Get the image dimensions
        height, width = sweep.get_image_size()
        length = sweep.get_array_range()[1]

        # Loop through all spots
        bbox = []
        for s in spots:

            # Loop through all pixels in spot and find the min/max
            # x/y/z coordinates
            b = [width, -1, height, -1, length, -1]
            for i in s:
                x, y, z = coords[i]
                if x < b[0]: b[0] = x
                if x > b[1]: b[1] = x
                if y < b[2]: b[2] = y
                if y > b[3]: b[3] = y
                if z < b[4]: b[4] = z
                if z > b[5]: b[5] = z

            # Make the bounding box and check it is valid
            b[1] += 1
            b[3] += 1
            b[5] += 1
            assert(b[0] < b[1] and b[2] < b[3] and b[4] < b[5])

            # Add to list
            bbox.append(tuple(b))

        # Return the bounding boxes
        return bbox

    def _calculate_centroids(self, coords, intensity, spots):
        '''Calculate the spot centroids.

        Params:
            coords The spot coords
            intensity The spot intensities
            spots The pixel-spot mapping

        Returns:
            (centroid position, centroid variance)

        '''
        from scitbx.array_family import flex

        # Loop through all the spots
        centroid_pos = flex.vec3_double()
        centroid_var = flex.vec3_double()
        for s in spots:

            # Get pixel coords and values
            pixel_coords = [map(lambda x: x + 0.5, coords[i]) for i in s]
            pixel_values = flex.double([intensity[i] for i in s])
            pixel_x, pixel_y, pixel_z = zip(*pixel_coords)

            # Calculate the centroid and variance
            xc = flex.mean_and_variance(flex.double(pixel_x), pixel_values)
            yc = flex.mean_and_variance(flex.double(pixel_y), pixel_values)
            zc = flex.mean_and_variance(flex.double(pixel_z), pixel_values)

            # Add the centroid and variance
            centroid_pos.append((xc.mean(), yc.mean(), zc.mean()))
            centroid_var.append((xc.gsl_stats_wvariance(),
                                 yc.gsl_stats_wvariance(),
                                 zc.gsl_stats_wvariance()))

        # Return the centroid and variance
        return centroid_pos, centroid_var

    def _create_reflection_list(self, coords, values, spots, bbox, cpos, cvar):
        '''Create a reflection list from the spot data.

        Params:
            coords The pixel coordinates
            values The pixel values
            spots The pixel->spot mapping
            bbox The spot bounding boxes
            cpos The centroid position
            cvar The centroid variance

        Returns:
            A list of reflections

        '''
        from dials.model.data import Reflection, ReflectionList
        from dials.algorithms.integration import allocate_reflection_profiles

        # Ensure the lengths are ok
        assert(len(spots) == len(bbox))
        assert(len(spots) == len(cpos))
        assert(len(spots) == len(cvar))

        # Create the reflection list
        reflection_list = ReflectionList(len(spots))
        for i in range(len(spots)):
            r = Reflection()
            r.bounding_box = bbox[i]
            r.centroid_position = cpos[i]
            r.centroid_variance = cvar[i]
            reflection_list[i] = r

        # Allocate memory for the reflection profiles
        allocate_reflection_profiles(reflection_list,
            shoebox_default=0, shoebox_mask_default=0)

        # Set the pixel and mask values
        for i in range(len(spots)):
            r = reflection_list[i]
            xs, xf, ys, yf, zs, zf = r.bounding_box
            for s in spots[i]:
                x, y, z = coords[s]
                x = x - xs
                y = y - ys
                z = z - zs
                r.shoebox[z, y, x] = values[s]
                r.shoebox_mask[z, y, x] = 1

        # Return the reflection list
        return reflection_list
