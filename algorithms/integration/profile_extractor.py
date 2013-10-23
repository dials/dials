#
# copy_reflection_profiles.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class ProfileExtractor(object):
    ''' A class to extract the profiles from the sweep '''

    def __init__(self, sweep, crystal, mask=None, gain=None, dark=None):
        ''' Initialise the class with the sweep etc

        Params:
            sweep The sweep to process
            crystal The crystal to process
            mask The detector mask
            gain The gain map
            dark The dark map

        '''
        from dials.algorithms import shoebox
        from scitbx.array_family import flex

        # Save the sweep
        self.sweep = sweep

        # Ensure image is a tuple
        image = sweep[0]
        if not isinstance(image, tuple):
            image = (image,)

        # Get the mask in tuple of masks form
        if mask:
            if not isinstance(mask, tuple):
                mask = (mask,)
        else:
            mask = tuple([im >= 0 for im in image])

        # Get the gain in tuple of gains form
        if gain:
            if not isinstance(gain, tuple):
                gain = (gain,)
        else:
            gain = tuple([flex.double(flex.grid(im.all()), 1) for im in image])

        # Get the dark in tuple of darks form
        if dark:
            if not isinstance(dark, tuple):
                dark = (dark,)
        else:
            dark = tuple([flex.double(flex.grid(im.all()), 0) for im in image])

        # Set the mask, gain and dark maps
        self.mask = mask
        self.gain = gain
        self.dark = dark

    def __call__(self, panels, bboxes):
        ''' Extract the profiles from the sweep

        Params:
            panels The panel numbers
            bboxes The bounding boxes

        Returns:
            The reflection list

        '''
        from dials.util.command_line import ProgressBar
        from dials.algorithms.shoebox import Extractor

        # Create the class to set all the shoebox pixels
        extractor = Extractor(panels, bboxes, self.mask, self.gain, self.dark)

        # For each image in the sweep, get the reflections predicted to have
        # been recorded on the image and copy the pixels from the image to
        # the reflection profile image. Update a progress bar as we go along
        progress = ProgressBar(title = "Extracting reflections")
        first_array_index = self.sweep.get_array_range()[0]
        for index, image in enumerate(self.sweep, start=first_array_index):

            # Ensure the image is a tuple
            if not isinstance(image, tuple):
                image = (image,)

            # Loop through all the images and add to the extractor
            for panel, im in enumerate(image):
                extractor.add_image(panel, index, im)

            # Update the progress bar
            progress.update(100*(index-first_array_index+1) / len(self.sweep))

        # Get the shoeboxes from the extractor
        shoeboxes = extractor.shoeboxes()

        # Finish the progress bar and return the profiles
        progress.finished("Extracted {0} profiles".format(len(shoeboxes)))
        return shoeboxes


class PartialProfileExtractor(object):
    ''' A class to extract the profiles from the sweep '''

    def __init__(self, sweep):
        ''' Initialise the class with the sweep etc

        Params:
            sweep The sweep to process

        '''
        self.sweep = sweep

    def __call__(self, panels, bboxes):
        ''' Extract the profiles from the sweep

        Params:
            panels The panel numbers
            bboxes The bounding boxes

        Returns:
            The reflection list

        '''
        from dials.util.command_line import ProgressBar
        from dials.algorithms.shoebox import PartialExtractor

        # Set the number of panels
        image = self.sweep[0]
        if isinstance(image, tuple):
            npanels = len(image)
        else:
            npanels = 1

        # Get the z range
        zrange = self.sweep.get_array_range()

        # Create the class to set all the shoebox pixels
        extractor = PartialExtractor(panels, bboxes, zrange, npanels)

        # For each image in the sweep, get the reflections predicted to have
        # been recorded on the image and copy the pixels from the image to
        # the reflection profile image. Update a progress bar as we go along
        progress = ProgressBar(title = "Extracting frames %d -> %d" % zrange)
        first_array_index = self.sweep.get_array_range()[0]
        for index, image in enumerate(self.sweep, start=first_array_index):

            # Ensure the image is a tuple
            if not isinstance(image, tuple):
                image = (image,)

            # Loop through all the images and add to the extractor
            for panel, im in enumerate(image):
                extractor.add_image(panel, index, im)

            # Update the progress bar
            progress.update(100*(index-first_array_index+1) / len(self.sweep))

        # Get the shoeboxes from the extractor
        shoeboxes = extractor.shoeboxes()
        indices = extractor.shoebox_indices()

        # Finish the progress bar and return the profiles
        progress.finished("Extracted %d profiles from frames %d -> %d"
            % ((len(shoeboxes),) + zrange))

        # Return the indices and shoeboxes
        return (indices, shoeboxes)
