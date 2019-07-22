#
# finder.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import logging
import math
import os

from dials.util import Sorry
import libtbx

logger = logging.getLogger(__name__)

_no_multiprocessing_on_windows = (
    "\n"
    + "*" * 80
    + "\n"
    + "Multiprocessing is not available on windows. Setting nproc = 1, njobs = 1"
    + "\n"
    + "*" * 80
    + "\n"
)


class Result(object):
    """
    A class to hold the result from spot finding on an image.

    When doing multi processing, we can process the result of
    each thread as it comes in instead of waiting for all results.
    The purpose of this class is to allow us to set the pixel list
    to None after each image to lower memory usage.
    """

    def __init__(self, pixel_list):
        """
        Set the pixel list
        """
        self.pixel_list = pixel_list


class ExtractPixelsFromImage(object):
    """
    A class to extract pixels from a single image
    """

    def __init__(
        self,
        imageset,
        threshold_function,
        mask,
        region_of_interest,
        max_strong_pixel_fraction,
        compute_mean_background,
    ):
        """
        Initialise the class

        :param imageset: The imageset to extract from
        :param threshold_function: The function to threshold with
        :param mask: The image mask
        :param region_of_interest: A region of interest to process
        :param max_strong_pixel_fraction: The maximum fraction of pixels allowed
        """
        self.threshold_function = threshold_function
        self.imageset = imageset
        self.mask = mask
        self.region_of_interest = region_of_interest
        self.max_strong_pixel_fraction = max_strong_pixel_fraction
        self.compute_mean_background = compute_mean_background
        if self.mask is not None:
            detector = self.imageset.get_detector()
            assert len(self.mask) == len(detector)
        self.first = True

    def __call__(self, index):
        """
        Extract strong pixels from an image

        :param index: The index of the image
        """
        from dials.model.data import PixelList
        from dxtbx.imageset import ImageSweep
        from dials.array_family import flex

        # Parallel reading of HDF5 from the same handle is not allowed. Python
        # multiprocessing is a bit messed up and used fork on linux so need to
        # close and reopen file.
        if self.first:
            if self.imageset.reader().is_single_file_reader():
                self.imageset.reader().nullify_format_instance()
            self.first = False

        # Get the frame number
        if isinstance(self.imageset, ImageSweep):
            frame = self.imageset.get_array_range()[0] + index
        else:
            ind = self.imageset.indices()
            if len(ind) > 1:
                assert all(i1 + 1 == i2 for i1, i2 in zip(ind[0:-1], ind[1:-1]))
            frame = ind[index]

        # Create the list of pixel lists
        pixel_list = []

        # Get the image and mask
        image = self.imageset.get_corrected_data(index)
        mask = self.imageset.get_mask(index)

        # Set the mask
        if self.mask is not None:
            assert len(self.mask) == len(mask)
            mask = tuple(m1 & m2 for m1, m2 in zip(mask, self.mask))

        logger.debug(
            "Number of masked pixels for image %i: %i"
            % (index, sum(m.count(False) for m in mask))
        )

        # Add the images to the pixel lists
        num_strong = 0
        average_background = 0
        for im, mk in zip(image, mask):
            if self.region_of_interest is not None:
                x0, x1, y0, y1 = self.region_of_interest
                height, width = im.all()
                assert x0 < x1, "x0 < x1"
                assert y0 < y1, "y0 < y1"
                assert x0 >= 0, "x0 >= 0"
                assert y0 >= 0, "y0 >= 0"
                assert x1 <= width, "x1 <= width"
                assert y1 <= height, "y1 <= height"
                im_roi = im[y0:y1, x0:x1]
                mk_roi = mk[y0:y1, x0:x1]
                tm_roi = self.threshold_function.compute_threshold(im_roi, mk_roi)
                threshold_mask = flex.bool(im.accessor(), False)
                threshold_mask[y0:y1, x0:x1] = tm_roi
            else:
                threshold_mask = self.threshold_function.compute_threshold(im, mk)

            # Add the pixel list
            plist = PixelList(frame, im, threshold_mask)
            pixel_list.append(plist)

            # Get average background
            if self.compute_mean_background:
                background = im.as_1d().select((mk & ~threshold_mask).as_1d())
                average_background += flex.mean(background)

            # Add to the spot count
            num_strong += len(plist)

        # Make average background
        average_background /= len(image)

        # Check total number of strong pixels
        if self.max_strong_pixel_fraction < 1:
            num_image = 0
            for im in image:
                num_image += len(im)
            max_strong = int(math.ceil(self.max_strong_pixel_fraction * num_image))
            if num_strong > max_strong:
                raise RuntimeError(
                    """
          The number of strong pixels found (%d) is greater than the
          maximum allowed (%d). Try changing spot finding parameters
        """
                    % (num_strong, max_strong)
                )

        # Print some info
        if self.compute_mean_background:
            logger.info(
                "Found %d strong pixels on image %d with average background %f"
                % (num_strong, frame + 1, average_background)
            )
        else:
            logger.info("Found %d strong pixels on image %d" % (num_strong, frame + 1))

        # Return the result
        return Result(pixel_list)


class ExtractPixelsFromImage2DNoShoeboxes(ExtractPixelsFromImage):
    """
    A class to extract pixels from a single image
    """

    def __init__(
        self,
        imageset,
        threshold_function,
        mask,
        region_of_interest,
        max_strong_pixel_fraction,
        compute_mean_background,
        min_spot_size,
        max_spot_size,
        filter_spots,
    ):
        """
        Initialise the class

        :param imageset: The imageset to extract from
        :param threshold_function: The function to threshold with
        :param mask: The image mask
        :param region_of_interest: A region of interest to process
        :param max_strong_pixel_fraction: The maximum fraction of pixels allowed
        """
        super(ExtractPixelsFromImage2DNoShoeboxes, self).__init__(
            imageset,
            threshold_function,
            mask,
            region_of_interest,
            max_strong_pixel_fraction,
            compute_mean_background,
        )

        # Save some stuff
        self.min_spot_size = min_spot_size
        self.max_spot_size = max_spot_size
        self.filter_spots = filter_spots

    def __call__(self, index):
        """
        Extract strong pixels from an image

        :param index: The index of the image
        """
        from dials.model.data import PixelListLabeller

        # Initialise the pixel labeller
        num_panels = len(self.imageset.get_detector())
        pixel_labeller = [PixelListLabeller() for p in range(num_panels)]

        # Call the super function
        result = super(ExtractPixelsFromImage2DNoShoeboxes, self).__call__(index)

        # Add pixel lists to the labeller
        assert len(pixel_labeller) == len(result.pixel_list), "Inconsistent size"
        for plabeller, plist in zip(pixel_labeller, result.pixel_list):
            plabeller.add(plist)

        # Create shoeboxes from pixel list
        converter = PixelListToReflectionTable(
            self.min_spot_size, self.max_spot_size, self.filter_spots, False
        )
        reflections, _ = converter(self.imageset, pixel_labeller)

        # Delete the shoeboxes
        del reflections["shoeboxes"]

        # Return the reflections
        return [reflections]


class ExtractSpotsParallelTask(object):
    """
    Execute the spot finder task in parallel

    We need this external class so that we can pickle it for cluster jobs
    """

    def __init__(self, function):
        """
        Initialise with the function to call
        """
        self.function = function

    def __call__(self, task):
        """
        Call the function with th task and save the IO
        """
        from dials.util import log

        log.config_simple_cached()
        result = self.function(task)
        handlers = logging.getLogger("dials").handlers
        assert len(handlers) == 1, "Invalid number of logging handlers"
        return result, handlers[0].messages()


class PixelListToShoeboxes(object):
    """
    A helper class to convert pixel list to shoeboxes
    """

    def __init__(self, min_spot_size, max_spot_size, write_hot_pixel_mask):
        """
        Initialize
        """
        self.min_spot_size = min_spot_size
        self.max_spot_size = max_spot_size
        self.write_hot_pixel_mask = write_hot_pixel_mask

    def __call__(self, imageset, pixel_labeller):
        """
        Convert the pixel list to shoeboxes
        """
        from dxtbx.imageset import ImageSweep
        from dials.array_family import flex

        # Extract the pixel lists into a list of reflections
        shoeboxes = flex.shoebox()
        spotsizes = flex.size_t()
        hotpixels = tuple(flex.size_t() for i in range(len(imageset.get_detector())))
        if isinstance(imageset, ImageSweep):
            twod = False
        else:
            twod = True
        for i, (p, hp) in enumerate(zip(pixel_labeller, hotpixels)):
            if p.num_pixels() > 0:
                creator = flex.PixelListShoeboxCreator(
                    p,
                    i,  # panel
                    0,  # zrange
                    twod,  # twod
                    self.min_spot_size,  # min_pixels
                    self.max_spot_size,  # max_pixels
                    self.write_hot_pixel_mask,
                )
                shoeboxes.extend(creator.result())
                spotsizes.extend(creator.spot_size())
                hp.extend(creator.hot_pixels())
        logger.info("")
        logger.info("Extracted {} spots".format(len(shoeboxes)))

        # Get the unallocated spots and print some info
        selection = shoeboxes.is_allocated()
        shoeboxes = shoeboxes.select(selection)
        ntoosmall = (spotsizes < self.min_spot_size).count(True)
        ntoolarge = (spotsizes > self.max_spot_size).count(True)
        assert ntoosmall + ntoolarge == selection.count(False)
        logger.info(
            "Removed %d spots with size < %d pixels" % (ntoosmall, self.min_spot_size)
        )
        logger.info(
            "Removed %d spots with size > %d pixels" % (ntoolarge, self.max_spot_size)
        )

        # Return the shoeboxes
        return shoeboxes, hotpixels


class ShoeboxesToReflectionTable(object):
    """
    A class to filter shoeboxes and create reflection table
    """

    def __init__(self, filter_spots):
        """
        Initialise the reflection table creator
        """
        self.filter_spots = filter_spots

    def __call__(self, imageset, shoeboxes):
        """
        Filter shoeboxes and create reflection table
        """
        from dials.array_family import flex

        # Calculate the spot centroids
        centroid = shoeboxes.centroid_valid()
        logger.info("Calculated {} spot centroids".format(len(shoeboxes)))

        # Calculate the spot intensities
        intensity = shoeboxes.summed_intensity()
        logger.info("Calculated {} spot intensities".format(len(shoeboxes)))

        # Create the observations
        observed = flex.observation(shoeboxes.panels(), centroid, intensity)

        # Filter the reflections and select only the desired spots
        flags = self.filter_spots(
            None, sweep=imageset, observations=observed, shoeboxes=shoeboxes
        )
        observed = observed.select(flags)
        shoeboxes = shoeboxes.select(flags)

        # Return as a reflection list
        return flex.reflection_table(observed, shoeboxes)


class PixelListToReflectionTable(object):
    """
    Helper class to convert the pixel list to reflection table
    """

    def __init__(
        self, min_spot_size, max_spot_size, filter_spots, write_hot_pixel_mask
    ):
        """
        Initialise the converter
        """

        # Setup the pixel list to shoebox converter
        self.pixel_list_to_shoeboxes = PixelListToShoeboxes(
            min_spot_size, max_spot_size, write_hot_pixel_mask
        )

        # Setup the reflection table converter
        self.shoeboxes_to_reflection_table = ShoeboxesToReflectionTable(filter_spots)

    def __call__(self, imageset, pixel_labeller):
        """
        Convert to reflection table
        """
        shoeboxes, hot_pixels = self.pixel_list_to_shoeboxes(imageset, pixel_labeller)

        return self.shoeboxes_to_reflection_table(imageset, shoeboxes), hot_pixels


class ExtractSpots(object):
    """
    Class to find spots in an image and extract them into shoeboxes.
    """

    def __init__(
        self,
        threshold_function=None,
        mask=None,
        region_of_interest=None,
        max_strong_pixel_fraction=0.1,
        compute_mean_background=False,
        mp_method=None,
        mp_nproc=1,
        mp_njobs=1,
        mp_chunksize=1,
        min_spot_size=1,
        max_spot_size=20,
        filter_spots=None,
        no_shoeboxes_2d=False,
        min_chunksize=50,
        write_hot_pixel_mask=False,
    ):
        """
        Initialise the class with the strategy

        :param threshold_function: The image thresholding strategy
        :param mask: The mask to use
        :param mp_method: The multi processing method
        :param nproc: The number of processors
        :param max_strong_pixel_fraction: The maximum number of strong pixels
        """
        # Set the required strategies
        self.threshold_function = threshold_function
        self.mask = mask
        self.mp_method = mp_method
        self.mp_chunksize = mp_chunksize
        self.mp_nproc = mp_nproc
        self.mp_njobs = mp_njobs
        self.max_strong_pixel_fraction = max_strong_pixel_fraction
        self.compute_mean_background = compute_mean_background
        self.region_of_interest = region_of_interest
        self.min_spot_size = min_spot_size
        self.max_spot_size = max_spot_size
        self.filter_spots = filter_spots
        self.no_shoeboxes_2d = no_shoeboxes_2d
        self.min_chunksize = min_chunksize
        self.write_hot_pixel_mask = write_hot_pixel_mask

    def __call__(self, imageset):
        """
        Find the spots in the imageset

        :param imageset: The imageset to process
        :return: The list of spot shoeboxes
        """
        if not self.no_shoeboxes_2d:
            return self._find_spots(imageset)
        else:
            return self._find_spots_2d_no_shoeboxes(imageset)

    def _compute_chunksize(self, nimg, nproc, min_chunksize):
        """
        Compute the chunk size for a given number of images and processes
        """
        chunksize = int(math.ceil(nimg / nproc))
        remainder = nimg % (chunksize * nproc)
        test_chunksize = chunksize - 1
        while test_chunksize >= min_chunksize:
            test_remainder = nimg % (test_chunksize * nproc)
            if test_remainder <= remainder:
                chunksize = test_chunksize
                remainder = test_remainder
            test_chunksize -= 1
        return chunksize

    def _find_spots(self, imageset):
        """
        Find the spots in the imageset

        :param imageset: The imageset to process
        :return: The list of spot shoeboxes
        """
        from dials.model.data import PixelListLabeller
        from dials.util.mp import batch_multi_node_parallel_map

        # Change the number of processors if necessary
        mp_nproc = self.mp_nproc
        mp_njobs = self.mp_njobs
        if os.name == "nt" and (mp_nproc > 1 or mp_njobs > 1):
            logger.warning(_no_multiprocessing_on_windows)
            mp_nproc = 1
            mp_njobs = 1
        if mp_nproc * mp_njobs > len(imageset):
            mp_nproc = min(mp_nproc, len(imageset))
            mp_njobs = int(math.ceil(len(imageset) / mp_nproc))

        mp_method = self.mp_method
        mp_chunksize = self.mp_chunksize

        if mp_chunksize == libtbx.Auto:
            mp_chunksize = self._compute_chunksize(
                len(imageset), mp_njobs * mp_nproc, self.min_chunksize
            )
            logger.info("Setting chunksize=%i" % mp_chunksize)

        len_by_nproc = int(math.floor(len(imageset) / (mp_njobs * mp_nproc)))
        if mp_chunksize > len_by_nproc:
            mp_chunksize = len_by_nproc
        if mp_chunksize == 0:
            mp_chunksize = 1
        assert mp_nproc > 0, "Invalid number of processors"
        assert mp_njobs > 0, "Invalid number of jobs"
        assert mp_njobs == 1 or mp_method is not None, "Invalid cluster method"
        assert mp_chunksize > 0, "Invalid chunk size"

        # The extract pixels function
        function = ExtractPixelsFromImage(
            imageset=imageset,
            threshold_function=self.threshold_function,
            mask=self.mask,
            max_strong_pixel_fraction=self.max_strong_pixel_fraction,
            compute_mean_background=self.compute_mean_background,
            region_of_interest=self.region_of_interest,
        )

        # The indices to iterate over
        indices = list(range(len(imageset)))

        # Initialise the pixel labeller
        num_panels = len(imageset.get_detector())
        pixel_labeller = [PixelListLabeller() for p in range(num_panels)]

        # Do the processing
        logger.info("Extracting strong pixels from images")
        if mp_njobs > 1:
            logger.info(
                " Using %s with %d parallel job(s) and %d processes per node\n"
                % (mp_method, mp_njobs, mp_nproc)
            )
        else:
            logger.info(" Using multiprocessing with %d parallel job(s)\n" % (mp_nproc))
        if mp_nproc > 1 or mp_njobs > 1:

            def process_output(result):
                for message in result[1]:
                    logger.log(message.levelno, message.msg)
                assert len(pixel_labeller) == len(
                    result[0].pixel_list
                ), "Inconsistent size"
                for plabeller, plist in zip(pixel_labeller, result[0].pixel_list):
                    plabeller.add(plist)
                result[0].pixel_list = None

            batch_multi_node_parallel_map(
                func=ExtractSpotsParallelTask(function),
                iterable=indices,
                nproc=mp_nproc,
                njobs=mp_njobs,
                cluster_method=mp_method,
                chunksize=mp_chunksize,
                callback=process_output,
            )
        else:
            for task in indices:
                result = function(task)
                assert len(pixel_labeller) == len(
                    result.pixel_list
                ), "Inconsistent size"
                for plabeller, plist in zip(pixel_labeller, result.pixel_list):
                    plabeller.add(plist)
                    result.pixel_list = None

        # Create shoeboxes from pixel list
        converter = PixelListToReflectionTable(
            self.min_spot_size,
            self.max_spot_size,
            self.filter_spots,
            self.write_hot_pixel_mask,
        )
        return converter(imageset, pixel_labeller)

    def _find_spots_2d_no_shoeboxes(self, imageset):
        """
        Find the spots in the imageset

        :param imageset: The imageset to process
        :return: The list of spot shoeboxes
        """
        from dials.array_family import flex
        from dials.util.mp import batch_multi_node_parallel_map

        # Change the number of processors if necessary
        mp_nproc = self.mp_nproc
        mp_njobs = self.mp_njobs
        if os.name == "nt" and (mp_nproc > 1 or mp_njobs > 1):
            logger.warning(_no_multiprocessing_on_windows)
            mp_nproc = 1
            mp_njobs = 1
        if mp_nproc * mp_njobs > len(imageset):
            mp_nproc = min(mp_nproc, len(imageset))
            mp_njobs = int(math.ceil(len(imageset) / mp_nproc))

        mp_method = self.mp_method
        mp_chunksize = self.mp_chunksize

        if mp_chunksize == libtbx.Auto:
            mp_chunksize = self._compute_chunksize(
                len(imageset), mp_njobs * mp_nproc, self.min_chunksize
            )
            logger.info("Setting chunksize=%i" % mp_chunksize)

        len_by_nproc = int(math.floor(len(imageset) / (mp_njobs * mp_nproc)))
        if mp_chunksize > len_by_nproc:
            mp_chunksize = len_by_nproc
        assert mp_nproc > 0, "Invalid number of processors"
        assert mp_njobs > 0, "Invalid number of jobs"
        assert mp_njobs == 1 or mp_method is not None, "Invalid cluster method"
        assert mp_chunksize > 0, "Invalid chunk size"

        # The extract pixels function
        function = ExtractPixelsFromImage2DNoShoeboxes(
            imageset=imageset,
            threshold_function=self.threshold_function,
            mask=self.mask,
            max_strong_pixel_fraction=self.max_strong_pixel_fraction,
            compute_mean_background=self.compute_mean_background,
            region_of_interest=self.region_of_interest,
            min_spot_size=self.min_spot_size,
            max_spot_size=self.max_spot_size,
            filter_spots=self.filter_spots,
        )

        # The indices to iterate over
        indices = list(range(len(imageset)))

        # The resulting reflections
        reflections = flex.reflection_table()

        # Do the processing
        logger.info("Extracting strong spots from images")
        if mp_njobs > 1:
            logger.info(
                " Using %s with %d parallel job(s) and %d processes per node\n"
                % (mp_method, mp_njobs, mp_nproc)
            )
        else:
            logger.info(" Using multiprocessing with %d parallel job(s)\n" % (mp_nproc))
        if mp_nproc > 1 or mp_njobs > 1:

            def process_output(result):
                for message in result[1]:
                    logger.log(message.levelno, message.msg)
                reflections.extend(result[0][0])
                result[0][0] = None

            batch_multi_node_parallel_map(
                func=ExtractSpotsParallelTask(function),
                iterable=indices,
                nproc=mp_nproc,
                njobs=mp_njobs,
                cluster_method=mp_method,
                chunksize=mp_chunksize,
                callback=process_output,
            )
        else:
            for task in indices:
                reflections.extend(function(task)[0])

        # Return the reflections
        return reflections, None


class SpotFinder(object):
    """
    A class to do spot finding and filtering.
    """

    def __init__(
        self,
        threshold_function=None,
        mask=None,
        region_of_interest=None,
        max_strong_pixel_fraction=0.1,
        compute_mean_background=False,
        mp_method=None,
        mp_nproc=1,
        mp_njobs=1,
        mp_chunksize=1,
        mask_generator=None,
        filter_spots=None,
        scan_range=None,
        write_hot_mask=True,
        hot_mask_prefix="hot_mask",
        min_spot_size=1,
        max_spot_size=20,
        no_shoeboxes_2d=False,
        min_chunksize=50,
    ):
        """
        Initialise the class.

        :param find_spots: The spot finding algorithm
        :param filter_spots: The spot filtering algorithm
        :param scan_range: The scan range to find spots over
        """

        # Set the filter and some other stuff
        self.threshold_function = threshold_function
        self.mask = mask
        self.region_of_interest = region_of_interest
        self.max_strong_pixel_fraction = max_strong_pixel_fraction
        self.compute_mean_background = compute_mean_background
        self.mask_generator = mask_generator
        self.filter_spots = filter_spots
        self.scan_range = scan_range
        self.write_hot_mask = write_hot_mask
        self.hot_mask_prefix = hot_mask_prefix
        self.min_spot_size = min_spot_size
        self.max_spot_size = max_spot_size
        self.mp_method = mp_method
        self.mp_chunksize = mp_chunksize
        self.mp_nproc = mp_nproc
        self.mp_njobs = mp_njobs
        self.no_shoeboxes_2d = no_shoeboxes_2d
        self.min_chunksize = min_chunksize

    def __call__(self, experiments):
        """
        Do the spot finding.

        :param experiments: The experiments to process
        :return: The observed spots
        """
        from dials.array_family import flex
        import six.moves.cPickle as pickle
        from dxtbx.format.image import ImageBool

        # Loop through all the imagesets and find the strong spots
        reflections = flex.reflection_table()
        for i, experiment in enumerate(experiments):

            imageset = experiment.imageset

            # Find the strong spots in the sweep
            logger.info("-" * 80)
            logger.info("Finding strong spots in imageset %d" % i)
            logger.info("-" * 80)
            logger.info("")
            table, hot_mask = self._find_spots_in_imageset(imageset)
            table["id"] = flex.int(table.nrows(), i)
            reflections.extend(table)

            # Write a hot pixel mask
            if self.write_hot_mask:
                if not imageset.external_lookup.mask.data.empty():
                    for m1, m2 in zip(hot_mask, imageset.external_lookup.mask.data):
                        m1 &= m2.data()
                    imageset.external_lookup.mask.data = ImageBool(hot_mask)
                else:
                    imageset.external_lookup.mask.data = ImageBool(hot_mask)
                imageset.external_lookup.mask.filename = "%s_%d.pickle" % (
                    self.hot_mask_prefix,
                    i,
                )

                # Write the hot mask
                with open(imageset.external_lookup.mask.filename, "wb") as outfile:
                    pickle.dump(hot_mask, outfile, protocol=pickle.HIGHEST_PROTOCOL)

        # Set the strong spot flag
        reflections.set_flags(
            flex.size_t_range(len(reflections)), reflections.flags.strong
        )

        # Check for overloads
        reflections.is_overloaded(experiments)

        # Return the reflections
        return reflections

    def _find_spots_in_imageset(self, imageset):
        """
        Do the spot finding.

        :param imageset: The imageset to process
        :return: The observed spots
        """
        from dials.array_family import flex
        from dxtbx.imageset import ImageSweep

        # The input mask
        mask = self.mask_generator.generate(imageset)
        if self.mask is not None:
            mask = tuple(m1 & m2 for m1, m2 in zip(mask, self.mask))

        # Set the spot finding algorithm
        extract_spots = ExtractSpots(
            threshold_function=self.threshold_function,
            mask=mask,
            region_of_interest=self.region_of_interest,
            max_strong_pixel_fraction=self.max_strong_pixel_fraction,
            compute_mean_background=self.compute_mean_background,
            mp_method=self.mp_method,
            mp_nproc=self.mp_nproc,
            mp_njobs=self.mp_njobs,
            mp_chunksize=self.mp_chunksize,
            min_spot_size=self.min_spot_size,
            max_spot_size=self.max_spot_size,
            filter_spots=self.filter_spots,
            no_shoeboxes_2d=self.no_shoeboxes_2d,
            min_chunksize=self.min_chunksize,
            write_hot_pixel_mask=self.write_hot_mask,
        )

        # Get the max scan range
        if isinstance(imageset, ImageSweep):
            max_scan_range = imageset.get_array_range()
        else:
            max_scan_range = (0, len(imageset))

        # Get list of scan ranges
        if not self.scan_range or self.scan_range[0] is None:
            scan_range = [(max_scan_range[0] + 1, max_scan_range[1])]
        else:
            scan_range = self.scan_range

        # Get spots from bits of scan
        hot_pixels = tuple(flex.size_t() for i in range(len(imageset.get_detector())))
        reflections = flex.reflection_table()
        for j0, j1 in scan_range:
            # Make sure we were asked to do something sensible
            if j1 < j0:
                raise Sorry("Scan range must be in ascending order")
            elif j0 < max_scan_range[0] or j1 > max_scan_range[1]:
                raise Sorry(
                    "Scan range must be within image range {}..{}".format(
                        max_scan_range[0] + 1, max_scan_range[1]
                    )
                )

            logger.info("\nFinding spots in image {0} to {1}...".format(j0, j1))
            j0 -= 1
            if isinstance(imageset, ImageSweep):
                j0 -= imageset.get_array_range()[0]
                j1 -= imageset.get_array_range()[0]
            r, h = extract_spots(imageset[j0:j1])
            reflections.extend(r)
            if h is not None:
                for h1, h2 in zip(hot_pixels, h):
                    h1.extend(h2)

        # Find hot pixels
        hot_mask = self._create_hot_mask(imageset, hot_pixels)

        # Return as a reflection list
        return reflections, hot_mask

    def _create_hot_mask(self, imageset, hot_pixels):
        """
        Find hot pixels in images
        """
        from dials.array_family import flex

        # Write the hot mask
        if self.write_hot_mask:

            # Create the hot pixel mask
            hot_mask = tuple(
                flex.bool(flex.grid(p.get_image_size()[::-1]), True)
                for p in imageset.get_detector()
            )
            num_hot = 0
            if hot_pixels > 0:
                for hp, hm in zip(hot_pixels, hot_mask):
                    for i in range(len(hp)):
                        hm[hp[i]] = False
                    num_hot += len(hp)
            logger.info("Found %d possible hot pixel(s)" % num_hot)

        else:
            hot_mask = None

        # Return the hot mask
        return hot_mask
