#!/usr/bin/env python
#
# dials.inspect_output.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner


def filter_indices_by_valid(refl):
    from math import sqrt
    index = []
    for i, r in enumerate(refl):
        if r.is_valid():
            index.append(i)
    return index

def filter_indices_by_i_over_sigma(refl, min_i_over_sigma):
    from math import sqrt
    index = []
    for i, r in enumerate(refl):
        if r.is_valid() and r.corrected_intensity_variance > 0:
            ios = r.corrected_intensity / sqrt(r.corrected_intensity_variance)
            if ios > min_i_over_sigma:
                index.append(i)
    return index


def display_centroid_deviations(refl, index):
    ''' Display centroid deviations. '''
    from math import sqrt
    from matplotlib import pylab

    print 'Displaying differences in observed-predicted position'

    # Get the predicted and observed xyz and the difference
    xyz_prd = [refl[i].image_coord_px + (refl[i].frame_number,) for i in index]
    xyz_obs = [refl[i].centroid_position for i in index]
    #xyz_dif = [sqrt(sum((x0-x1)**2 for x0, x1 in zip(*xyz)))
    xyz_dif = [(x0-x1 for x0, x1 in zip(*xyz)) for xyz in zip(xyz_prd, xyz_obs)]

    # Get the individual components
    x, y, z = zip(*xyz_prd)
    xd, yd, zd = zip(*xyz_dif)

    # Plot the output
    pylab.subplot(3, 1, 1)
    pylab.title('Difference between predicted-observed position in x, y, z')
    pylab.scatter(x, xd)
    pylab.axhline(sum(xd) / len(xd))
    pylab.subplot(3, 1, 2)
    pylab.scatter(y, yd)
    pylab.axhline(sum(yd) / len(yd))
    pylab.subplot(3, 1, 3)
    pylab.scatter(z, zd)
    pylab.axhline(sum(zd) / len(zd))
    pylab.show()


def display_profile_correlation(refl, reference, index):
    ''' Display the correlation coefficient between profiles and
    their references.'''
    from scitbx.array_family import flex
    from matplotlib import pylab
    from math import sqrt

    print 'Displaying profile correlations coefficients'

    # Get the correlations for each reflection
    ios = []
    cc = []
    for i in index:

        # Get the profile and reference
        r = refl[i]
        prof_a = r.transformed_shoebox
        prof_b = reference.profile(r.image_coord_px + (r.frame_number,))
        assert(len(prof_a) == len(prof_b))

        # Calculate the correlation
        n = len(prof_a)
        mv_a = flex.mean_and_variance(prof_a.as_1d())
        mv_b = flex.mean_and_variance(prof_b.as_1d())
        ma, sa = mv_a.mean(), mv_a.unweighted_sample_standard_deviation()
        mb, sb = mv_b.mean(), mv_b.unweighted_sample_standard_deviation()
        R = (1.0/(n-1.0)) * flex.sum((prof_a-ma) * (prof_b-mb) / (sa*sb))
        cc.append(R)

        # Calculate the I/sig(I)
        ios.append(r.corrected_intensity / sqrt(r.corrected_intensity_variance))

    # Bin by resolution to get mean and sdev in each bin
    ios_max = int(max(ios)) + 1
    binned_cc = [[] for i in range(ios_max)]
    for i in range(len(ios)):
        idx = int(ios[i])
        binned_cc[idx].append(cc[i])
    binned_ios = [i for i in range(len(binned_cc)) if len(binned_cc[i]) > 0]
    binned_cc = [c for c in binned_cc if len(c) > 0]
    meancc = [sum(c) / len(c) for c in binned_cc if len(c) > 0]
    sdevcc = [sqrt(sum((x - m)**2 for x in c) / len(c))
        for m, c in zip(meancc, binned_cc)]

    # Plot correllation coeffcient vs I/sig(I)
    pylab.subplot(2, 1, 1)
    pylab.title('Correlation coefficients between profile and reference vs I/sig(I)')
    pylab.ylim([0.0, 1.0])
    pylab.scatter(ios, cc)
    pylab.subplot(2, 1, 2)
    pylab.ylim([0.0, 1.0])
    pylab.scatter(binned_ios, meancc)
    pylab.errorbar(binned_ios, meancc, yerr=sdevcc)
    pylab.show()


def display_reference_correlation(reference):
    ''' Display correlations between reference profiles '''
    from scitbx.array_family import flex
    from matplotlib import pylab
    width = 9
    height = len(reference) // 9
    assert(width * height == len(reference))

    # Calculate the correllations between profiles on the same frame
    cc = flex.double(flex.grid(height, width))
    cc2 = flex.double(flex.grid(height, width))
    for j in range(height):
        for i in range(width):
            prof_a = reference.profile(j * width)
            prof_b = reference.profile(j * width + i)
            assert(len(prof_a) == len(prof_b))

            # Calculate the correlation
            n = len(prof_a)
            mv_a = flex.mean_and_variance(prof_a.as_1d())
            mv_b = flex.mean_and_variance(prof_b.as_1d())
            ma, sa = mv_a.mean(), mv_a.unweighted_sample_standard_deviation()
            mb, sb = mv_b.mean(), mv_b.unweighted_sample_standard_deviation()
            R = (1.0/(n-1.0)) * flex.sum((prof_a-ma) * (prof_b-mb) / (sa*sb))
            cc[j,i] = R

            prof_a = reference.profile(i)
            prof_b = reference.profile(j * width + i)
            assert(len(prof_a) == len(prof_b))

            # Calculate the correlation
            n = len(prof_a)
            mv_a = flex.mean_and_variance(prof_a.as_1d())
            mv_b = flex.mean_and_variance(prof_b.as_1d())
            ma, sa = mv_a.mean(), mv_a.unweighted_sample_standard_deviation()
            mb, sb = mv_b.mean(), mv_b.unweighted_sample_standard_deviation()
            R = (1.0/(n-1.0)) * flex.sum((prof_a-ma) * (prof_b-mb) / (sa*sb))
            cc2[j,i] = R


    # Plot the correlations
    y = [k // width for k in range(len(reference))]
    x = [k % width for k in range(len(reference))]

    pylab.suptitle('Correlation coefficients between reference profiles')
    pylab.subplot(1, 2, 1)
    im = pylab.imshow(cc.as_numpy_array().transpose(),
        vmin=0.0, vmax=1.0, interpolation='none')
    pylab.colorbar(im)
    pylab.subplot(1, 2, 2)
    im = pylab.imshow(cc2.as_numpy_array().transpose(),
        vmin=0.0, vmax=1.0, interpolation='none')
    pylab.colorbar(im)
    pylab.show()


def display_spots_per_frame(refl, index):
    ''' Display the number of spots per frame. '''

    from matplotlib import pylab

    x, y, z = zip(*[refl[i].centroid_position for i in index])
    count = [0 for i in range(int(max(z)) + 1)]
    for zz in z:
        zz = int(zz)
        count[zz] += 1

    x0, y0, y0, y1, z0, z1 = zip(*[refl[i].bounding_box for i in index])
    count2 = [0 for i in range(int(max(z1)))]
    for zz0, zz1 in zip(z0, z1):
        zz0 = int(zz0)
        zz1 = int(zz1)
        for zzz in range(zz0, zz1):
            count2[zzz] += 1

    pylab.subplot(2, 1, 1)
    pylab.title('Spot centroids per frame')
    pylab.scatter(range(len(count)), count)

    pylab.subplot(2, 1, 2)
    pylab.title('Spots recorded per frame')
    pylab.scatter(range(len(count2)), count2)
    pylab.show()


def display_spot_sizes(refl, index):
    ''' Display the spot sizes. '''
    from matplotlib import pylab

    bbox = [refl[i].bounding_box for i in index]
    all_sizes = [(b[1]-b[0])*(b[3]-b[2])*(b[5]-b[4]) for b in bbox]
    sizes = sorted(all_sizes)[0:int(0.9*len(all_sizes))]
    counts = [0 for i in range(max(sizes)+1)]
    for s in sizes:
        counts[s] += 1
    pylab.title('Spot sizes (excluding top 10%)')
    pylab.hist(sizes, bins=20)
    pylab.plot(counts, color='black', linewidth=2)
    pylab.show()


def display_reference_profiles(reference, index):
    from matplotlib import pylab
    from math import sqrt, ceil
    from scitbx.array_family import flex
    for i in index:
        profile = reference.profile(i)
        vmin = flex.min(profile)
        vmax = flex.max(profile)
        size = profile.all()
        nimage = size[2]
        nrow = int(sqrt(nimage))
        ncol = int(ceil(nimage / nrow))
        for j in range(size[2]):
            pylab.subplot(nrow, ncol, j)
            pylab.imshow(profile.as_numpy_array()[j], vmin=vmin, vmax=vmax,
                interpolation='none')
        pylab.show()


class Script(ScriptRunner):
    '''A class for running the script.'''

    def __init__(self):
        '''Initialise the script.'''
        from dials.util.command_line import parse_range_list_string

        # The script usage
        usage = "usage: %prog [options] [param.phil] "\
                "{sweep.json | image1.file [image2.file ...]}"

        # Initialise the base class
        ScriptRunner.__init__(self, usage=usage)

        self.config().add_option(
            '--centroid-deviations',
            dest = 'centroid_deviations',
            action = 'store_true', default = False,
            help = 'Display a graph of position-centroid deviations')

        self.config().add_option(
            '--profile-correlations',
            dest = 'profile_correlations',
            action = 'store_true', default = False,
            help = 'Display a graph of correlation coefficients')

        self.config().add_option(
            '--reference-correlations',
            dest = 'reference_correlations',
            action = 'store_true', default = False,
            help = 'Display a graph of correlation coefficients')

        self.config().add_option(
            '--spots-per-frame',
            dest = 'spots_per_frame',
            action = 'store_true', default = False,
            help = 'Display a graph of spots per frame')

        self.config().add_option(
            '--spot-sizes',
            dest = 'spot_sizes',
            action = 'store_true', default = False,
            help = 'Display a graph of spots sizes')

        def range_callback(option, opt_str, value, parser, args=None, kwargs=None):
            setattr(parser.values, option.dest, parse_range_list_string(value))

        self.config().add_option(
            '--reference-profiles',
            dest = 'reference_profiles',
            type = 'string', action = 'callback',
            callback = range_callback,
            help = 'Display a set of reference profiles')

        self.config().add_option(
            '--min-i-over-sigma',
            dest = 'min_i_over_sigma',
            type = 'float', default = 0.0,
            help = 'The minimum I/sig(I) to use')

    def main(self, params, options, args):
        '''Execute the script.'''
        from dials.util.command_line import Importer

        # Try importing the command line arguments
        importer = Importer(args)

        # Display centroid deviations
        if options.centroid_deviations:
            refl = importer.reflections
            display_centroid_deviations(refl,
              filter_indices_by_i_over_sigma(refl, options.min_i_over_sigma))

        # Display profile correlation coefficients
        if options.profile_correlations:
            refl = importer.reflections
            reference = importer.reference
            display_profile_correlation(refl, reference,
              filter_indices_by_i_over_sigma(refl, options.min_i_over_sigma))

        # Display the reference correlation
        if options.reference_correlations:
            reference = importer.reference
            display_reference_correlation(reference)

        # Display the spots per frame
        if options.spots_per_frame:
            refl = importer.reflections
            display_spots_per_frame(refl,
              filter_indices_by_valid(refl))

        # Display the spot sizes
        if options.spot_sizes:
            refl = importer.reflections
            display_spot_sizes(refl, filter_indices_by_valid(refl))

        # Display the reference profiles
        if options.reference_profiles:
            reference = importer.reference
            display_reference_profiles(reference, options.reference_profiles)

if __name__ == '__main__':
    script = Script()
    script.run()
