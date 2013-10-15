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
    pylab.scatter(x, xd)
    pylab.subplot(3, 1, 2)
    pylab.scatter(y, yd)
    pylab.subplot(3, 1, 3)
    pylab.scatter(z, zd)
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
    pylab.ylim([0.0, 1.0])
    pylab.scatter(ios, cc)
    pylab.subplot(2, 1, 2)
    pylab.ylim([0.0, 1.0])
    pylab.scatter(binned_ios, meancc)
    pylab.errorbar(binned_ios, meancc, yerr=sdevcc)
    pylab.show()


class Script(ScriptRunner):
    '''A class for running the script.'''

    def __init__(self):
        '''Initialise the script.'''

        # The script usage
        usage = "usage: %prog [options] [param.phil] "\
                "{sweep.json | image1.file [image2.file ...]}"

        # Initialise the base class
        ScriptRunner.__init__(self, usage=usage)

        self.config().add_option(
            '--show-centroid-deviations',
            dest = 'show_centroid_deviations',
            action = 'store_true', default = False,
            help = 'Display a graph of position-centroid deviations')

        self.config().add_option(
            '--show-profile-cross-correlation',
            dest = 'show_profile_correlation',
            action = 'store_true', default = False,
            help = 'Display a graph of correlation coefficients')

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
        if options.show_centroid_deviations:
            refl = importer.reflections
            display_centroid_deviations(refl,
              filter_indices_by_i_over_sigma(refl, options.min_i_over_sigma))

        # Display profile correlation coefficients
        if options.show_profile_correlation:
            refl = importer.reflections
            reference = importer.reference
            display_profile_correlation(refl, reference,
              filter_indices_by_i_over_sigma(refl, options.min_i_over_sigma))


if __name__ == '__main__':
    script = Script()
    script.run()
