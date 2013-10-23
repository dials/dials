#!/usr/bin/env python
#
# dials.extract.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner

class Script(ScriptRunner):
    '''A class for running the script.'''

    def __init__(self):
        '''Initialise the script.'''

        # The script usage
        usage = "usage: %prog [options] [param.phil] "\
                "{sweep.json | image1.file [image2.file ...]}"

        # Initialise the base class
        ScriptRunner.__init__(self, usage=usage)

        # The block length
        self.config().add_option(
            '-n', '--num-blocks',
            dest = 'num_blocks',
            type = 'int', default = 1,
            help = 'Set the number of blocks')

        # Output filename option
        self.config().add_option(
            '-o', '--output-filename',
            dest = 'output_filename',
            type = 'string', default = 'extracted.tar',
            help = 'Set the filename for the extracted spots.')

    def main(self, params, options, args):
        '''Execute the script.'''
        from dials.model.serialize import load, dump
        from dials.util.command_line import Command
        from dials.util.command_line import Importer
        from dials.algorithms.integration import ReflectionPredictor
        from dials.algorithms.integration import BlockProfileExtractor
        from dials.algorithms.shoebox import BBoxCalculator
        from dials.model.data import ReflectionList
        from dials.array_family import flex

        # Try importing the command line arguments
        importer = Importer(args)
        if len(importer.imagesets) == 0 or len(importer.imagesets) == 0:
            self.config().print_help()
            return
        elif len(importer.imagesets) > 1:
            raise RuntimeError("Only one imageset can be processed at a time")

        sweep = importer.imagesets[0]
        crystal = importer.crystals[0]

        predict = ReflectionPredictor()
        predicted = predict(sweep, crystal)

        n_sigma = params.integration.shoebox.n_sigma

        # Create the bbox calculator
        compute_bbox = BBoxCalculator(
            sweep.get_beam(), sweep.get_detector(),
            sweep.get_goniometer(), sweep.get_scan(),
            n_sigma * sweep.get_beam().get_sigma_divergence(deg=False),
            n_sigma * crystal.get_mosaicity(deg=False))

        # Calculate the bounding boxes of all the reflections
        Command.start('Calculating bounding boxes')
        compute_bbox(predicted)
        Command.end('Calculated {0} bounding boxes'.format(len(predicted)))

        predicted = ReflectionList(sorted(predicted, key=lambda x: x.frame_number))

        extract = BlockProfileExtractor(sweep, options.num_blocks, predicted)

        for i in range(len(extract)):
            Command.start('Reading')
            extract[i]
            Command.end('Read')

#        from dials.model.serialize import partial_shoebox

#        reader = partial_shoebox.Reader('extracted.tar')
#
#        from time import time
#        st = time()
#        for r in range(len(reader)):
#            stt = time()
#            reader[r]
#            print time() - stt
#
#        print "Total Time: ", time() - st


if __name__ == '__main__':
    script = Script()
    script.run()
