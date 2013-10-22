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



def get_blocks(num_blocks, sweep_length):
    from math import ceil
    blocks = [0]
    block_length = int(ceil(sweep_length / num_blocks))
    for i in range(num_blocks):
        frame = (i + 1) * block_length
        if frame > sweep_length:
            frame = sweep_length
        blocks.append(frame)
    return blocks

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
        from dials.algorithms.integration import ReflectionExtractor2

        # Try importing the command line arguments
        importer = Importer(args)
        if len(importer.imagesets) == 0 or len(importer.imagesets) == 0:
            self.config().print_help()
            return
        elif len(importer.imagesets) > 1:
            raise RuntimeError("Only one imageset can be processed at a time")

        sweep = importer.imagesets[0]
        crystal = importer.crystals[0]

        # Get the blocks
        blocks = get_blocks(options.num_blocks, len(sweep))

        # Predict the reflections
        extract = ReflectionExtractor2(params.integration.shoebox.n_sigma)
        import tarfile
        import cPickle as pickle
        import StringIO
        from time import time

        table = []

        f = tarfile.open(options.output_filename, 'w')
        for i, (b0, b1) in enumerate(zip(blocks[:-1], blocks[1:])):
            print '\nExtracting frames %d to %d\n' % (b0, b1)

            extracted = extract(sweep[b0:b1], crystal)

            minz = min([r.frame_number for r in extracted])
            maxz = max([r.frame_number for r in extracted])

            print int(minz), int(maxz)

            Command.start('Writing file')
            data = StringIO.StringIO(pickle.dumps(extracted, protocol=pickle.HIGHEST_PROTOCOL))
            info = f.tarinfo()
            info.name = 'reflection_%d.p' % i
            info.size = data.len
            info.mtime = int(time())

            f.addfile(info, data)
            Command.end('wrote file')
            table.append((b0, b1))

        data = StringIO.StringIO(pickle.dumps(table, protocol=pickle.HIGHEST_PROTOCOL))
        info = f.tarinfo()
        info.name = 'table.p'
        info.size = data.len
        info.mtime = int(time())
        f.addfile(info, data)
        f.close()

if __name__ == '__main__':
    script = Script()
    script.run()
