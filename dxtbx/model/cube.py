#!/usr/bin/env python
# cube.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# A data cube implementation, with aims:
#
# - to be able to read in images (and unload images) silently
# - to be able to pass a 3 dimensional flex array corresponding to a
#   subset of frames, area of detector (shoebox)
# - to be instantiated with a filename template - or through the factory below
# - to behave a little like a flex array

class cube:

    def __init__(self, template):
        '''
        Instantiate cube instance based on template - will use the template
        to find matching files, keep a record of the files available. Template
        should use # in place of frame numbers.
        '''

        import glob
        self._template = template
        self._filenames = glob.glob(template.replace('#', '?'))

        if not self._filenames:
            raise RuntimeError, 'no files found matching %s' % template

        prefix = template.split('#')[0]
        suffix = template.split('#')[-1]

        assert(prefix + '#' * template.count('#') + suffix == template)

        self._format = '%s%%0%dd%s' % (prefix, template.count('#'), suffix)

        self._frames = [int(filename[len(prefix):-len(suffix)]) \
                        for filename in self._filenames]

        # this will contain a dictionary of flex arrays
        self._loaded_frames = { }

        # this will contain references to that dictionary - at the moment
        # keep 25 frames in memory
        from Queue import Queue
        self._cache = Queue(maxsize = 25)

        self._nfast = 0
        self._nslow = 0

        return

    def __repr__(self):
        return '%s\n' % self._template + \
            '%d -> %d' % (min(self._frames), max(self._frames))

    def _load(self, frame_number):

        if frame_number in self._loaded_frames:
            return

        if self._cache.full():
            del(self._loaded_frames[self._cache.get()])

        from iotbx.detectors import ImageFactory
        image = ImageFactory(self._format % frame_number)
        image.read()

        if self._nfast == self._nslow == 0:
            self._nfast = image.parameters['SIZE1']
            self._nslow = image.parameters['SIZE2']
        else:
            assert(self._nfast == image.parameters['SIZE1'])
            assert(self._nslow == image.parameters['SIZE2'])

        self._loaded_frames[frame_number] = image.get_raw_data()
        self._cache.put(frame_number)

        return

    def totals(self):
        '''
        Compute sum of pixel values on frames, to exercise system.
        '''

        totals = { }

        for frame in self._frames:
            self._load(frame)
            totals[frame] = sum(self._loaded_frames[frame])

        return totals

    def get(self, frame0, frame1, slow0, slow1, fast0, fast1):
        '''
        Return a flex array containing double-precision values from the
        region fast0 <= i < fast1, etc.
        '''

        # first validate input - only really useful after first frame loaded

        for frame in range(frame0, frame1):
            assert frame in self._frames

        self._load(frame0)

        assert(fast0 >= 0)
        assert(fast1 <= self._nfast)
        assert(slow0 >= 0)
        assert(slow1 <= self._nslow)

        from scitbx.array_family import flex
        grid = flex.grid((frame1 - frame0, slow1 - slow0, fast1 - fast0))
        shoebox = flex.double(grid)

        # FIXME need to find a faster way to code this: possibly through
        # some additional C++ code?

        for frame in range(frame0, frame1):
            self._load(frame)
            for s in range(slow0, slow1):
                for f in range(fast0, fast1):
                    shoebox[(frame - frame0, s - slow0, f - fast0)] = \
                        self._loaded_frames[frame][(s, f)]

        return shoebox

class cube_factory:
    '''A factory class for cube instances.'''

    @staticmethod
    def from_filename(filename):
        from cube_helpers import TR
        template = TR(filename)[0]
        return cube(template)

if __name__ == '__main__':
    import sys

    c = cube_factory.from_filename(sys.argv[1])
    totals = c.totals()

    for f in sorted(totals):
        print f, totals[f]
