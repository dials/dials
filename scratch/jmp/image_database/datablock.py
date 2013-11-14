
class ImageRecord(object):
    ''' Record for storing image metadata. '''
    def __init__(self, mtime=None, beam=None, detector=None,
                 goniometer=None, scan=None, template=None, index=None):
        self.mtime = mtime
        self.beam = beam
        self.detector = detector
        self.goniometer = goniometer
        self.scan = scan
        self.template = template
        self.index = index

    def clone(self, rhs):
        self.mtime = rhs.mtime
        self.beam = rhs.beam
        self.detector = rhs.detector
        self.goniometer = rhs.goniometer
        self.scan = rhs.scan
        self.template = rhs.template
        self.index = rhs.index


class DataBlock(object):
    ''' High level container for blocks of sweeps and imagesets. '''

    def __init__(self, filenames=None, format_class=None):
        ''' Instantiate from filenames or a format class.

        If no format is given then the format is found from the first
        filename. All files must use the same format, otherwise an
        exception will be raised.

        Params:
            filenames The filenames to use
            format_class The format class to use

        '''
        from dxtbx.format.Registry import Registry
        from collections import OrderedDict

        # Try to get a format class
        if format_class == None:
            if filenames is not None and len(filenames) > 0:
                format_class = Registry.find(filenames[0])
        self._format_class = format_class

        # Load any images
        self._images = OrderedDict()
        if filenames is not None:
            for f in filenames:
                self.append(f)

    def format_class(self):
        ''' Return the format class. '''
        return self._format_class

    def filenames(self):
        ''' Return the list of filenames. '''
        return list(self._images.iterkeys())

    def metadata(self):
        ''' Return the list of filenames and meta data. '''
        return list(self._images.iteritems())

    def extract_all(self):
        ''' Extract all the images as an image set. '''
        from dxtbx.imageset import ImageSetFactory
        return ImageSetFactory.make_imageset(
            self._images.keys(), self._format_class)

    def extract_stills(self):
        ''' Extract all the still images as an image set. '''
        from dxtbx.imageset import ImageSetFactory
        stills = [f for f, i in self._images.iteritems() if i.template is None]
        if len(stills) == 0:
            return None
        return ImageSetFactory.make_imageset(stills, self._format_class)

    def extract_sweeps(self):
        ''' Extract all the sweeps from the block. '''
        from dxtbx.imageset import ImageSetFactory
        from itertools import groupby
        sweeps = []

        # Get consecutive groups of scans (which should correspond to sweeps)
        groups = groupby(self._images.itervalues(), key=lambda r: r.scan)
        for scan, records in groups:
            records = list(records)
            if scan is not None and len(records) > 0:
                templates = [r.template for r in records]
                indices = [r.index for r in records]
                assert(templates[0] is not None and len(indices) > 0)
                assert(templates.count(templates[0]) == len(templates))
                assert(all(j == i+1 for i, j in zip(indices[:-1], indices[1:])))
                sweep = ImageSetFactory.make_sweep(
                    templates[0], indices, self._format_class)
                sweeps.append(sweep)

        # Return the list of sweeps
        return sweeps

    def append(self, filename):
        ''' Add another image to the block. The image must use the same
        format class, otherwise an exception will be raised. '''
        from os.path import abspath

        # Check the image is not already here and can be understood
        filename = abspath(filename)
        if filename in self._images:
            raise RuntimeError('%s already in data block' % filename)
        if not self.understand(filename):
            raise RuntimeError('cannot understand file: %s' % filename)

        # Get the meta data
        if len(self._images) > 0:
            last = self._images[next(reversed(self._images))]
        else:
            last = None
        self._images[filename] = self._create_record(filename, last)


    def understand(self, filename):
        ''' Check if the data block format understands the given file. This
        function checks the method resolution order for the format class and
        calls the understand methods of each parent down to the bottom level.
        Just calling the format class understand method directly can result
        in problems. This is really a workaround for a bug in the way that
        the format understand method works. '''
        from dxtbx.format.Format import Format
        mro = self._format_class.mro()[::-1]
        if len(mro) <= 2 or mro[0] != object or mro[1] != Format:
            return False
        for m in mro[2:]:
            if m.understand(filename) == False:
                return False
        return True

    def _create_record(self, filename, last=None):
        ''' Get any information about the image we can and create a record. '''
        from dxtbx.sweep_filenames import template_regex
        from os.path import getmtime

        # Read the image
        fmt = self._format_class(filename)

        # Get the meta data from the format
        try: b = fmt.get_beam()
        except Exception: b = None
        try: d = fmt.get_detector()
        except Exception: d = None
        try: g = fmt.get_goniometer()
        except Exception: g = None
        try: s = fmt.get_scan()
        except Exception: s = None

        # Get the template and index if possible
        if s is not None:
            template, index = template_regex(filename)
        else:
            template, index = None, None

        # Test against last item in list
        if last is not None:
            if last.beam == b:
                b = last.beam
            if last.detector == d:
                d = last.detector
            if last.goniometer == g:
                g = last.goniometer
            try:
                if last.template == template and last.index + 1 == index:
                    last.scan += s
                    s = last.scan
            except Exception:
                pass

        # Create the record and return
        return ImageRecord(
            mtime=getmtime(filename),
            beam=b, detector=d, goniometer=g, scan=s,
            template=template, index=index)

    def __getstate__(self):
        ''' Return the dict for pickling. '''
        return self.__dict__

    def __setstate__(self, state):
        ''' Update the dict from pickling. On reload, update any records for
        images that have changed since the class was pickled. '''
        from os.path import getmtime
        self.__dict__.update(state)
        last = None
        for filename, image in self._images.iteritems():
            if getmtime(filename) > image.mtime:
                image.clone(self._create_record(filename, last))
            last = image


class DataBlockFactory(object):
    ''' Class for creating DataBlock instances'''

    @staticmethod
    def from_filenames(filenames, verbose=False):
        ''' Create a list of data blocks from a list of filenames. '''
        return DataBlockFactory.create_list(filenames, verbose)

    @staticmethod
    def create_list(filenames, verbose=False):
        ''' Create a list of data blocks from a list of filenames. '''
        datablock_list = []
        for f in filenames:
            if verbose: print 'Loading file: %s' % f
            try:
                datablock_list[-1].append(f)
            except Exception:
                datablock_list.append(DataBlockFactory.create_single([f]))
        return datablock_list

    @staticmethod
    def create_single(filenames, verbose=False):
        ''' Create a single data blocks from a list of filenames. '''

        # Ensure we have a list of images
        if len(filenames) < 1:
            raise RuntimeError('Need at least 1 image to create a data block')

        # Create the datablock
        return DataBlock(filenames)


def get_all_working_images():
    return [
        ('ALS_1231', 'q315r_lyso_1_001.img'),
        ('ALS_501', 'als501_q4_1_001.img'),
        ('ALS_821', 'q210_lyso_1_101.img'),
        ('ALS_831', 'q315r_lyso_001.img'),
        ('APS_14BMC', 'q315_1_001.img'),
        ('APS_17ID', 'q210_1_001.img'),
        ('APS_19ID', 'q315_unbinned_a.0001.img'),
        ('APS_22ID', 'mar300.0001'),
        ('APS_23IDD', 'mar300_1_E1.0001'),
        ('APS_24IDC', 'pilatus_1_0001.cbf'),
        ('APS_24IDC', 'q315_1_001.img'),
        ('CLS1_08ID1', 'mar225_2_E0_0001.img'),
        ('DESY_ID141', 'q210_2_001.img'),
        ('ESRF_BM14', 'mar165_001.mccd'),
        ('ESRF_BM14', 'mar225_1_001.mccd'),
        ('ESRF_ID231', 'q315r_7_001.img'),
        ('RAXIS-HTC', 'test1_lysozyme_0111060001.osc'),
        ('SLS_X06SA', 'mar225_2_001.img'),
        ('SLS_X06SA', 'pilatus6m_1_00001.cbf'),
        ('SRS_101', 'mar225_001.img'),
        ('SRS_142', 'q4_1_001.img'),
        ('SSRL_bl111', 'mar325_1_001.mccd'),
        ('xia2', 'merge2cbf_averaged_0001.cbf'),
#        ('XDS', 'XPARM.XDS'),
#        ('XDS', 'INTEGRATE.HKL'),
#        ('XDS', 'XDS_ASCII.HKL')
        ]


if __name__ == '__main__':

    from dxtbx import model
    from glob import glob
    from os.path import join
    import pickle
    from time import time
    path1 = '/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data'
    path2 = '/home/upc86896/Data/TRP_M1S3_2_'
    path3 = '/home/upc86896/Projects/cctbx/sources/dials_regression/image_examples/'
    filenames3 = [join(path3, *p) for p in get_all_working_images()]

    filenames1 = sorted(glob(join(path1, 'centroid_*.cbf')))
    filenames2 = sorted(glob(join(path2, '*.cbf')))
    filenames = filenames1 + filenames2 + filenames3
    print len(filenames)

#    st = time()

    datablock_list = DataBlockFactory.from_filenames(filenames, verbose=True)

#    print time() - st

    print "Dumping"
    pickle.dump(datablock_list, open('temp.p', 'wb'))
#
    print "Loading"
    datablock_list2 = pickle.load(open('temp.p', 'rb'))
    print "Loaded"

#    print datablock_list2[0].filenames()
#    print datablock_list2[0].metadata()

    imageset = datablock_list2[0].extract_all()
#    print len(imageset)
#
#
#    for d in datablock_list:
#        print type(d.extract_stills())



#    for d in datablock_list2:
#        d.extract_sweeps()

#    imageset = datablock_list2[0].extract_sweeps()[0]
#    for i in imageset[0:2]:
#        pass
    pickle.dump(imageset, open("imageset.p", "wb"))
    pickle.load(open("imageset.p", "rb"))
#    print len(imageset)
#    for i in imageset:
#        print len(i)

#    datablock = datablock_list[0]
#    images = datablock._images[9:]
#    for i, j in zip(images[0:], images[1:]):
#        assert(i.scan == j.scan)
