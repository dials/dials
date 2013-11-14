
class ImageRecord(object):
    def __init__(self, filename=None, mtime=None, beam=None,
                 detector=None, goniometer=None, scan=None):
        self.filename = filename
        self.mtime = mtime
        self.beam = beam
        self.detector = detector
        self.goniometer = goniometer
        self.scan = scan


class DataBlock(object):

    def __init__(self, format_class):
        self.format_class = format_class
        self._images = []

    def append(self, filename):
        from os.path import getmtime

        # Try to read the file
        assert(self.understand(filename))
        fmt = self.format_class(filename)

        # Get the meta data
        b, d, g, s = self._extract_metadata(fmt)

        # Create the image record
        record = ImageRecord(
            filename=filename,
            mtime=getmtime(filename),
            beam=b, detector=d,
            goniometer=g, scan=s)

        # Add the image record
        self._images.append(record)

    def understand(self, filename):
        from dxtbx.format.Format import Format
        mro = self.format_class.mro()[::-1]
        if len(mro) <= 2 or mro[0] != object or mro[1] != Format:
            return False
        for m in mro[2:]:
            if m.understand(filename) == False:
                return False
        return True

    def _extract_metadata(self, fmt):

        # Get the meta data from the format
        try: b = fmt.get_beam()
        except Exception: b = None
        try: d = fmt.get_detector()
        except Exception: d = None
        try: g = fmt.get_goniometer()
        except Exception: g = None
        try: s = fmt.get_scan()
        except Exception: s = None

        # Test against last item in list
        if len(self._images) > 0:
            last = self._images[-1]
            if last.beam == b:
                b = last.beam
            if last.detector == d:
                d = last.detector
            if last.goniometer == g:
                g = last.goniometer
            try:
                last.scan += s
                s = last.scan
            except Exception:
                pass

        # Return the meta data
        return b, d, g, s

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, state):
        from os.path import getmtime
        self.__dict__.update(state)
        for image in self._images:
            if getmtime(image.filename) > image.mtime:
                pass


class DataBlockFactory(object):

    @staticmethod
    def create_list(filenames, verbose=False):
        datablock_list = []
        for f in filenames:
            if verbose: print 'Loading file: %s' % f
            try:
                datablock_list[-1].append(f)
            except Exception:
                datablock_list.append(DataBlockFactory.create([f]))
        return datablock_list

    @staticmethod
    def create(filenames, verbose=False):
        from dxtbx.format.Registry import Registry

        # Ensure we have a list of images
        if len(filenames) < 1:
            raise RuntimeError('Need at least 1 image to create data block')

        # Create the datablock
        datablock = DataBlock(Registry.find(filenames[0]))
        for f in filenames:
            if verbose: print 'Loading file: %s' % f
            datablock.append(f)
        return datablock


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
        ('XDS', 'XPARM.XDS'),
        ('XDS', 'INTEGRATE.HKL'),
        ('XDS', 'XDS_ASCII.HKL')
        ]


if __name__ == '__main__':


    from glob import glob
    from os.path import join
    import cPickle as pickle
    from time import time
    path1 = '/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data'
    path2 = '/home/upc86896/Data/TRP_M1S3_2_'
    path3 = '/home/upc86896/Projects/cctbx/sources/dials_regression/image_examples/'
    filenames3 = [join(path3, *p) for p in get_all_working_images()]

    filenames1 = sorted(glob(join(path1, 'centroid_*.cbf')))
    filenames2 = sorted(glob(join(path2, '*.cbf')))
    filenames = filenames1 + filenames2 + filenames3
    print len(filenames)

    st = time()

    datablock_list = DataBlockFactory.create_list(filenames, verbose=False)

    print time() - st

    print "Dumping"
    pickle.dump(datablock_list, open('temp.p', 'wb'))

    print "Loading"
    datablock_list2 = pickle.load(open('temp.p', 'rb'))
    print "Loaded"

#    datablock = datablock_list[0]
#    images = datablock._images[9:]
#    for i, j in zip(images[0:], images[1:]):
#        assert(i.scan == j.scan)
