import sys
import h5py
import numpy



from dials.model.experiment import Beam, Goniometer, ScanData, Panel, Detector
from dials.model.data import Reflection, ReflectionList

class NexusWriter(object):
    '''Write all our data to Nexus format.'''

    def __init__(self, filename):
        '''Try to open the given nexus file.'''
        import nxs
        self._handle = h5py.open('temp.h5', 'w')

    def set_beam(self, beam):
        '''Set the beam.'''
        self._make_and_open_group('entry', 'NXentry')
        self._make_and_open_group('instrument','NXinstrument')
        self._write_beam(beam)
        self._handle.closegroup() # NXinstrument
        self._handle.closegroup() # NXentry

    def set_detector(self, detector):
        '''Set the detector.'''
        self._make_and_open_group('entry', 'NXentry')
        self._make_and_open_group('instrument','NXinstrument')
        self._write_detector(detector)
        self._handle.closegroup() # NXinstrument
        self._handle.closegroup() # NXentry

    def set_goniometer(self, goniometer):
        '''Set the goniometer.'''
        self._make_and_open_group('entry', 'NXentry')
        self._make_and_open_group('instrument','NXinstrument')
        self._write_goniometer(goniometer)
        self._handle.closegroup() # NXinstrument
        self._handle.closegroup() # NXentry

    def set_scan(self, scan):
        '''Set the scan.'''
        self._make_and_open_group('entry', 'NXentry')
        self._make_and_open_group('instrument','NXinstrument')
        self._write_scan(scan)
        self._handle.closegroup() # NXinstrument
        self._handle.closegroup() # NXentry

    def set_reflections(self, reflections):
        '''Set the reflections.'''
        self._write_reflection_list(reflections)

    def close(self):
        '''Close the nexus file.'''
        self._handle.close()
        del(self._handle)

    def _make_and_open_group(self, name, nxclass):
        '''Try to open the group. If not found, make then open.'''
        try:
            self._handle[name]
        except KeyError:
            self._handle.makegroup(name, nxclass)
            self._handle.opengroup(name, nxclass)

    def _write_attribute(self, name, value):
        '''Write an attribute to the Nexus file.'''
        import json
        try:
            self._handle.putattr(name, value)
        except TypeError:
            self._handle.putattr(name, json.dumps(value))

    def _write_beam(self, beam):
        '''Open the beam group and write the beam attributes.'''
        self._make_and_open_group('beam', 'NXbeam')
        self._write_attribute('wavelength', beam.get_wavelength())
        self._write_attribute('direction',  beam.get_direction())
        self._handle.closegroup()

    def _write_goniometer(self, gonio):
        '''Open the goniometer group and write the goniometer attributes.'''
        self._make_and_open_group('goniometer', 'NXgoniometer')
        self._write_attribute('rotation_axis', gonio.get_rotation_axis())
        self._handle.closegroup()

    def _write_scan(self, scan):
        '''Open the scan group and write the scan attributes.'''
        self._make_and_open_group('scan', 'NXscan')
        self._write_attribute('image_range',   scan.get_image_range())
        self._write_attribute('oscillation',   scan.get_oscillation())
        self._write_attribute('exposure_time', scan.get_exposure_time())
        self._handle.closegroup()

    def _write_panel(self, panel):
        '''Open the panel group and write the panel attributes.'''
        self._make_and_open_group('panel', 'NXpanel')
        self._write_attribute('type',          panel.get_type())
        self._write_attribute('fast_axis',     panel.get_fast_axis())
        self._write_attribute('slow_axis',     panel.get_slow_axis())
        self._write_attribute('origin',        panel.get_origin())
        self._write_attribute('image_size',    panel.get_image_size())
        self._write_attribute('pixel_size',    panel.get_pixel_size())
        self._write_attribute('trusted_range', panel.get_trusted_range())
        self._handle.closegroup()

    def _write_detector(self, detector):
        '''Open the detector group and write the detector attributes.'''
        self._make_and_open_group('detector', 'NXdetector')
        self._write_attribute('num_panels', len(detector))
        for panel in detector:
            self._write_panel(panel)
        self._handle.closegroup()

    def _write_data(self, name, value):
        '''Write a flex array to the dataset.'''
        value_numpy = value.as_numpy_array()
        if len(value) > 0:
            self._handle.makedata(name, value_numpy.dtype, value_numpy.shape)
            self._handle.opendata(name)
            self._handle.putdata(value_numpy)
            self._handle.closedata()

    def _write_reflection(self, reflection, index):
        '''Open the reflection group and write the reflection data.'''
        self._make_and_open_group('reflection_{0}'.format(index), 'NXdata')
        self._write_attribute('miller_index',      reflection.miller_index)
        self._write_attribute('rotation_angle',    reflection.rotation_angle)
        self._write_attribute('beam_vector',       reflection.beam_vector)
        self._write_attribute('image_coord_mm',    reflection.image_coord_mm)
        self._write_attribute('image_coord_px',    reflection.image_coord_px)
        self._write_attribute('frame_number',      reflection.frame_number)
        self._write_attribute('panel_number',      reflection.panel_number)
        self._write_attribute('bounding_box',      reflection.bounding_box)
        self._write_attribute('centroid_position', reflection.centroid_position)
        self._write_attribute('centroid_variance', reflection.centroid_variance)
        self._write_data('shoebox', reflection.shoebox)
        self._write_data('shoebox_mask', reflection.shoebox_mask)
        self._write_data('transformed_shoebox', reflection.transformed_shoebox)
        self._handle.closegroup()

    def _write_reflection_list(self, reflection_list):
        '''Open the reflection list group and write the reflection data.'''
        self._make_and_open_group('reflection_list', 'NXentry')
        for index, reflection in enumerate(reflection_list):
            self._write_reflection(reflection, index)
        self._handle.closegroup()


class NexusReader(object):
    '''Read all our data from Nexus format.'''

    def __init__(self, filename):
        '''Try to open the given nexus file.'''
        import nxs
        self._handle = nxs.open('temp.h5', 'r')

    def close(self):
        '''Close the nexus file.'''
        self._handle.close()
        del(self._handle)

    def get_beam(self):
        '''Get the beam.'''
        self._handle.opengroup('entry', 'NXentry')
        self._handle.opengroup('instrument','NXinstrument')
        beam = self._read_beam()
        self._handle.closegroup() # NXinstrument
        self._handle.closegroup() # NXentry
        return beam

    def get_detector(self):
        '''Get the detector.'''
        self._handle.opengroup('entry', 'NXentry')
        self._handle.opengroup('instrument','NXinstrument')
        detector = self._read_detector()
        self._handle.closegroup() # NXinstrument
        self._handle.closegroup() # NXentry
        return detector

    def get_goniometer(self):
        '''Get the goniometer.'''
        self._handle.opengroup('entry', 'NXentry')
        self._handle.opengroup('instrument','NXinstrument')
        goniometer = self._read_goniometer()
        self._handle.closegroup() # NXinstrument
        self._handle.closegroup() # NXentry
        return goniometer

    def get_scan(self):
        '''Get the scan.'''
        self._handle.opengroup('entry', 'NXentry')
        self._handle.opengroup('instrument','NXinstrument')
        scan = self._read_scan()
        self._handle.closegroup() # NXinstrument
        self._handle.closegroup() # NXentry
        return scan

    def get_reflections(self):
        '''Get the reflections.'''
        return self._read_reflection_list()

    def _read_beam(self):
        '''Open the beam group and write the beam attributes.'''
        self._handle.open_group('beam', 'NXbeam')
        beam = Beam()
        beam.set_wavelength(self._read_attribute('wavelength'))
        beam.set_direction(self._read_attribute('direction'))
        self._handle.closegroup()
        return beam

    def _read_goniometer(self):
        '''Open the goniometer group and write the goniometer attributes.'''
        self._handle.open_group('goniometer', 'NXgoniometer')
        gonio = Goniometer()
        gonio.set_rotation_axis(self._read_attribute('rotation_axis'))
        self._handle.closegroup()
        return gonio

    def _read_scan(self):
        '''Open the scan group and write the scan attributes.'''
        self._handle.open_group('scan', 'NXscan')
        scan = ScanData()
        scan.set_image_range(self._read_attribute('image_range'))
        scan.set_oscillation(self._read_attribute('oscillation'))
        scan.set_exposure_time(self._read_attribute('exposure_time'))
        self._handle.closegroup()
        return scan

    def _read_panel(self):
        '''Open the panel group and write the panel attributes.'''
        self._handle.open_group('panel', 'NXpanel')
        panel = Panel()
        panel.set_type(self._read_attribute('type',          panel.get_type())
        panel.set_type(self._read_attribute('fast_axis',     panel.get_fast_axis())
        panel.set_type(self._read_attribute('slow_axis',     panel.get_slow_axis())
        self._read_attribute('origin',        panel.get_origin())
        self._read_attribute('image_size',    panel.get_image_size())
        self._read_attribute('pixel_size',    panel.get_pixel_size())
        self._read_attribute('trusted_range', panel.get_trusted_range())
        self._handle.closegroup()
        return panel

    def _read_detector(self, detector):
        '''Open the detector group and write the detector attributes.'''
        self._handle.open_group('detector', 'NXdetector')
        self._write_attribute('num_panels', len(detector))
        for panel in detector:
            self._read_panel(panel)
        self._handle.closegroup()
        return detector

from scitbx.array_family import flex
r1 = Reflection()
r1.shoebox = flex.int(flex.grid(5, 5, 5))
r2 = Reflection()
r2.shoebox = flex.int(flex.grid(5, 5, 5))
rlist = ReflectionList()
rlist.append(r1)
rlist.append(r2)




writer = NexusWriter('temp.h5')
writer.set_beam(Beam())
writer.set_detector(Detector(Panel()))
writer.set_goniometer(Goniometer())
writer.set_scan(ScanData())
writer.set_reflections(rlist)
writer.close()

reader = NexusReader('temp.h5')
reader.get_beam()
reader.get_detector()
reader.get_goniometer()
reader.get_scan()
reader.get_reflections()
reader.close()

#make_and_open_group(nf, 'entry', 'NXentry')
#make_and_open_group(nf, 'instrument','NXinstrument')

#write_beam(nf, Beam())
#write_goniometer(nf, Goniometer())
#write_scan(nf, ScanData())
#write_detector(nf, Detector(Panel()))

#nf.closegroup() # NXinstrument
#nf.closegroup() # NXentry

#from scitbx.array_family import flex
#r1 = Reflection()
#r1.shoebox = flex.int(flex.grid(5, 5, 5))
#r2 = Reflection()
#r2.shoebox = flex.int(flex.grid(5, 5, 5))
#rlist = ReflectionList()
#rlist.append(r1)
#rlist.append(r2)
#write_reflection_list(nf, rlist)

#nf.close()
