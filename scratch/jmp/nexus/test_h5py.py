#!/usr/bin/env python
#
# NexusFile.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

class H5PYEncoder(object):
    '''Encoder base class.'''
    
    def flex_array(self, arr):
        '''Encode a flex array for a dataset.'''
        if len(arr):
            return arr.as_numpy_array()
        else:
            return []
    
    def encode(self, obj, handle):
        '''Encode the object in the HDF5 file.'''
        raise RuntimeError("Overload!")


class H5PYDecoder(object):
    '''Decoder base class.'''
    
    def flex_array(self, dset, flex_type):
        '''Decode a flex array from a dataset.'''
        import numpy
        if len(dset):
            return flex_type(numpy.array(dset))
        else:
            return flex_type()
    
    def decode(self, handle):
        '''Decode the object from the HDF5 file.'''
        raise RuntimeError("Overload!")
               

class BeamEncoder(H5PYEncoder):
    '''Encoder for the beam model.'''

    def encode(self, obj, handle):
        '''Encode the beam model.'''
        group = handle.create_group('experiment/beam')
        group.attrs['direction']  = obj.get_direction()
        group.attrs['wavelength'] = obj.get_wavelength()
                

class BeamDecoder(H5PYDecoder):
    '''Decoder for the beam model.'''
    
    def decode(self, handle):
        '''Decode the beam model.'''
        from dials.model.experiment import Beam
        group = handle['experiment/beam']
        direction  = group.attrs['direction']
        wavelength = group.attrs['wavelength']
        return Beam(direction, wavelength)
        

class GoniometerEncoder(H5PYEncoder):
    '''Encoder for the goniometer model.'''
    
    def encode(self, obj, handle):
        '''Encode the goniometer model.'''
        group = handle.create_group('experiment/goniometer')
        group.attrs['rotation_axis']  = obj.get_rotation_axis()
        group.attrs['fixed_rotation'] = obj.get_fixed_rotation()


class GoniometerDecoder(H5PYDecoder):
    '''Decoder for the goniometer model.'''
    
    def decode(self, handle):
        '''Decode the goniometer model.'''
        from dials.model.experiment import Goniometer
        group = handle['experiment/goniometer']
        rotation_axis  = group.attrs['rotation_axis']
        fixed_rotation = group.attrs['fixed_rotation']
        return Goniometer(rotation_axis, fixed_rotation)


class DetectorEncoder(H5PYEncoder):
    '''Encoder for the detector model.'''
    
    def encode(self, obj, handle):
        '''Encode the detector model.'''
        group = handle.create_group('experiment/detector')
        group.attrs['num_panels'] = len(obj)
        for i, p in enumerate(obj): 
            self.encode_panel(p, i, group)

    def encode_panel(self, p, i, handle):
        '''Encode the panel model.'''
        group = handle.create_group('panel_{0}'.format(i))
        group.attrs['type']          = p.get_type()
        group.attrs['fast_axis']     = p.get_fast_axis()
        group.attrs['slow_axis']     = p.get_slow_axis()
        group.attrs['origin']        = p.get_origin()
        group.attrs['image_size']    = p.get_image_size()
        group.attrs['pixel_size']    = p.get_pixel_size()
        group.attrs['trusted_range'] = p.get_trusted_range()


class DetectorDecoder(H5PYDecoder):
    '''Decoder for the detector model.'''
    
    def decode(self, handle):
        '''Decode the detector model.'''
        from dials.model.experiment import Detector
        import numpy
        group = handle['experiment/detector']
        panel_list = []
        for i in range(group.attrs['num_panels']):
            panel_list.append(self.decode_panel(i, group))
        return Detector(panel_list)

    def decode_panel(self, i, handle):
        '''Decode the panel model.'''
        from dials.model.experiment import Panel
        group = handle['panel_{0}'.format(i)]
        stype         = group.attrs['type']
        fast_axis     = group.attrs['fast_axis']
        slow_axis     = group.attrs['slow_axis']
        origin        = group.attrs['origin']
        pixel_size    = group.attrs['pixel_size']
        image_size    = map(int, group.attrs['image_size'])
        trusted_range = group.attrs['trusted_range']
        return Panel(stype, fast_axis, slow_axis, origin, 
            pixel_size, image_size, trusted_range)
               

class ScanEncoder(H5PYEncoder):
    '''Encoder for the scan model.'''
    
    def encode(self, obj, handle):
        '''Encode the scan model.'''
        group = handle.create_group('experiment/scan')
        group.attrs['image_range']   = obj.get_image_range()
        group.attrs['oscillation']   = obj.get_oscillation()
        group.attrs['exposure_time'] = obj.get_exposure_time()
        
        
class ScanDecoder(H5PYDecoder):
    '''Decoder for the scan model.'''
    
    def decode(self, handle):
        '''Decode the scan model.'''
        from dials.model.experiment import ScanData
        group = handle['experiment/scan']
        image_range   = map(int, group.attrs['image_range'])
        oscillation   = group.attrs['oscillation']
        exposure_time = group.attrs['exposure_time']
        return ScanData(image_range, oscillation, exposure_time)
                

class ReflectionListEncoder(H5PYEncoder):
    '''Encoder for the reflection data.'''

    def encode(self, obj, handle):
        '''Encode the reflection data.'''    
    
        group = handle.create_group('data/reflections')
        group.attrs['num_reflections'] = len(obj)
        for i, r in enumerate(obj):
            self.encode_reflection(r, i, group)
      
    def encode_reflection(self, r, i, handle):
        '''Encode the single reflectio's data.'''
        group = handle.create_group('r_{0}'.format(i))
        group.attrs['miller_index']      = r.miller_index
        group.attrs['rotation_angle']    = r.rotation_angle
        group.attrs['beam_vector']       = r.beam_vector
        group.attrs['image_coord_mm']    = r.image_coord_mm
        group.attrs['image_coord_px']    = r.image_coord_px
        group.attrs['frame_number']      = r.frame_number
        group.attrs['panel_number']      = r.panel_number
        group.attrs['bounding_box']      = r.bounding_box
        group.attrs['centroid_position'] = r.centroid_position
        group.attrs['centroid_variance'] = r.centroid_variance
        
        # Create datasets for the shoebox, mask and transformed shoebox
        group.create_dataset('shoebox', 
            data=self.flex_array(r.shoebox))
        group.create_dataset('shoebox_mask', 
            data=self.flex_array(r.shoebox_mask))
        group.create_dataset('transformed_shoebox',
            data=self.flex_array(r.transformed_shoebox))


class ReflectionListDecoder(H5PYDecoder):
    '''Decoder for the reflection data.'''
    
    def decode(self, handle):
        '''Decode the reflection data.'''    
        from dials.model.data import ReflectionList    
        group = handle['data/reflections']
        num_reflections = group.attrs['num_reflections']
        reflections = ReflectionList(int(num_reflections))
        for i in range(num_reflections):
            reflections[i] = self.decode_reflection(i, group)
        return reflections
            
    def decode_reflection(self, i, handle):
        '''Decode the single reflection's data.'''
        from dials.model.data import Reflection
        group = handle['r_{0}'.format(i)]        
        r = Reflection()
        r.miller_index      = map(int, group.attrs['miller_index'])
        r.rotation_angle    = group.attrs['rotation_angle']
        r.beam_vector       = group.attrs['beam_vector']
        r.image_coord_mm    = group.attrs['image_coord_mm']
        r.image_coord_px    = group.attrs['image_coord_px']
        r.frame_number      = group.attrs['frame_number']
        r.panel_number      = int(group.attrs['panel_number'])
        r.bounding_box      = map(int, group.attrs['bounding_box'])
        r.centroid_position = group.attrs['centroid_position']
        r.centroid_variance = group.attrs['centroid_variance']
       
        # Get the shoebox, mask and transformed shoebox
        r.shoebox = self.flex_array(group['shoebox'], flex.int)
        r.shoebox_mask = self.flex_array(group['shoebox_mask'], flex.int)
        r.transformed_shoebox = self.flex_array(
            group['transformed_shoebox'], flex.double)
        return r


class NexusFile(object):
    '''Interface to Nexus file.'''

    def __init__(self, filename, mode='a'):
        '''Open the file with the given mode.'''
        import h5py
        self._handle = h5py.File(filename, mode)
       
    def close(self):
        '''Close the file.'''
        self._handle.close()
        del(self._handle)        

    def set_model(self, model, encoder):
        '''Set the model data using the supplied encoder.'''
        encoder.encode(model, self._handle)
        
    def get_model(self, decoder):
        '''Get the model data using the supplied decoder.'''
        return decoder.decode(self._handle)  

    def get_beam(self):
        '''Get the beam model.'''
        return self.get_model(BeamDecoder())
        
    def set_beam(self, beam):
        '''Set the beam model.'''
        self.set_model(beam, BeamEncoder())

    def get_goniometer(self):
        '''Get the goniometer model.'''
        return self.get_model(GoniometerDecoder())
        
    def set_goniometer(self, gonio):
        '''Set the goniometer model.'''
        self.set_model(gonio, GoniometerEncoder())
 
    def get_detector(self):
        '''Get the detector model.'''
        return self.get_model(DetectorDecoder())       

    def set_detector(self, detector):
        '''Set the detector model.'''
        self.set_model(detector, DetectorEncoder())
    
    def get_scan(self):
        '''Get the scan model.'''
        return self.get_model(ScanDecoder())
    
    def set_scan(self, scan):
        '''Set the scan model.'''
        self.set_model(scan, ScanEncoder())
    
    def get_reflections(self):
        '''Get the reflection data.'''
        return self.get_model(ReflectionListDecoder())
        
    def set_reflections(self, reflections):
        '''Set the reflection data.'''
        self.set_model(reflections, ReflectionListEncoder())
    
    
if __name__ == '__main__':
    from dials.model.experiment import Beam
    from dials.model.experiment import Goniometer
    from dials.model.experiment import Panel, Detector
    from dials.model.experiment import ScanData
    from dials.model.data import Reflection, ReflectionList
    from scitbx.array_family import flex

    beam = Beam((1, 2, 3))
    gonio = Goniometer((1, 2, 3))
    p = Panel()
    p.set_frame((1, 0, 0), (0, 1, 0), (0, 0, 1))
    detector = Detector(p)
    scan = ScanData()

    reflections = ReflectionList()
    r = Reflection()
    r.shoebox = flex.int(flex.grid(5, 5, 5))
    reflections.append(r)
    reflections.append(Reflection())
    reflections.append(Reflection())
    reflections.append(Reflection())

    print len(reflections[0].shoebox)

    writer = NexusFile('temp.h5', 'w')
    writer.set_beam(beam)
    writer.set_goniometer(gonio)
    writer.set_detector(detector)
    writer.set_reflections(reflections)
    writer.set_scan(scan)
    writer.close()

    reader = NexusFile('temp.h5', 'r')
    print reader.get_beam()
    print reader.get_goniometer()
    print reader.get_detector()
    print reader.get_scan()
    reflections = reader.get_reflections()
    reader.close()

    print reflections[0].shoebox.all()
    for r in reflections: print r

