from __future__ import division
#!/usr/bin/env python
#
# nexus.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

class H5PYEncoder(object):
    '''Encoder base class.'''

    def encode(self, obj, handle):
        '''Encode the object in the HDF5 file.'''
        raise RuntimeError("Overload!")


class H5PYDecoder(object):
    '''Decoder base class.'''

    def decode(self, handle):
        '''Decode the object from the HDF5 file.'''
        raise RuntimeError("Overload!")


class ReflectionListEncoder(H5PYEncoder):
    '''Encoder for the reflection data.'''

    def encode(self, reflections, handle):
        '''Encode the reflection data.'''
        import numpy

        # Create the reflection data group
        group = handle.create_group('entry/data_processing')
        group.attrs['num_reflections'] = len(reflections)

        # Create datasets for each reflection property
        self.create_property(group, reflections, 'miller_index')
        self.create_property(group, reflections, 'rotation_angle')
        self.create_property(group, reflections, 'beam_vector')
        self.create_property(group, reflections, 'image_coord_mm')
        self.create_property(group, reflections, 'image_coord_px')
        self.create_property(group, reflections, 'frame_number')
        self.create_property(group, reflections, 'panel_number')
        self.create_property(group, reflections, 'bounding_box')
        self.create_property(group, reflections, 'centroid_position')
        self.create_property(group, reflections, 'centroid_variance')

        # Create the reflection profile arrays
        self.create_profile(group, reflections, 'shoebox')
        self.create_profile(group, reflections, 'shoebox_mask')
        self.create_profile(group, reflections, 'transformed_shoebox')

    def create_property(self, g, rl, name):
        '''Extract the properties from the reflection list and save.

        Params:
            g The group
            rl The reflection list
            name The name of the property

        '''
        g[name] = [getattr(r, name) for r in rl]

    def create_profile(self, g, rl, name):
        '''Create a profile group and save each profile in its own dataset

        Params:
            g The group
            rl The reflection list
            name The profile name

        '''
        from scitbx.array_family import flex
        subgroup = g.create_group(name)
        for i, r in enumerate(rl):
            dataset_name  = '{0}'.format(i)
            dataset_value = getattr(r, name)
            if len(dataset_value) > 0:
                subgroup[dataset_name] = dataset_value.as_numpy_array()


class ReflectionListDecoder(H5PYDecoder):
    '''Decoder for the reflection data.'''

    def decode(self, handle):
        '''Decode the reflection data.'''
        from dials.model.data import ReflectionList
        from scitbx.array_family import flex
        import numpy

        # Get the group containing the reflection data
        g = handle['entry/data_processing']

        # Create the list of reflections
        rl = ReflectionList(int(g.attrs['num_reflections']))

        # Extract the reflection properties
        self.extract_values(rl, g, 'miller_index',      lambda x: map(int, x))
        self.extract_values(rl, g, 'rotation_angle',    lambda x: float(x))
        self.extract_values(rl, g, 'beam_vector',       lambda x: map(float, x))
        self.extract_values(rl, g, 'image_coord_mm',    lambda x: map(float, x))
        self.extract_values(rl, g, 'image_coord_px',    lambda x: map(float, x))
        self.extract_values(rl, g, 'frame_number',      lambda x: float(x))
        self.extract_values(rl, g, 'panel_number',      lambda x: int(x))
        self.extract_values(rl, g, 'bounding_box',      lambda x: map(int, x))
        self.extract_values(rl, g, 'centroid_position', lambda x: map(float, x))
        self.extract_values(rl, g, 'centroid_variance', lambda x: map(float, x))

        # Extract the reflection profiles
        self.extract_profile(g, rl, 'shoebox', flex.int)
        self.extract_profile(g, rl, 'shoebox_mask', flex.int)
        self.extract_profile(g, rl, 'transformed_shoebox', flex.double)

        # Return the list of reflections
        return rl

    def extract_values(self, rl, g, name, setter=lambda x: x):
        '''Extract a series of values from a dataset.

        Params:
            rl The reflection list
            g The group
            name The name of value to extract
            setter A function for casting the dataset value

        '''
        import numpy

        # Check the length of the dataset is valid
        assert(len(g[name]) == len(rl))

        # Extract the dataset values
        for i, x in enumerate(numpy.array(g[name])):
            setattr(rl[i], name, setter(x));

    def extract_profile(self, g, rl, name, flex_type):
        '''Extract the reflection profiles.

        Params:
            g The group
            rl The reflection list
            name The name of the profile
            flex_type The type of the flex array

        '''
        import numpy
        subgroup = g[name]
        for k, v in subgroup.iteritems():
            if len(v.shape) > 0:
                setattr(rl[int(k)], name, flex_type(numpy.array(v)))


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

    def set_data(self, data, encoder):
        '''Set the model data using the supplied encoder.'''
        encoder.encode(data, self._handle)

    def get_data(self, decoder):
        '''Get the model data using the supplied decoder.'''
        return decoder.decode(self._handle)

    def set_reflections(self, reflections):
        '''Set the reflection data.'''
        self.set_data(reflections, ReflectionListEncoder())

    def get_reflections(self):
        '''Get the reflection data.'''
        return self.get_data(ReflectionListDecoder())


if __name__ == '__main__':
    from dials.model.data import Reflection, ReflectionList
    from scitbx.array_family import flex
    r = Reflection()
    r.shoebox = flex.int(flex.grid(10, 10, 10))

    reflections = ReflectionList()
    reflections.append(r)
    reflections.append(r)
    reflections.append(r)
    reflections.append(r)

    writer = NexusFile('temp.h5', 'w')
    writer.set_reflections(reflections)
    writer.close()

    reader = NexusFile('temp.h5', 'r')
    reflections = reader.get_reflections()
    reader.close()

    for r in reflections:
        print r
        print r.shoebox.all()
