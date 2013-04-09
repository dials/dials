#!/usr/bin/env python
#
# NexusFile.py
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
        group = handle.create_group('reflection_data')
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
        g.create_dataset(name, data=[getattr(r, name) for r in rl])

    def create_profile(self, g, rl, name):
        '''Create a profile dataset. The reflection profiles are 3D arrays
        which vary in size depending on the reflection. They are therefore
        saved in the HDF file as a single 1D array containing the data from
        all reflections, an offset inot the array and the shape of the array.

        Params:
            g The group
            rl The reflection list
            name The profile name

        '''
        import numpy

        # Create the data offset and shape arrays
        data = numpy.zeros(shape=0, dtype=numpy.int32)
        offset = numpy.zeros(shape=len(rl) + 1, dtype=numpy.int32)
        shape = numpy.zeros(shape=(len(rl),3), dtype=numpy.int32)

        # Loop through all the reflections. Get the array from the struct
        # and append the data to the data array. Also add the offset and
        # shape to their respective arrays
        for i, r in enumerate(rl):
            flex_array = getattr(r, name)
            data = numpy.append(data, flex_array.as_1d().as_numpy_array())
            offset[i+1] = offset[i] + len(flex_array)
            if len(flex_array):
                shape[i,:] = flex_array.all()
            else:
                shape[i,:] = [0, 0, 0]

        # Create datasets for the profile information
        sub_group = g.create_group(name)
        sub_group.create_dataset('offset', data=offset)
        sub_group.create_dataset('data',   data=data)
        sub_group.create_dataset('shape', data=shape)


class ReflectionListDecoder(H5PYDecoder):
    '''Decoder for the reflection data.'''

    def decode(self, handle):
        '''Decode the reflection data.'''
        from dials.model.data import ReflectionList
        from scitbx.array_family import flex
        import numpy

        # Get the group containing the reflection data
        g = handle['reflection_data']

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
        self.extract_profiles(g, rl, 'shoebox', flex.int)
        self.extract_profiles(g, rl, 'shoebox_mask', flex.int)
        self.extract_profiles(g, rl, 'transformed_shoebox', flex.double)

        # Return the list of reflections
        return reflections

    def extract_values(self, rl, g, name, setter=lambda x: x):
        '''Extract a series of values from a dataset.

        Params:
            rl The reflection list
            g The group
            name The name of value to extract
            setter A function for casting the dataset value

        '''
        # Check the length of the dataset is valid
        assert(len(g[name]) == len(rl))

        # Extract the dataset values
        for i, x in enumerate(g[name]):
            setattr(rl[i], name, setter(x));

    def extract_profiles(self, g, rl, name, flex_type):
        '''Extract the reflection profiles.

        Params:
            g The group
            rl The reflection list
            name The name of the profile
            flex_type The type of the flex array

        '''
        # Get the array of offsets and data
        sub_group = g[name]
        data = sub_group['data']
        offset = sub_group['offset']
        shape = sub_group['shape']

        # Loop through the array of reflections and extract a profile
        for i in range(len(rl)):
            setattr(rl[i], name, self.extract_single(
                data, offset, shape, i, flex_type))

    def extract_single(self, data, offset, shape, index, flex_type):
        '''The profiles are stored in the form of a 1D array containing
        all the data values and an array containing the offset for each
        reflection within the data array. Extract the profile for the
        given reflection and return it as a flex array

        Params:
            data The value dataset
            offset The offset dataset
            index The index of the reflection
            flex_type The flex array type
            shape The shape of array

        Returns:
            A flex array of flex_type containing the profile data

        '''
        import numpy

        # Get the index range
        first, last = offset[index], offset[index+1]

        # Check the offset and return the array
        if last > first:
            array = numpy.array(data[first:last])
            array.shape = shape[index]
            return flex_type(array)
        else:
            assert((shape[index] == (0, 0, 0)).all())
            return flex_type()


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


def create_reflection():
    from dials.model.data import Reflection, ReflectionList
    from scitbx.array_family import flex
    r = Reflection()
    r.miller_index = (1, 2, 3)
    r.beam_vector = (4,5, 6)
    r.shoebox = flex.int(flex.grid(5, 5, 5))
    return r

if __name__ == '__main__':
    from dials.model.data import Reflection, ReflectionList
    from scitbx.array_family import flex

    reflections = ReflectionList()
    reflections.append(create_reflection())
    reflections.append(create_reflection())
    reflections.append(create_reflection())
    reflections.append(create_reflection())

    writer = NexusFile('temp.h5', 'w')
    writer.set_reflections(reflections)
    writer.close()

    reader = NexusFile('temp.h5', 'r')
    reflections = reader.get_reflections()
    reader.close()

    for r in reflections:
        print r
        print r.shoebox.all()
