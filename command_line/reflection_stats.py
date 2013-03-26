#
# reflection_stats.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

class ReflectionStats(object):
    '''A class to calculate and print reflection stats.'''
    
    def __init__(self, reflections, adjacency_list=None):
        '''Init the class.
        
        Params:
            reflections The list of reflections
            adjacency_list The adjacency list
        
        '''
        # Save the reflections and adjacency list
        self.reflections = reflections
        self.adjacency_list = adjacency_list
    
        # Compute the stats
        self.compute()
            
    def compute(self):
        '''Compute the statistics'''
        self.stats = list()
        
        # Append # reflections
        self.stats.append(("Num reflections", len(self.reflections)))
        
        # Calculate coordinate stas
        self._compute_coord_stats()
        
        # Calculate bounding box stats
        self._compute_bbox_stats()

        # If adjacency list is given, calculate overlap stats
        if self.adjacency_list:
            self._compute_overlap_stats()    
            
        # Compute centroid statistics
        self._compute_centroid_stats()
    
    def _value_to_string(self, value):
        '''Convert values to string. '''
        
        # If is float then return with a certain precision
        if isinstance(value, float):
            return "{0:.2f}".format(value)
            
        # If a tuple or list then convert each to certain precision
        if isinstance(value, tuple) or isinstance(value, list):
            value_str_tuple = map(self._value_to_string, value)
            value_str = "(" + value_str_tuple[0]
            for v in value_str_tuple[1:]:
                value_str += ", " + v
            value_str += ")"
            return value_str
        
        # Default return as string
        return str(value)        
            
    def __str__(self):
        '''Get the statistics as a string.''' 

        # Calculate max name length for alignment
        max_len = 0
        for name, value in self.stats:
            if len(name) > max_len: 
                max_len = len(name)
        
        # Create the string from the name/value pairs
        stat_str  = '\n'
        stat_str += 'Reflection statistics\n'
        stat_str += '---------------------\n'
        for name, value in self.stats:
            stat_str += '{0:<{1}} : {2}\n'.format(
              name, 
              max_len, 
              self._value_to_string(value))
        
        # Return the string
        return stat_str
        
    def _compute_coord_stats(self):
        '''Compute stats about reflection coordinates.'''
        
        # Extract x, y, z coords
        x = [r.image_coord_px[0] for r in self.reflections]
        y = [r.image_coord_px[1] for r in self.reflections]
        z = [r.frame_number for r in self.reflections]

        # Add min, max to stats
        self.stats.append(("Min (x, y, z)", (min(x), min(y), min(z))))
        self.stats.append(("Max (x, y, z)", (max(x), max(y), max(z))))
        
    def _compute_bbox_stats(self):
        '''Compute bounding box statistics.'''
        import numpy

        # Extract the bounding boxes
        bounding_boxes = [r.bounding_box for r in self.reflections]

        # Calculate the min, max, mean pixels in bounding box
        bbox_count = [(s[1] - s[0]) * (s[3] - s[2]) * (s[5] - s[4]) 
            for s in bounding_boxes]
        min_count = numpy.min(bbox_count)
        max_count = numpy.max(bbox_count)
        med_count = int(numpy.median(bbox_count))

        # Calculate the min, max and median range for x, y, z
        min_range, max_range, med_range = [], [], []
        for ind in ((0, 1), (2, 3), (4, 5)):
            bbox_range = [s[ind[1]] - s[ind[0]] for s in bounding_boxes]
            min_range.append(numpy.min(bbox_range))
            max_range.append(numpy.max(bbox_range))
            med_range.append(numpy.median(bbox_range))

        # Make ranges into tuples
        min_range = tuple(min_range)
        max_range = tuple(max_range)
        med_range = tuple(med_range)

        # Add the stats to the list
        self.stats.append(("Min bounding box count", min_count))
        self.stats.append(("Max bounding box count", max_count))
        self.stats.append(("Median bounding box count", med_count))        
        self.stats.append(("Min bounding box range (x, y, z)", min_range))
        self.stats.append(("Max bounding box range (x, y, z)", max_range))
        self.stats.append(("Median bounding box range (x, y, z)", med_range))
        
    def _compute_overlap_stats(self):
        '''Compute the overlap statistics.'''

        # Loop through all the edges
        min_overlap_x, min_overlap_y, min_overlap_z = 999999, 999999, 999999
        max_overlap_x, max_overlap_y, max_overlap_z = 0, 0, 0
        min_common_pixels = 999999
        max_common_pixels = 0
        for e in self.adjacency_list.edges():
            v1, v2 = self.adjacency_list[e]
            r1, r2 = self.reflections[v1], self.reflections[v2]
            s1, s2 = r1.bounding_box, r2.bounding_box

            # Z overlap
            if s1[0] < s2[0]:
              overlap_z = s1[1] - s2[0]
            else:
              overlap_z = s2[1] - s1[0]

            # Y overlap
            if s1[2] < s2[2]:
              overlap_y = s1[3] - s2[2]
            else:
              overlap_y = s2[3] - s1[2]

            # X overlap
            if s1[4] < s2[4]:
              overlap_x = s1[5] - s2[4]
            else:
              overlap_x = s2[5] - s1[4]

            # calculate the common pixels
            opixels = overlap_x * overlap_y * overlap_z

            # Set min overlap
            if overlap_x < min_overlap_x: min_overlap_x = overlap_x
            if overlap_y < min_overlap_y: min_overlap_y = overlap_y
            if overlap_z < min_overlap_z: min_overlap_z = overlap_z

            # Set max overlap
            if overlap_x > max_overlap_x: max_overlap_x = overlap_x
            if overlap_y > max_overlap_y: max_overlap_y = overlap_y
            if overlap_z > max_overlap_z: max_overlap_z = overlap_z

            # Set min/max common pixels
            if opixels < min_common_pixels: min_common_pixels = opixels
            if opixels > max_common_pixels: max_common_pixels = opixels

        min_overlap_xyz = (min_overlap_x, min_overlap_y, min_overlap_z)
        max_overlap_xyz = (max_overlap_x, max_overlap_y, max_overlap_z)

        # Find the maximum number of overlaps for a reflection
        max_overlap = 0
        for v1 in self.adjacency_list.vertices():
            overlaps = [v2 for v2 in self.adjacency_list.adjacent_vertices(v1)]
            num_overlap = len(overlaps)
            if num_overlap > max_overlap:
                max_overlap = num_overlap                 

        # Add some statistics to the list
        self.stats.append(("Num overlaps", self.adjacency_list.num_edges()))
        self.stats.append(("Max overlap (x, y, z)", max_overlap_xyz))
        self.stats.append(("Min overlap (x, y, z)", min_overlap_xyz))
        self.stats.append(("Max common pixels", max_common_pixels))
        self.stats.append(("Min common pixels", min_common_pixels))
        self.stats.append(("Max number of overlaps", max_overlap))

    def _compute_centroid_stats(self):
        '''Compute some centroid stats.'''
        import numpy
        from math import sqrt
        
        # Get the centroid variance and calculate min/max
        var = [r.centroid_variance for r in self.reflections]
        var_x = [vx for vx, vy, vz in var]
        var_y = [vy for vx, vy, vz in var]
        var_z = [vz for vx, vy, vz in var]       
        min_var = (numpy.min(var_x), numpy.min(var_y), numpy.min(var_z))
        max_var = (numpy.max(var_x), numpy.max(var_y), numpy.max(var_z))
        mean_var = (numpy.mean(var_x), numpy.mean(var_y), numpy.mean(var_z))
        
        # Calculate the difference between the centroid and image coord
        diff = []
        for r in self.reflections:
            diff.append((r.centroid_position[0] - r.image_coord_px[0],
                         r.centroid_position[1] - r.image_coord_px[1],
                         r.centroid_position[2] - r.frame_number))
        
        print diff
        
        # Calculate distance
        dist = [sqrt(d[0]**2 + d[1]**2 + d[2]**2) for d in diff]
        min_diff = numpy.min(dist)
        max_diff = numpy.max(dist)
        mean_diff = numpy.mean(dist)
        
        # Add stats
        self.stats.append(("Min centroid variance (x, y, z)", min_var))
        self.stats.append(("Max centroid variance (x, y, z)", max_var))
        self.stats.append(("Mean centroid variance (x, y, z)", mean_var))
        self.stats.append(("Min centroid-predicted", min_diff))
        self.stats.append(("Max centroid-predicted", max_diff))
        self.stats.append(("Mean centroid-predicted", mean_diff))
