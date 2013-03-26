import exceptions
class CentroidException(exceptions.Exception):
    pass

class centroid_interface_prototype(object):
    '''Centroid calculation interface: all of the stuff in here is in pixels.'''

    def __init__(self, reflections):
        '''Initialise the algorithm with the list of reflections.'''

        # Save the reflection list
        self._reflections = reflections

        # Calculate the centroids
        self._compute_centroids()

    def get_reflections(self):
        '''Return the list of reflections'''
        return self._reflections

    def compute_shoebox_centroid(self, shoebox):
        '''Overload me: shoebox has the form of the subset of data which
        should contain all pixels identified within the reflection bounding
        box.'''

        raise RuntimeError, 'overload me'

    def _compute_centroids(self):
        ''' Compute the centroids. For each reflection, call the overloaded
        compute_centroid_from_bbox method.'''

        # Loop through all the refflections and try to calculate the
        # centroid. If a CentroidException is encountered, then pass.
        for i, ref in enumerate(self._reflections):

            try:
                
                # Compute the centroid
                f, r, c, sf, sr, sc = self.compute_shoebox_centroid(
                    ref.shoebox)

                # FIXME observe that bounding_box here has the elements in
                # the reversed order to the one we are aiming for e.g.
                # not frame, then row, then column instead the reverse

                # Add the bounding box offset to the centroid position
                f += ref.bounding_box[0]
                r += ref.bounding_box[2]
                c += ref.bounding_box[4]

                # Add centroid data to reflection
                ref.centroid_position = (c, r, f)
                ref.centroid_variance = (sc, sr, sf)

                # Copy reflection back into array
                self._reflections[i] = ref

            except CentroidException, e:
                continue

        return
