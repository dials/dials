import exceptions
class CentroidException(exceptions.Exception):
    pass

class centroid_interface_prototype(object):
    '''Centroid calculation interface: all of the stuff in here is in pixels.'''

    def __init__(self, reflections):
        self._reflections = reflections
        self._compute_centroids()

        return

    def get_reflections(self):
        return self._reflections

    def compute_centroid_from_bbox(self, shoebox):
        '''Overload me: shoebox has the form of the subset of data which
        should contain all pixels identified within the reflection bounding
        box.'''

        raise RuntimeError, 'overload me'

    def _compute_centroids(self):

        self._centroids = { }

        for i, ref in enumerate(self._reflections):

            try:
                f, r, c, sf, sr, sc = self.compute_centroid_from_bbox(
                    ref.shoebox)

                # FIXME observe that bounding_box here has the elements in
                # the reversed order to the one we are aiming for e.g.
                # not frame, then row, then column instead the reverse

                bounding_box = ref.bounding_box
                f += ref.bounding_box[4]
                r += ref.bounding_box[2]
                c += ref.bounding_box[0]

                # Add centroid data to reflection
                ref.centroid_position = (c, r, f)
                ref.centroid_variance = (sc, sr, sf)

                # Copy reflection back into array
                self._reflections[i] = ref
                
            except CentroidException, e:
                continue

        return
