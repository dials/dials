
#ifndef DIALS_INTEGRATION_REFLECTION_MASK_H
#define DIALS_INTEGRATION_REFLECTION_MASK_H

#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>

namespace dials { namespace integration {

/** flex array of vec3 doubles */
typedef scitbx::af::flex <scitbx::vec3 <double> >::type flex_vec3_double;

/** A class representing the reflection mask */
class ReflectionMask {
public:

    /** Default constructor */
    ReflectionMask() {}

    /**
     * Initialise the reflection mask to the given size.
     * @param size The size of the mask (same as detector size)
     * @param roi_size The region of interest about the reflection point
     */
    ReflectionMask(scitbx::vec3 <int> size, 
                   scitbx::vec3 <int> roi_size)
        : mask_(scitbx::af::flex_grid <> (size[0], size[1], size[2])),
                size_(size),
                roi_size_(roi_size) {}

    /**
     * Create the reflection mask. Set the values of the reflection mask in the
     * region of interest around each of the given reflection points to the 
     * reflection index. Exclude pixels from the reflection's mask if they are
     * closer to a neighbouring reflection.
     * @param reflection_xyz
     */
    void create(const flex_vec3_double &reflection_xyz) {

        // Initialise mask to -1
        for (int i = 0; i < mask_.size(); ++i) {
            mask_[i] = -1;
        }

        // Loop through all the reflection detector coordinates given. For each
        // reflection, loop through the mask elements within the reflection's 
        // region of interest. If the mask element value is -1 (i.e. currently
        // unset) then set it to the reflection index. If the mask value is
        // already set to another reflection index, then calculate the distance
        // from the mask point to the currently set reflection xyz point and
        // compare it to the distance between the mask point and the new 
        // reflection xyz point. If the new distance is lower than the old,
        // then set the mask value to the new reflection index.
        int mask_stride_x = size_[2];
        int mask_stride_y = mask_stride_x * size_[1];	
        for (int index = 0; index < reflection_xyz.size(); ++index) {
            scitbx::vec3 <double> xyz = reflection_xyz[index];
            int x0 = (int)xyz[0] - roi_size_[2];
            int x1 = (int)xyz[0] + roi_size_[2];
            int y0 = (int)xyz[1] - roi_size_[1];
            int y1 = (int)xyz[1] + roi_size_[1];
            int z0 = (int)xyz[2] - roi_size_[0];
            int z1 = (int)xyz[2] + roi_size_[0];
            for (int k = z0; k <= z1; ++k) {
                for (int j = y0; j <= y1; ++j) {
                    for (int i = x0; i <= x1; ++i) {
                        int mask_index = i + j * mask_stride_x + k * mask_stride_y;
                        int curr_index = mask_[mask_index];
                        if (curr_index == -1) {
                            mask_[mask_index] = index;
                        } else {
                            scitbx::vec3 <double> point(i, j, k);
                            int distance = (xyz - point).length();
                            int distance_curr = (reflection_xyz[curr_index] - 
                                                 point).length();
                            if (distance < distance_curr) {
                                mask_[mask_index] = index;
                            }
                        }
                    }
                }
            } 
        }
    }

    /** Get the mask array */
    scitbx::af::flex_int get_mask() {
        return mask_;
    }

private:

    scitbx::af::flex_int mask_;
    scitbx::vec3 <int> size_;
    scitbx::vec3 <int> roi_size_;
};

}} // namespace dials::integration

#endif // DIALS_INTEGRATION_REFLECTION_MASK_H
