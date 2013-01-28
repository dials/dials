
#ifndef DIALS_INTEGRATION_REFLECTION_MASK_H
#define DIALS_INTEGRATION_REFLECTION_MASK_H

#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/flex_types.h>
#include "../array_family/array_types.h"
#include "../error.h"

namespace dials { namespace integration {

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
    ReflectionMask(scitbx::vec3 <int> mask_size)
        : mask_(scitbx::af::flex_grid <> (
            mask_size[0], 
            mask_size[1], 
            mask_size[2]), 
            -1),
          size_(mask_size) {}

    /** Reset the mask values to -1 */
    void reset_mask() {
        for (int i = 0; i < mask_.size(); ++i) {
            mask_[i] = -1;
        }
    }

    /**
     * Create the reflection mask. Set the values of the reflection mask in the
     * region of interest around each of the given reflection points to the 
     * reflection index. Exclude pixels from the reflection's mask if they are
     * closer to a neighbouring reflection.
     * 
     * @todo The pixels comprising the mask of a reflection are not guarenteed
     *       to be contingious, this should be fixed.
     *
     * @param image_volume_coords The image volume coordinates
     * @param region_of_interest The regions of interest
     * @returns The status True/False for each reflection
     */
    scitbx::af::flex_bool create(
                const dials::af::flex_vec3_double &image_volume_coords,
                const dials::af::flex_tiny6_int &region_of_interest) 
    {
        // Ensure array sizes match
        DIALS_ASSERT(image_volume_coords.size() == region_of_interest.size());

        // Initialise mask to -1
        this->reset_mask();
        
        // Create an array for the status
        scitbx::af::flex_bool status(region_of_interest.size());
        
        // Loop through all the reflection detector coordinates given. For each
        // reflection, loop through the mask elements within the reflection's 
        // region of interest. If the mask element value is -1 (i.e. currently
        // unset) then set it to the reflection index. If the mask value is
        // already set to another reflection index, then calculate the distance
        // from the mask point to the currently set reflection xyz point and
        // compare it to the distance between the mask point and the new 
        // reflection xyz point. If the new distance is lower than the old,
        // then set the mask value to the new reflection index.
        for (int index = 0; index < image_volume_coords.size(); ++index) {
            scitbx::vec3 <double> xyz = image_volume_coords[index];
            scitbx::af::tiny <double, 6> roi = region_of_interest[index];
            if (is_roi_valid(roi)) {
                status[index] = true;
                for (int k = roi[4]; k < roi[5]; ++k) {
                    for (int j = roi[2]; j < roi[3]; ++j) {
                        for (int i = roi[0]; i < roi[1]; ++i) {
                            int curr_index = mask_(k, j, i);
                            if (curr_index == -1) {
                                mask_(k, j, i) = index;
                            } else {
                                scitbx::vec3 <double> point(i, j, k);
                                scitbx::vec3 <double> curr_xyz = 
                                    image_volume_coords[curr_index];
                                int distance = (xyz - point).length();
                                int distance_curr = (curr_xyz - point).length();
                                if (distance < distance_curr) {
                                    mask_(k, j, i) = index;
                                }
                            }
                        }
                    }
                } 
            } else {
                status[index] = false;
            }
        }
        
        // Return the status
        return status;
    }

    /** Get the mask array */
    scitbx::af::flex_int get_mask() {
        return mask_;
    }

private:

    /** Check the roi is valid */
    bool is_roi_valid(scitbx::af::tiny <double, 6> roi) {
        return roi[0] >= 0 && roi[1] < size_[2] &&
               roi[2] >= 0 && roi[3] < size_[1] &&
               roi[4] >= 0 && roi[5] < size_[0];
    }

    scitbx::af::flex_int mask_;
    scitbx::vec3 <int> size_;
};

}} // namespace dials::integration

#endif // DIALS_INTEGRATION_REFLECTION_MASK_H
