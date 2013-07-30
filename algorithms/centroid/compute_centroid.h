/*
 * centroid3d.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_CENTROID_CENTROID3D_H
#define DIALS_ALGORITHMS_CENTROID_CENTROID3D_H

#include <dials/algorithms/image/centroid/centroid_masked_image.h>
#include <dials/model/data/reflection.h>
#include <dials/algorithms/shoebox/mask_code.h>

namespace dials { namespace algorithms {

  using scitbx::af::int6;
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /**
   * Class to calculate centroid for reflections
   */
  class ComputeCentroid {
  public:

    /** Initialise the class */
    ComputeCentroid() {}

    /**
     * Calculate the centroids for the whole reflection list
     * @param reflections The reflection list
     */
    void operator()(ReflectionList &reflections) const {
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        if (reflections[i].is_valid()) {
          try {
            this->operator()(reflections[i]);
          } catch(dials::error) {
            reflections[i].set_valid(false);
          }
        }
      }
    }

    /**
     * Calculate the centroid for the reflection
     * @param reflection The reflection
     */
    void operator()(Reflection &reflection) const {

      // Get the shoebox and mask
      flex_double shoebox = reflection.get_shoebox();
      flex_int shoebox_mask = reflection.get_shoebox_mask();

      // Create the mask from the shoebox mask
      flex_int mask(shoebox_mask.accessor(), 0);
      if (reflection.is_strong()) {
        for (std::size_t i = 0; i < mask.size(); ++i) {
          if (shoebox_mask[i] & shoebox::Strong) {
            mask[i] = 1;
          }
        }
      } else {
        for (std::size_t i = 0; i < mask.size(); ++i) {
          mask[i] = shoebox_mask[i] != 0;
        }
      }

      // Get the centroid offset
      int6 bbox = reflection.get_bounding_box();
      vec3<double> offset(bbox[0], bbox[2], bbox[4]);

      // Compute the centroid and set the stuff
      CentroidMaskedImage3d centroid(shoebox, mask);
      reflection.set_centroid_position(offset + centroid.mean());
      reflection.set_centroid_variance(centroid.standard_error_sq());
      reflection.set_centroid_sq_width(centroid.unbiased_variance());
    }
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_CENTROID_CENTROID3D_H */
