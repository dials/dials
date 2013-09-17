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
#include <dials/array_family/scitbx_shared_and_versa.h>

namespace dials { namespace algorithms {

  using scitbx::af::int6;
  using dials::model::Reflection;
  using dials::algorithms::CentroidMaskedImage3d;

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
    void operator()(af::ref<Reflection> reflections) const {
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
      af::const_ref<double, af::c_grid<3> > shoebox =
        reflection.get_shoebox().const_ref();
      af::const_ref<int, af::c_grid<3> > shoebox_mask =
        reflection.get_shoebox_mask().const_ref();

      // Create the mask from the shoebox mask
      af::versa<bool, af::c_grid<3> > mask(shoebox_mask.accessor(), false);
      if (reflection.is_strong()) {
        for (std::size_t i = 0; i < mask.size(); ++i) {
          if (shoebox_mask[i] & shoebox::Strong) {
            mask[i] = true;
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
      CentroidMaskedImage3d centroid(shoebox, mask.const_ref());
      reflection.set_centroid_position(offset + centroid.mean());
      reflection.set_centroid_variance(centroid.standard_error_sq());
      reflection.set_centroid_sq_width(centroid.unbiased_variance());
    }
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_CENTROID_CENTROID3D_H */
