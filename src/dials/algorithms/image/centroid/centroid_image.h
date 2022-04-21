/*
 * centroid_image.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_IMAGE_H
#define DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_IMAGE_H

#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/tiny_types.h>
#include <dials/error.h>
#include "centroid_points.h"

namespace dials { namespace algorithms {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int2;
  using scitbx::af::int3;
  using scitbx::af::product;

  /**
   * A class to calculate the centroid of a 2D image.
   */
  template <typename FloatType = double, typename CoordType = vec2<double> >
  class CentroidImage2d : public CentroidPoints<FloatType, CoordType> {
  public:
    // Useful typedefs
    typedef FloatType pixel_type;
    typedef CentroidPoints<FloatType, CoordType> centroid_algorithm_type;
    typedef typename centroid_algorithm_type::value_type value_type;
    typedef typename centroid_algorithm_type::coord_type coord_type;
    typedef typename centroid_algorithm_type::matrix_type matrix_type;

    /**
     * Initialise the algorithm
     * @param image The image pixels
     */
    CentroidImage2d(const af::const_ref<FloatType, af::c_grid<2> > &image)
        : centroid_algorithm_type(image.as_1d(),
                                  generate_coords(image.accessor()).const_ref()) {}

  private:
    /**
     * Generate coordinates.
     * @param size The size of the image
     */
    af::shared<coord_type> generate_coords(af::c_grid<2> size) {
      // Check the sizes
      DIALS_ASSERT(size.all_gt(0));

      // Put all the image coordinates into the array
      af::shared<coord_type> coords(product(size.const_ref()),
                                    af::init_functor_null<coord_type>());
      for (std::size_t j = 0, k = 0; j < size[0]; ++j) {
        for (std::size_t i = 0; i < size[1]; ++i, ++k) {
          coords[k] = coord_type(i + 0.5, j + 0.5);
        }
      }

      // Return the array
      return coords;
    }
  };

  /**
   * A class to calculate the centroid of a 3D image.
   */
  template <typename FloatType = double, typename CoordType = vec3<double> >
  class CentroidImage3d : public CentroidPoints<FloatType, CoordType> {
  public:
    // Useful typedefs
    typedef FloatType pixel_type;
    typedef CentroidPoints<FloatType, CoordType> centroid_algorithm_type;
    typedef typename centroid_algorithm_type::value_type value_type;
    typedef typename centroid_algorithm_type::coord_type coord_type;
    typedef typename centroid_algorithm_type::matrix_type matrix_type;

    /**
     * Initialise the algorithm
     * @param image The image pixels
     */
    CentroidImage3d(const af::const_ref<FloatType, af::c_grid<3> > &image)
        : centroid_algorithm_type(image.as_1d(),
                                  generate_coords(image.accessor()).const_ref()) {}

  private:
    /**
     * Generate coordinates.
     * @param size The size of the image
     */
    af::shared<coord_type> generate_coords(af::c_grid<3> size) {
      // Check the sizes
      DIALS_ASSERT(size.all_gt(0));

      // Put all the image coordinates into the array
      af::shared<coord_type> coords(product(size.const_ref()),
                                    af::init_functor_null<coord_type>());
      for (std::size_t k = 0, l = 0; k < size[0]; ++k) {
        for (std::size_t j = 0; j < size[1]; ++j) {
          for (std::size_t i = 0; i < size[2]; ++i, ++l) {
            coords[l] = coord_type(i + 0.5, j + 0.5, k + 0.5);
          }
        }
      }

      // Return the array
      return coords;
    }
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_IMAGE_H */
