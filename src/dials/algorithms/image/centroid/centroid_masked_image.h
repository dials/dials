/*
 * centroid_masked_image.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_MASKED_IMAGE_H
#define DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_MASKED_IMAGE_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
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
  class CentroidMaskedImage2d : public CentroidPoints<FloatType, CoordType> {
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
    CentroidMaskedImage2d(const af::const_ref<FloatType, af::c_grid<2> > &image,
                          const af::const_ref<bool, af::c_grid<2> > &mask)
        : centroid_algorithm_type(select_pixels(image, mask).const_ref(),
                                  generate_coords(image, mask).const_ref()) {}

  private:
    /**
     * Get the mask indices.
     * @param size The size of the image
     */
    af::shared<FloatType> select_pixels(
      const af::const_ref<FloatType, af::c_grid<2> > &image,
      const af::const_ref<bool, af::c_grid<2> > &mask) {
      // Check the sizes
      DIALS_ASSERT(mask.accessor().all_eq(image.accessor()));
      DIALS_ASSERT(mask.accessor().all_gt(0));

      // Put all the image coordinates into the array
      af::shared<FloatType> pixels(image.size(), af::init_functor_null<FloatType>());
      std::size_t count = 0;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if (mask[i]) {
          pixels[count++] = image[i];
        }
      }
      pixels.resize(count);

      // Return the array
      return pixels;
    }

    /**
     * Generate coordinates.
     * @param size The size of the image
     */
    af::shared<coord_type> generate_coords(
      const af::const_ref<FloatType, af::c_grid<2> > &image,
      const af::const_ref<bool, af::c_grid<2> > &mask) {
      // Check the sizes
      DIALS_ASSERT(mask.accessor().all_eq(image.accessor()));
      DIALS_ASSERT(mask.accessor().all_gt(0));

      // Put all the image coordinates into the array
      af::shared<coord_type> coords(image.size(), af::init_functor_null<coord_type>());
      std::size_t count = 0;
      for (std::size_t j = 0; j < mask.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < mask.accessor()[1]; ++i) {
          if (mask(j, i)) {
            coords[count++] = coord_type(i + 0.5, j + 0.5);
          }
        }
      }
      coords.resize(count);

      // Return the array
      return coords;
    }
  };

  /**
   * A class to calculate the centroid of a 3D image.
   */
  template <typename FloatType = double, typename CoordType = vec3<double> >
  class CentroidMaskedImage3d : public CentroidPoints<FloatType, CoordType> {
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
    CentroidMaskedImage3d(const af::const_ref<FloatType, af::c_grid<3> > &image,
                          const af::const_ref<bool, af::c_grid<3> > &mask)
        : centroid_algorithm_type(select_pixels(image, mask).const_ref(),
                                  generate_coords(image, mask).const_ref()) {}

  private:
    /**
     * Get the mask indices.
     * @param size The size of the image
     */
    af::shared<FloatType> select_pixels(
      const af::const_ref<FloatType, af::c_grid<3> > &image,
      const af::const_ref<bool, af::c_grid<3> > &mask) {
      // Check the sizes
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_gt(0));

      // Put all the image coordinates into the array
      af::shared<FloatType> pixels(image.size(), af::init_functor_null<FloatType>());
      std::size_t count = 0;
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if (mask[i]) {
          pixels[count++] = image[i];
        }
      }
      pixels.resize(count);

      // Return the array
      return pixels;
    }

    /**
     * Generate coordinates.
     * @param size The size of the image
     */
    af::shared<coord_type> generate_coords(
      const af::const_ref<FloatType, af::c_grid<3> > &image,
      const af::const_ref<bool, af::c_grid<3> > &mask) {
      // Check the sizes
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_gt(0));

      // Put all the image coordinates into the array
      af::shared<coord_type> coords(image.size(), af::init_functor_null<coord_type>());
      std::size_t count = 0;
      for (std::size_t k = 0; k < image.accessor()[0]; ++k) {
        for (std::size_t j = 0; j < image.accessor()[1]; ++j) {
          for (std::size_t i = 0; i < image.accessor()[2]; ++i) {
            if (mask(k, j, i)) {
              coords[count++] = coord_type(i + 0.5, j + 0.5, k + 0.5);
            }
          }
        }
      }

      coords.resize(count);

      // Return the array
      return coords;
    }
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_MASKED_IMAGE_H */
