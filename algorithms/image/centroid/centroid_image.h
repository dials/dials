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

#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/error.h>
#include "centroid_points.h"

namespace dials { namespace algorithms {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::product;
  using scitbx::af::flex_double;

  /**
   * A class to calculate the centroid of a 2D image.
   */
  class CentroidImage2d : public CentroidPoints<vec2<double> > {
  public:

    // Useful typedefs
    typedef typename CentroidPoints::coord_type coord_type;
    typedef typename CentroidPoints::matrix_type matrix_type;
    typedef typename CentroidPoints::flex_type flex_type;

    /**
     * Initialise the algorithm
     * @param image The image pixels
     */
    CentroidImage2d(const flex_double &image)
      : CentroidPoints(
          image,
          generate_coords(
            image.accessor().all())) {}

  private:

    /**
     * Generate coordinates.
     * @param size The size of the image
     */
    flex_type generate_coords(const flex_double::index_type &size) {

      // Check the sizes
      DIALS_ASSERT(size.size() == 2);
      DIALS_ASSERT(size[0] > 0 && size[1] > 0);

      // Put all the image coordinates into the array
      flex_type coords(product(size));
      for (std::size_t j = 0, k = 0; j < size[0]; ++j) {
        for (std::size_t i = 0; i < size[1]; ++i, ++k) {
          coords[k] = vec2<double>(i + 0.5, j + 0.5);
        }
      }

      // Return the array
      return coords;
    }
  };

  /**
   * A class to calculate the centroid of a 3D image.
   */
  class CentroidImage3d : public CentroidPoints<vec3<double> > {
  public:

    // Useful typedefs
    typedef typename CentroidPoints::coord_type coord_type;
    typedef typename CentroidPoints::matrix_type matrix_type;
    typedef typename CentroidPoints::flex_type flex_type;

    /**
     * Initialise the algorithm
     * @param image The image pixels
     */
    CentroidImage3d(const flex_double &image)
      : CentroidPoints(
          image,
          generate_coords(
            image.accessor().all())) {}

  private:

    /**
     * Generate coordinates.
     * @param size The size of the image
     */
    flex_type generate_coords(const flex_double::index_type &size) {

      // Check the sizes
      DIALS_ASSERT(size.size() == 3);
      DIALS_ASSERT(size[0] > 0 && size[1] > 0 && size[2] > 0);

      // Put all the image coordinates into the array
      flex_type coords(product(size));
      for (std::size_t k = 0, l = 0; k < size[0]; ++k) {
        for (std::size_t j = 0; j < size[1]; ++j) {
          for (std::size_t i = 0; i < size[2]; ++i, ++l) {
            coords[l] = vec3<double>(i + 0.5, j + 0.5, k + 0.5);
          }
        }
      }

      // Return the array
      return coords;
    }
  };

}}

#endif /* DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_IMAGE_H */
