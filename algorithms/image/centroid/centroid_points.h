/*
 * centroid_list.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_LIST_H
#define DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_LIST_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/misc_functions.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  // Useful imports
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::fn::pow2;
  using scitbx::af::sum;
  using scitbx::af::sum_sq;
  using scitbx::af::tiny;
  using scitbx::af::flex;
  using scitbx::af::flex_double;

  /**
   * Class to calculate the centroid of a list of coordinates
   */
  template <typename CoordType>
  class CentroidPoints {
  public:

    // Get the dimensions
    static const std::size_t DIM = CoordType::fixed_size;

    // Useful typedefs
    typedef CoordType coord_type;
    typedef tiny<double, DIM*DIM> matrix_type;
    typedef typename flex<coord_type>::type flex_type;

    /**
     * Calculate the centroid.
     * @param coord The list of coordinates.
     * @param value The list of values
     */
    CentroidPoints(const flex_double &pixels, const flex_type &coords)
      : sum_pixels_(sum(pixels.const_ref())),
        sum_pixels_sq_(sum_sq(pixels.const_ref())),
        sum_pixels_coords_(0.0),
        sum_pixels_delta_sq_(0.0),
        sum_pixels_delta_cross_(0.0) {

      // Check the size of the input
      DIALS_ASSERT(DIM > 1);
      DIALS_ASSERT(coords.size() > 0);
      DIALS_ASSERT(coords.size() == pixels.size());
      DIALS_ASSERT(sum_pixels_ > 0);

      // Calculate the sum of pixel * coords
      for (std::size_t i = 0; i < coords.size(); ++i) {
        sum_pixels_coords_ += pixels[i] * coords[i];
      }

      // Calculate the centroid position
      coord_type m = mean();

      // Calculate the sum of pixels * (coord - mean)**2
      for (std::size_t i = 0; i < coords.size(); ++i) {
        sum_pixels_delta_sq_ += pixels[i] * pow2c(coords[i] - m);
      }

      // Calculate the sum of pixels * (coordA - meanA) * (coordB - meanB)
      // for the cross terms of the covariance matrix
      for (std::size_t i = 0; i < coords.size(); ++i) {
        for (std::size_t j = 0, l = 0; j < DIM - 1; ++j) {
          for (std::size_t k = j + 1; k < DIM; ++k, ++l) {
            sum_pixels_delta_cross_[l] += pixels[i] *
              (coords[i][j] - m[j]) * (coords[i][k] - m[k]);
          }
        }
      }
    }

    /** @returns The sum of the pixel counts. */
    double sum_pixels() const {
      return sum_pixels_;
    }

    /** @returns The sum of the pixels squared */
    double sum_pixels_sq() const {
      return sum_pixels_sq_;
    }

    /** @returns The sum of the pixels * coordinates */
    coord_type sum_pixels_coords() const {
      return sum_pixels_coords_;
    }

    /** @returns The sum of the pixels x (coords - mean)**2 */
    coord_type sum_pixels_delta_sq() const {
      return sum_pixels_delta_sq_;
    }

    /** @returns The sum of pixels x (coordsA - meanA) x (coordsB - meanB) */
    coord_type sum_pixels_delta_cross() const {
      return sum_pixels_delta_cross_;
    }

    /** @returns The centroid position. */
    coord_type mean() const {
      return sum_pixels_coords_ / sum_pixels_;
    }

    /** @returns The biased variance. */
    coord_type variance() const {
      return sum_pixels_delta_sq_ / sum_pixels_;
    }

    /** @returns The unbiased variance. */
    coord_type unbiased_variance() const {
      DIALS_ASSERT(pow2(sum_pixels_) > sum_pixels_sq_);
      return sum_pixels_delta_sq_ * sum_pixels_ /
        (pow2(sum_pixels_) - sum_pixels_sq_);
    }

    /** @returns The biased standard error on the mean squared. */
    coord_type standard_error_sq() const {
      return variance() / sum_pixels_ + 1.0;
    }

    /** @returns The unbiased standard error on the mean squared. */
    coord_type unbiased_standard_error_sq() const {
      return unbiased_variance() / sum_pixels_ + 1.0;
    }

    /** @returns The covariance matrix. */
    matrix_type covariance_matrix() const {
      DIALS_ASSERT(pow2(sum_pixels_) > sum_pixels_sq_);

      // Calculate the scaling of covariance terms.
      double scale = sum_pixels_ / (pow2(sum_pixels_) - sum_pixels_sq_);

      // Calculate the diagonal and cross terms
      coord_type variance = scale * sum_pixels_delta_sq_;
      coord_type covariance = scale * sum_pixels_delta_cross_;

      // Create the covariance matrix
      matrix_type matrix;
      for (std::size_t j = 0, k = 0; j < DIM-1; ++j) {
        matrix[j + j * DIM] = variance[j];
        for (std::size_t i = j+1; i < DIM; ++i, ++k) {
          matrix[i + j * DIM] = covariance[k];
          matrix[j + i * DIM] = covariance[k];
        }
      }

      // Return the covariance matrix
      return matrix;
    }

  private:

    coord_type pow2c(const coord_type &x) const {
      coord_type r;
      for (std::size_t i = 0; i < DIM; ++i) {
        r[i] = x[i] * x[i];
      }
      return r;
    }

    double sum_pixels_;
    double sum_pixels_sq_;
    coord_type sum_pixels_coords_;
    coord_type sum_pixels_delta_sq_;
    coord_type sum_pixels_delta_cross_;
  };
}}

#endif /* DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_LIST_H */
