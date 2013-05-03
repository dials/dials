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
#include <dials/error.h>

namespace dials { namespace algorithms {

  // Useful imports
  using scitbx::af::sum;
  using scitbx::af::sum_sq;
  using scitbx::fn::pow2;
  using scitbx::af::tiny;
  using scitbx::af::flex;
  using scitbx::af::flex_double;

  /**
   * Class to calculate the centroid of a list of coordinates
   */
  template <typename CoordType>
  class CentroidList {
  public:

    // Useful typedefs
    typedef CoordType coord_type;
    typedef typename flex<coord_type>::type flex_type;

    /**
     * Calculate the centroid.
     * @param coord The list of coordinates.
     * @param value The list of values
     */
    CentroidList(const flex_double &pixels, const flex_type &coords)
      : sum_pixels_(sum(pixels.const_ref())),
        sum_pixels_sq_(sum_sq(pixels.const_ref())),
        sum_pixels_coords_(0.0),
        sum_pixels_delta_sq_(0.0),
        sum_pixels_delta_cross_(0.0) {

      // Check the size of the input
      DIALS_ASSERT(coords.size() > 0);
      DIALS_ASSERT(coords[0].size() > 1);
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
        sum_pixels_delta_sq_ += pixels[i] * pow2(coords[i] - m);
      }

      // Calculate the sum of pixels * (coordA - meanA) * (coordB - meanB)
      // for the cross terms of the covariance matrix
      for (std::size_t i = 0; i < coords.size(); ++i) {
        for (std::size_t j = 0, l = 0; j < coords[0].size() - 1; ++j) {
          for (std::size_t k = j + 1; k < coords[0].size(); ++k, ++l) {
            sum_pixels_delta_cross_[l] += pixels[i] *
              (coords[i][j] - m[j]) * (coords[i][k] - m[k]);
          }
        }
      }
    }

    double sum_pixels() const {
      return sum_pixels_;
    }

    double sum_pixels_sq() const {
      return sum_pixels_sq_;
    }

    coord_type sum_pixels_coords() const {
      return sum_pixels_coords_;
    }

    coord_type sum_pixels_delta_sq() const {
      return sum_pixels_delta_sq_;
    }

    coord_type sum_pixels_delta_cross() const {
      return sum_pixels_delta_cross_;
    }

    coord_type mean() const {
      return sum_pixels_coords_ / sum_pixels_;
    }

    coord_type biased_variance() const {
      return sum_pixels_delta_sq_ / sum_pixels_;
    }

    coord_type unbiased_variance() const {
      DIALS_ASSERT(pow2(sum_pixels_) > sum_pixels_sq_);
      return sum_pixels_delta_sq_ * sum_pixels_ /
        (pow2(sum_pixels) - sum_pixels_sq_);
    }

    coord_type biased_standard_error_sq() const {
      return biased_variance() / sum_pixels_ + 1.0;
    }

    coord_type unbiased_standard_error_sq() const {
      return unbiased_variance() / sum_pixels_ + 1.0;
    }

    tiny<double, 9> covariance_matrix() const {
      DIALS_ASSERT(pow2(sum_pixels_) > sum_pixels_sq_);
      double scale = sum_pixels_ / (pow2(sum_pixels) - sum_pixels_sq_);

      coord_type variance = unbiased_variance();
      coord_type covariance = scale * sum_pixels_delta_cross_;

      const int DIM = 3;

      tiny<double, 9> matrix;
      for (std::size_t j = 0, k = 0; j < DIM-1; ++j) {
        matrix[j + j * DIM] = variance[j];
        for (std::size_t i = j+1; i < DIM; ++i, ++k) {
          matrix[i + j * DIM] = covariance[k];
          matrix[j + i * DIM] = covariance[k];
        }
      }

      return matrix;
    }

  private:

    double sum_pixels_;
    double sum_pixels_sq_;
    coord_type sum_pixels_coords_;
    coord_type sum_pixels_delta_sq_;
    coord_type sum_pixels_delta_cross_;
  };
}}

#endif /* DIALS_ALGORITHMS_IMAGE_CENTROID_CENTROID_LIST_H */
