/*
 * summation.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_SUMMATION_H
#define DIALS_ALGORITHMS_INTEGRATION_SUMMATION_H

#include <scitbx/array_family/tiny_algebra.h>
#include <dials/algorithms/image/centroid/centroid_image.h>
#include <dials/algorithms/image/centroid/centroid_masked_image.h>
#include <dials/algorithms/image/centroid/centroid_points.h>

namespace dials { namespace algorithms {

  using scitbx::af::sqrt;

  typedef flex< vec3<double> >::type flex_vec3_double;

  /**
   * Class to perform summation integration.
   */
  class IntegrateBySummation {
  public:

    // Useful typedefs
    typedef vec3<double> coord_type;
    typedef tiny<double, 9> matrix_type;

    /**
     * Perform the integration on a 3D image.
     * @param pixels The 3D image.
     */
    IntegrateBySummation(const flex_double &pixels) {
      CentroidImage3d centroid(pixels);
      intensity_         = centroid.sum_pixels();
      centroid_          = centroid.mean();
      variance_          = centroid.unbiased_variance();
      standard_error_sq_ = centroid.unbiased_standard_error_sq();
      covariance_matrix_ = centroid.covariance_matrix();
    }

    /**
     * Perform the integration on a 3D image with a mask.
     * @param pixels The 3D image.
     * @param mask The corresponding mask
     */
    IntegrateBySummation(const flex_double &pixels, const flex_int &mask) {
      CentroidMaskedImage3d centroid(pixels, mask);
      intensity_         = centroid.sum_pixels();
      centroid_          = centroid.mean();
      variance_          = centroid.unbiased_variance();
      standard_error_sq_ = centroid.unbiased_standard_error_sq();
      covariance_matrix_ = centroid.covariance_matrix();
    }

    /**
     * Perform the integration on a set of 3D points
     * @param pixels The image pixels.
     * @param coords The image coordinates.
     */
    IntegrateBySummation(const flex_double &pixels,
        const flex_vec3_double &coords) {
      CentroidPoints< vec3<double> > centroid(pixels, coords);
      intensity_         = centroid.sum_pixels();
      centroid_          = centroid.mean();
      variance_          = centroid.unbiased_variance();
      standard_error_sq_ = centroid.unbiased_standard_error_sq();
      covariance_matrix_ = centroid.covariance_matrix();
    }

    /** @return The integrated intensity. */
    double intensity() const {
      return intensity_;
    }

    /** @return The spot centroid. */
    coord_type centroid() const {
      return centroid_;
    }

    /** @return The centroid variance. */
    coord_type variance() const {
      return variance_;
    }

    /** @return The centroid standard error squared. */
    coord_type standard_error_sq() const {
      return standard_error_sq_;
    }

    /** @return The centroid standard error. */
    coord_type standard_error() const {
      return sqrt(standard_error_sq_.as_tiny());
    }

    /** @return The covariance matrix. */
    matrix_type covariance_matrix() const {
      return covariance_matrix_;
    }

  private:

    double intensity_;
    coord_type centroid_;
    coord_type variance_;
    coord_type standard_error_sq_;
    matrix_type covariance_matrix_;
  };

}}

#endif /* DIALS_ALGORITHMS_INTEGRATION_SUMMATION_H */
