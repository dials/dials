#ifndef DIALS_ALGORITHMS_INTEGRATION_SUMMATION_H
#define DIALS_ALGORITHMS_INTEGRATION_SUMMATION_H

#include <scitbx/array_family/tiny_algebra.h>
#include <dials/algorithms/image/centroid/centroid_image.h>
#include <dials/algorithms/image/centroid/centroid_masked_image.h>
#include <dials/algorithms/image/centroid/centroid_points.h>

namespace dials { namespace algorithms {

  using scitbx::af::sqrt;

  typedef flex< vec3<double> >::type flex_vec3_double;

  class IntegrateBySummation {
  public:

    typedef vec3<double> coord_type;
    typedef tiny<double, 9> matrix_type;

    IntegrateBySummation(const flex_double &pixels) {
      CentroidImage3d centroid(pixels);
      intensity_         = centroid.sum_pixels();
      centroid_          = centroid.mean();
      variance_          = centroid.unbiased_variance();
      standard_error_sq_ = centroid.unbiased_standard_error_sq();
      covariance_matrix_ = centroid.covariance_matrix();
    }

    IntegrateBySummation(const flex_double &pixels, const flex_int &mask) {
      CentroidMaskedImage3d centroid(pixels, mask);
      intensity_         = centroid.sum_pixels();
      centroid_          = centroid.mean();
      variance_          = centroid.unbiased_variance();
      standard_error_sq_ = centroid.unbiased_standard_error_sq();
      covariance_matrix_ = centroid.covariance_matrix();
    }

    IntegrateBySummation(const flex_double &pixels,
        const flex_vec3_double &coords) {
      CentroidPoints< vec3<double> > centroid(pixels, coords);
      intensity_         = centroid.sum_pixels();
      centroid_          = centroid.mean();
      variance_          = centroid.unbiased_variance();
      standard_error_sq_ = centroid.unbiased_standard_error_sq();
      covariance_matrix_ = centroid.covariance_matrix();
    }

    double intensity() const {
      return intensity_;
    }

    coord_type centroid() const {
      return centroid_;
    }

    coord_type variance() const {
      return variance_;
    }

    coord_type standard_error_sq() const {
      return standard_error_sq_;
    }

    coord_type standard_error() const {
      return sqrt(standard_error_sq_.as_tiny());
    }

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
