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

#include <omptbx/omp_or_stubs.h>
#include <algorithm>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <dials/model/data/reflection.h>
#include <dials/algorithms/image/centroid/centroid_image.h>
#include <dials/algorithms/image/centroid/centroid_masked_image.h>
#include <dials/algorithms/image/centroid/centroid_points.h>

namespace dials { namespace algorithms {

  using scitbx::af::int6;
  using scitbx::af::sqrt;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

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

    /** @return The variance on the intensity */
    double variance() const {
      return 0.0;
    }

    /** @return the standard deviation on the intensity */
    double standard_deviation() const {
      return std::sqrt(variance());
    }

    /** @return The spot centroid. */
    coord_type centroid() const {
      return centroid_;
    }

    /** @return The centroid variance. */
    coord_type centroid_variance() const {
      return variance_;
    }

    /** @return The centroid standard error squared. */
    coord_type centroid_standard_error_sq() const {
      return standard_error_sq_;
    }

    /** @return The centroid standard error. */
    coord_type centroid_standard_error() const {
      return sqrt(standard_error_sq_.as_tiny());
    }

    /** @return The covariance matrix. */
    matrix_type centroid_covariance_matrix() const {
      return covariance_matrix_;
    }

  private:

    double intensity_;
    coord_type centroid_;
    coord_type variance_;
    coord_type standard_error_sq_;
    matrix_type covariance_matrix_;
  };


  /**
   * A class to do 3D summation integration
   */
  class Summation3d {
  public:

    typedef IntegrateBySummation integrator;

    /** Init the algorithm. */
    Summation3d() {}

    /**
     * Integrate a set of pixels with a mask
     * @param pixels The array of pixels
     * @param background The background pixels
     * @param mask The mask
     * @returns The integrator struct
     */
    integrator operator()(const flex_double &pixels,
                          const flex_double &background,
        const flex_int &mask) const {
      flex_double pixels_double(pixels.accessor());
      for (std::size_t i = 0; i < pixels.size(); ++i) {
        pixels_double[i] = (double)(pixels[i] - background[i]);
      }
      return integrator(pixels_double, mask);
    }

    /**
     * Integrate a reflection
     * @param r The reflection container
     */
    void operator()(Reflection &r) const {

      // Integrate the reflection
      integrator result = this->operator()(
        r.get_shoebox(),
        r.get_shoebox_background(),
        r.get_shoebox_mask());

      // Get the centroid offset
      int6 bbox = r.get_bounding_box();
      vec3<double> offset(bbox[0], bbox[2], bbox[4]);

      // Put data back into reflection container
      r.set_centroid_position(offset + result.centroid());
      r.set_centroid_variance(result.centroid_standard_error_sq());
      r.set_centroid_sq_width(result.centroid_variance());
      r.set_intensity(result.intensity());
      r.set_intensity_variance(result.variance());
    }

    /**
     * Integrate a list of reflections
     * @param reflections The reflection list
     */
    void operator()(ReflectionList &reflections) const {
      #pragma omp parallel for
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        try {
          if (reflections[i].get_status() == 0) {
            this->operator()(reflections[i]);
          }
        } catch (dials::error) {
          reflections[i].set_valid(false);
        }
      }
    }
  };

}}

#endif /* DIALS_ALGORITHMS_INTEGRATION_SUMMATION_H */
