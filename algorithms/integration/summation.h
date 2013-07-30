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
#include <dials/algorithms/shoebox/mask_code.h>

namespace dials { namespace algorithms {

  using scitbx::af::int6;
  using scitbx::af::sqrt;
  using scitbx::af::flex_bool;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  typedef flex< vec3<double> >::type flex_vec3_double;

  /**
   * Class to sum the intensity in 3D
   */
  class SumIntensity3d {
  public:

    /**
     * Integrate the intensity
     * @param signal The signal to integrate
     * @param background The background to the signal
     */
    SumIntensity3d(const flex_double &signal,
                   const flex_double &background)
    {
      // Check both arrays are the same size
      DIALS_ASSERT(signal.size() == background.size());

      // Calculate the signal and background intensity
      signal_intensity_ = 0.0;
      background_intensity_ = 0.0;
      for (std::size_t i = 0; i < signal.size(); ++i) {
        signal_intensity_ += signal[i];
        background_intensity_ += background[i];
      }

      // Set the signal and background variance
      signal_variance_ = signal_intensity_;
      background_variance_ = background_intensity_;
    }

    /**
     * Integrate the intensity
     * @param signal The signal to integrate
     * @param background The background to the signal
     * @param mask The mask to the signal
     */
    SumIntensity3d(const flex_double &signal,
                   const flex_double &background,
                   const flex_bool &mask)
    {
      // Check both arrays are the same size
      DIALS_ASSERT(signal.size() == background.size());

      // Calculate the signal and background intensity
      signal_intensity_ = 0.0;
      background_intensity_ = 0.0;
      for (std::size_t i = 0; i < signal.size(); ++i) {
        if (mask[i]) {
          signal_intensity_ += signal[i];
          background_intensity_ += background[i];
        }
      }

      // Set the signal and background variance
      signal_variance_ = signal_intensity_;
      background_variance_ = background_intensity_;
    }

    /**
     * @returns The reflection intensity
     */
    double intensity() const {
      return signal_intensity() - background_intensity();
    }

    /**
     * @returns the signal intensity
     */
    double signal_intensity() const {
      return signal_intensity_;
    }

    /**
     * @returns the background intensity
     */
    double background_intensity() const {
      return background_intensity_;
    }

    /**
     * @returns the variance on the integrated intensity
     */
    double variance() const {
      return signal_variance() + background_variance();
    }

    /**
     * @returns the variance on the signal intensity
     */
    double signal_variance() const {
      return signal_variance_;
    }

    /**
     * @returns the variance on the background intensity
     */
    double background_variance() const {
      return background_variance_;
    }

  private:

    double signal_intensity_;
    double signal_variance_;
    double background_intensity_;
    double background_variance_;
  };

  /**
   * Class to perform summation integration.
   */
  class IntegrateBySummation {
  public:

    /**
     * Perform the integration on a 3D image.
     * @param pixels The 3D image.
     * @param background The background
     */
    IntegrateBySummation(const flex_double &pixels,
                         const flex_double &background) {

      // Create an array with pixels - background
      flex_double pixels_m_background(subtract_background(pixels, background));

      // Calculate the itensity and sigma
      SumIntensity3d isum(pixels, background);
      intensity_ = isum.intensity();
      ivariance_ = isum.variance();
    }

    /**
     * Perform the integration on a 3D image with a mask.
     * @param pixels The 3D image.
     * @param background The pixel background
     * @param mask The corresponding mask
     */
    IntegrateBySummation(const flex_double &pixels,
                         const flex_double &background,
                         const flex_bool &mask) {

      // Create an array with pixels - background
      flex_double pixels_m_background(subtract_background(pixels, background));

      // Calculate the itensity and sigma
      SumIntensity3d isum(pixels, background, mask);
      intensity_ = isum.intensity();
      ivariance_ = isum.variance();
    }

    /** @return The integrated intensity. */
    double intensity() const {
      return intensity_;
    }

    /** @return The variance on the intensity */
    double variance() const {
      return ivariance_;
    }

    /** @return the standard deviation on the intensity */
    double standard_deviation() const {
      return std::sqrt(variance());
    }

  private:

    flex_double subtract_background(const flex_double &pixels,
        const flex_double &background) {
      // Check the sizes are the same
      DIALS_ASSERT(pixels.size() == background.size());

      // Create an array with pixels - background
      flex_double pixels_m_background(pixels.accessor());
      for (std::size_t i = 0; i < pixels.size(); ++i) {
        pixels_m_background[i] = pixels[i] - background[i];
      }

      // Return the subtracted pixels
      return pixels_m_background;
    }

    double intensity_;
    double ivariance_;
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
                          const flex_bool &mask) const {
      return integrator(pixels, background, mask);
    }

    /**
     * Integrate a reflection
     * @param r The reflection container
     */
    void operator()(Reflection &r) const {

      flex_int shoebox_mask = r.get_shoebox_mask();
      flex_bool mask(shoebox_mask.accessor());
      for (std::size_t i = 0; i < mask.size(); ++i) {
        mask[i] = (shoebox_mask[i] & shoebox::Valid) ? true : false;
      }

      // Integrate the reflection
      integrator result = this->operator()(
        r.get_shoebox(),
        r.get_shoebox_background(),
        mask);

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
          if (reflections[i].is_valid()) {
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
