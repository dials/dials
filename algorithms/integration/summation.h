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

  using scitbx::af::int3;
  using scitbx::af::int6;
  using scitbx::af::sqrt;
  using dials::model::Reflection;

  /**
   * Class to sum the intensity in 3D
   */
  class Summation {
  public:

    /**
     * Integrate the intensity
     * @param signal The signal to integrate
     * @param background The background to the signal
     */
    Summation(const af::const_ref<double> &signal,
               const af::const_ref<double> &background)
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
    Summation(const af::const_ref<double> &signal,
               const af::const_ref<double> &background,
               const af::const_ref<bool> &mask)
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
     * @returns the variance on the integrated intensity
     */
    double variance() const {
      return signal_variance() + background_variance();
    }

    /**
     * @returns the standard deviation on the intensity
     */
    double standard_deviation() const {
      return std::sqrt(variance());
    }

    /**
     * @returns the signal intensity
     */
    double signal_intensity() const {
      return signal_intensity_;
    }

    /**
     * @returns the variance on the signal intensity
     */
    double signal_variance() const {
      return signal_variance_;
    }

    /**
     * @returns the standard deviation on the signal intensity
     */
    double signal_standard_deviation() const {
      return std::sqrt(signal_variance());
    }

    /**
     * @returns the background intensity
     */
    double background_intensity() const {
      return background_intensity_;
    }

    /**
     * @returns the variance on the background intensity
     */
    double background_variance() const {
      return background_variance_;
    }

    /**
     * @returns the standard deviation on the background intensity
     */
    double background_standard_deviation() const {
      return std::sqrt(background_variance());
    }

  private:

    double signal_intensity_;
    double signal_variance_;
    double background_intensity_;
    double background_variance_;
  };


  /**
   * A class to do 3D summation integration
   */
  class Summation3d {
  public:

    typedef Summation integrator;

    /** Init the algorithm. */
    Summation3d() {}

    /**
     * Integrate a set of pixels with a mask
     * @param pixels The array of pixels
     * @param background The background pixels
     * @param mask The mask
     * @returns The integrator struct
     */
    integrator operator()(
        const af::const_ref< double, af::c_grid<3> > &pixels,
        const af::const_ref< double, af::c_grid<3> > &background,
        const af::const_ref< bool, af::c_grid<3> > &mask) const {
      return integrator(pixels.as_1d(), background.as_1d(), mask.as_1d());
    }

    /**
     * Integrate a reflection
     * @param r The reflection container
     */
    void operator()(Reflection &r) const {

      af::const_ref< int, af::c_grid<3> > shoebox_mask =
        r.get_shoebox_mask().const_ref();
      af::versa< bool, af::c_grid<3> > mask(shoebox_mask.accessor());
      for (std::size_t i = 0; i < mask.size(); ++i) {
        mask[i] = (shoebox_mask[i] & shoebox::Valid) ? true : false;
      }

      // Integrate the reflection
      integrator result = this->operator()(
        r.get_shoebox().const_ref(),
        r.get_shoebox_background().const_ref(),
        mask.const_ref());

      r.set_intensity(result.intensity());
      r.set_intensity_variance(result.variance());
    }

    /**
     * Integrate a list of reflections
     * @param reflections The reflection list
     */
    void operator()(af::ref<Reflection> reflections) const {
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
