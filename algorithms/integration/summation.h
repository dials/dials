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

#include <algorithm>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::int3;
  using scitbx::af::int6;
  using scitbx::af::sqrt;

  /**
   * Class to sum the intensity in 3D
   */
  template <typename FloatType = double>
  class Summation {
  public:

    /**
     * Perform the summation integration
     * @param signal The signal array
     * @param background The background array
     * @param n_background The number of counts used to calculate background
     */
    Summation(const af::const_ref<FloatType> &signal,
              const af::const_ref<FloatType> &background,
              std::size_t n_background) {
      init(signal, background, n_background);
    }

    /**
     * Perform the summation integration
     * @param signal The signal array
     * @param background The background array
     * @param mask The mask array
     * @param n_background The number of counts used to calculate background
     */
    Summation(const af::const_ref<FloatType> &signal,
              const af::const_ref<FloatType> &background,
              const af::const_ref<bool> &mask,
              std::size_t n_background) {
      init(signal, background, mask, n_background);
    }

    /**
     * Perform the summation integration
     * @param signal The signal array
     * @param background The background array
     * @param n_background The number of counts used to calculate background
     */
    Summation(const af::const_ref< FloatType, af::c_grid<2> > &signal,
              const af::const_ref< FloatType, af::c_grid<2> > &background,
              std::size_t n_background) {
      init(signal.as_1d(), background.as_1d(), n_background);
    }

    /**
     * Perform the summation integration
     * @param signal The signal array
     * @param background The background array
     * @param mask The mask array
     * @param n_background The number of counts used to calculate background
     */
    Summation(const af::const_ref< FloatType, af::c_grid<2> > &signal,
              const af::const_ref< FloatType, af::c_grid<2> > &background,
              const af::const_ref< bool, af::c_grid<2> > &mask,
              std::size_t n_background) {
      init(signal.as_1d(), background.as_1d(), mask.as_1d(), n_background);
    }

    /**
     * Perform the summation integration
     * @param signal The signal array
     * @param background The background array
     * @param n_background The number of counts used to calculate background
     */
    Summation(const af::const_ref< FloatType, af::c_grid<3> > &signal,
              const af::const_ref< FloatType, af::c_grid<3> > &background,
              std::size_t n_background) {
      init(signal.as_1d(), background.as_1d(), n_background);
    }

    /**
     * Perform the summation integration
     * @param signal The signal array
     * @param background The background array
     * @param mask The mask array
     * @param n_background The number of counts used to calculate background
     */
    Summation(const af::const_ref< FloatType, af::c_grid<3> > &signal,
              const af::const_ref< FloatType, af::c_grid<3> > &background,
              const af::const_ref< bool, af::c_grid<3> > &mask,
              std::size_t n_background) {
      init(signal.as_1d(), background.as_1d(), mask.as_1d(), n_background);
    }

    /**
     * @returns The reflection intensity
     */
    FloatType intensity() const {
      return signal_intensity() - background_intensity();
    }

    /**
     * @returns the variance on the integrated intensity
     */
    FloatType variance() const {
      return signal_variance() + background_variance();
    }

    /**
     * @returns the standard deviation on the intensity
     */
    FloatType standard_deviation() const {
      return std::sqrt(variance());
    }

    /**
     * @returns the signal intensity
     */
    FloatType signal_intensity() const {
      return signal_intensity_;
    }

    /**
     * @returns the variance on the signal intensity
     */
    FloatType signal_variance() const {
      return signal_variance_;
    }

    /**
     * @returns the standard deviation on the signal intensity
     */
    FloatType signal_standard_deviation() const {
      return std::sqrt(signal_variance());
    }

    /**
     * @returns the background intensity
     */
    FloatType background_intensity() const {
      return background_intensity_;
    }

    /**
     * @returns the variance on the background intensity
     */
    FloatType background_variance() const {
      double m_n = n_background_ > 0 ?
        (double)n_signal_ / (double)n_background_ : 0.0;
      return background_variance_ * m_n;
    }

    /**
     * @returns the standard deviation on the background intensity
     */
    FloatType background_standard_deviation() const {
      return std::sqrt(background_variance());
    }

    /**
     * @returns the number of background pixels
     */
    std::size_t n_background() const {
      return n_background_;
    }

    /**
     * @returns the number of signal pixels
     */
    std::size_t n_signal() const {
      return n_signal_;
    }

  private:

    /**
     * Integrate the intensity
     * @param signal The signal to integrate
     * @param background The background to the signal
     * @param n_background The number of counts used to calculate background
     */
    void init(const af::const_ref<FloatType> &signal,
              const af::const_ref<FloatType> &background,
              std::size_t n_background)
    {
      // Check both arrays are the same size
      DIALS_ASSERT(signal.size() == background.size());

      // Save the number of background pixels
      n_background_ = n_background;

      // Calculate the signal and background intensity
      n_signal_ = signal.size();
      signal_intensity_ = 0.0;
      background_intensity_ = 0.0;
      for (std::size_t i = 0; i < signal.size(); ++i) {
        signal_intensity_ += signal[i];
        background_intensity_ += background[i];
      }

      // Set the signal and background variance
      signal_variance_ = std::abs(signal_intensity_);
      background_variance_ = std::abs(background_intensity_);
    }

    /**
     * Integrate the intensity
     * @param signal The signal to integrate
     * @param background The background to the signal
     * @param mask The mask to the signal
     * @param n_background The number of counts used to calculate background
     */
    void init(const af::const_ref<FloatType> &signal,
              const af::const_ref<FloatType> &background,
              const af::const_ref<bool> &mask,
              std::size_t n_background)
    {
      // Check both arrays are the same size
      DIALS_ASSERT(signal.size() == background.size());
      DIALS_ASSERT(signal.size() == mask.size());

      // Save the number of background pixels
      n_background_ = n_background;

      // Calculate the signal and background intensity
      n_signal_ = 0;
      signal_intensity_ = 0.0;
      background_intensity_ = 0.0;
      for (std::size_t i = 0; i < signal.size(); ++i) {
        if (mask[i]) {
          signal_intensity_ += signal[i];
          background_intensity_ += background[i];
          n_signal_++;
        }
      }

      // Set the signal and background variance
      signal_variance_ = std::abs(signal_intensity_);
      background_variance_ = std::abs(background_intensity_);
    }

    FloatType signal_intensity_;
    FloatType signal_variance_;
    FloatType background_intensity_;
    FloatType background_variance_;
    std::size_t n_background_;
    std::size_t n_signal_;
  };

}}

#endif /* DIALS_ALGORITHMS_INTEGRATION_SUMMATION_H */
