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
      return sum_p_ - sum_b_;
    }

    /**
     * @returns the variance on the integrated intensity
     */
    FloatType variance() const {
      FloatType Is = intensity();
      FloatType Ib = sum_b_;
      FloatType m_n = n_background_ > 0 ?
        (FloatType)n_signal_ / (FloatType)n_background_ : 0.0;
      return std::abs(Is) + std::abs(Ib) * (1.0 + m_n);
    }

    /**
     * @returns the standard deviation on the intensity
     */
    FloatType standard_deviation() const {
      return std::sqrt(variance());
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
      sum_p_ = 0.0;
      sum_b_ = 0.0;
      for (std::size_t i = 0; i < signal.size(); ++i) {
        sum_p_ += signal[i];
        sum_b_ += background[i];
      }
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
      sum_p_ = 0.0;
      sum_b_ = 0.0;
      for (std::size_t i = 0; i < signal.size(); ++i) {
        if (mask[i]) {
          sum_p_ += signal[i];
          sum_b_ += background[i];
          n_signal_++;
        }
      }
    }

    FloatType sum_p_;
    FloatType sum_b_;
    std::size_t n_background_;
    std::size_t n_signal_;
  };

}}

#endif /* DIALS_ALGORITHMS_INTEGRATION_SUMMATION_H */
