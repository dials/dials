/*
 * bayesian_integrator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_BAYESIAN_INTEGRATOR_H
#define DIALS_ALGORITHMS_INTEGRATION_BAYESIAN_INTEGRATOR_H

#include <algorithm>
#include <boost/math/special_functions/gamma.hpp>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <dials/model/data/mask_code.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using dials::model::Background;
  using dials::model::BackgroundUsed;
  using dials::model::Foreground;
  using dials::model::Valid;
  using scitbx::af::int3;
  using scitbx::af::int6;
  using scitbx::af::sqrt;

  /**
   * Class to sum the intensity in 3D
   */
  template <typename FloatType = double>
  class BayesianIntegrator {
  public:
    /**
     * Perform the summation integration
     * @param signal The signal array
     * @param background The background array
     * @param mask The mask array
     */
    BayesianIntegrator(const af::const_ref<FloatType> &signal,
                       const af::const_ref<FloatType> &background,
                       const af::const_ref<int> &mask) {
      init(signal, background, mask);
    }

    /**
     * Perform the summation integration
     * @param signal The signal array
     * @param background The background array
     * @param mask The mask array
     */
    BayesianIntegrator(const af::const_ref<FloatType, af::c_grid<2> > &signal,
                       const af::const_ref<FloatType, af::c_grid<2> > &background,
                       const af::const_ref<int, af::c_grid<2> > &mask) {
      init(signal.as_1d(), background.as_1d(), mask.as_1d());
    }

    /**
     * Perform the summation integration
     * @param signal The signal array
     * @param background The background array
     * @param mask The mask array
     */
    BayesianIntegrator(const af::const_ref<FloatType, af::c_grid<3> > &signal,
                       const af::const_ref<FloatType, af::c_grid<3> > &background,
                       const af::const_ref<int, af::c_grid<3> > &mask) {
      init(signal.as_1d(), background.as_1d(), mask.as_1d());
    }

    /**
     * @returns The reflection intensity
     */
    FloatType intensity() const {
      using boost::math::gamma_q;
      double C = sum_p_;
      double B = sum_b_;
      return gamma_q(C + 1, B) * ((C + 1) * gamma_q(C + 1, B) - B * gamma_q(C + 1, B));
    }

    /**
     * @returns the variance on the integrated intensity
     */
    FloatType variance() const {
      FloatType Is = intensity();
      FloatType Ib = sum_b_;
      FloatType m_n =
        n_background_ > 0 ? (FloatType)n_signal_ / (FloatType)n_background_ : 0.0;
      return std::abs(Is) + std::abs(Ib) * (1.0 + m_n);
    }

    /**
     * @returns The background
     */
    FloatType background() const {
      return sum_b_;
    }

    /**
     * @returns The background variance
     */
    FloatType background_variance() const {
      FloatType m_n =
        n_background_ > 0 ? (FloatType)n_signal_ / (FloatType)n_background_ : 0.0;
      return std::abs(sum_b_) * (1.0 + m_n);
    }

    /**
     * @returns the number of signal pixels
     */
    std::size_t n_signal() const {
      return n_signal_;
    }

    /**
     * @returns the number of background pixels
     */
    std::size_t n_background() const {
      return n_background_;
    }

    /**
     * @returns Was the algorithm successful
     */
    bool success() const {
      return success_;
    }

  private:
    /**
     * Integrate the intensity
     * @param signal The signal to integrate
     * @param background The background to the signal
     * @param mask The mask to the signal
     */
    void init(const af::const_ref<FloatType> &signal,
              const af::const_ref<FloatType> &background,
              const af::const_ref<int> &mask) {
      // Check both arrays are the same size
      DIALS_ASSERT(signal.size() == background.size());
      DIALS_ASSERT(signal.size() == mask.size());

      // Calculate the signal and background intensity
      int bg_code = Valid | Background | BackgroundUsed;
      success_ = true;
      sum_p_ = 0;
      sum_b_ = 0;
      n_signal_ = 0;
      n_background_ = 0;
      for (std::size_t i = 0; i < signal.size(); ++i) {
        if ((mask[i] & Foreground) == Foreground) {
          if ((mask[i] & Valid) == Valid) {
            sum_p_ += signal[i];
            sum_b_ += background[i];
            n_signal_++;
          } else {
            success_ = false;
          }
        } else if ((mask[i] & bg_code) == bg_code) {
          n_background_++;
        }
      }
    }

    FloatType sum_p_;
    FloatType sum_b_;
    std::size_t n_background_;
    std::size_t n_signal_;
    bool success_;
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_BAYESIAN_INTEGRATOR_H */
