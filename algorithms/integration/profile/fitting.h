/*
 * fitting.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_PROFILE_FITTING_H
#define DIALS_ALGORITHMS_INTEGRATION_PROFILE_FITTING_H

#include <boost/math/tools/minima.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::sum;
  using scitbx::af::flex_double;
  using boost::math::tools::brent_find_minima;

  /**
   * A class representing the profile model to minimize.
   */
  class ProfileModel {
  public:

    /**
     * Instantiate the model with the reflection profile
     * @param p The profile to fit to
     * @param c The contents of the pixels
     * @param b The background of the pixels
     */
    ProfileModel(const flex_double &p,
                 const flex_double &c,
                 const flex_double &b)
      : p_(p), c_(c), b_(b) {
      DIALS_ASSERT(p_.size() == c_.size());
      DIALS_ASSERT(p_.size() == b_.size());
    }

    /**
     * Evaluate the target function to minimize.
     * @param I the intensity.
     * @returns The value of the target function
     */
    double operator()(double I) const {
      double phi_I = 0.0;
      for (std::size_t j = 0; j < p_.size(); ++j) {
        double d = (c_[j] - b_[j] - I * p_[j]);
        double v = b_[j] + I * p_[j];
        if (v == 0.0) {
          continue;
        }
        phi_I += d * d / v;
      }
      return phi_I;
    }

    /**
     * Calculate the variance.
     * @param I the intensity
     * @returns The variance
     */
    double variance(double I) const {
      double var = 0.0;
      for (std::size_t j = 0; j < p_.size(); ++j) {
        var += b_[j] + I * p_[j];
      }
      return var;
    }

  private:
    flex_double p_, c_, b_;
  };


  /**
   * Class to fix the observed with the reference profile
   */
  class ProfileFitting {
  public:

    /**
     * Instantiate the fitting algorithm with the reflection profile
     * @param p The profile to fit to
     * @param c The contents of the pixels
     * @param b The background of the pixels
     */
    ProfileFitting(const flex_double &p,
                   const flex_double &c,
                   const flex_double &b,
                   int bits = 16,
                   std::size_t max_iter = 50)
      : model_(p, c, b) {

      // Set the maximum number of iterations and the bounds
      double xmin = 0;
      double xmax = sum(c.const_ref());

      // Find the minimum of the function
      std::pair<double, double> xf = brent_find_minima(
        model_, xmin, xmax, bits, max_iter);

      // Set the intensity
      intensity_ = xf.first;

      // Ensure the intensity is > 0 and less than max
      DIALS_ASSERT(intensity_ > xmin);
    }

    /**
     * @returns The intensity
     */
    double intensity() const {
      return intensity_;
    }

    /**
     * @returns the variance
     */
    double variance() const {
      return model_.variance(intensity_);
    }

  private:
    ProfileModel model_;
    double intensity_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_FITTING_H */
