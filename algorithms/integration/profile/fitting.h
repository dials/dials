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
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/vec2.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::sum;
  using scitbx::vec2;

  /**
   * Class to fix the observed with the reference profile
   */
  template <typename FloatType = double>
  class ProfileFitting {
  public:

    typedef FloatType float_type;

    /**
     * Instantiate the fitting algorithm with the reflection profile
     * @param p The profile to fit to
     * @param c The contents of the pixels
     * @param b The background of the pixels
     */
    ProfileFitting(const af::const_ref<FloatType, af::c_grid<3> > &p,
                   const af::const_ref<FloatType, af::c_grid<3> > &c,
                   const af::const_ref<FloatType, af::c_grid<3> > &b,
                   double eps = 1e-3,
                   std::size_t max_iter = 10)
    {
      // Check the input
      DIALS_ASSERT(p.size() == c.size());
      DIALS_ASSERT(p.size() == b.size());
      DIALS_ASSERT(eps > 0.0);
      DIALS_ASSERT(max_iter >= 1);

      // Iterate to calculate the intensity. Exit if intensity goes less
      // than zero or if the tolerance or number of iteration is reached.
      double I0 = sum(c);
      vec2<double> I(0.0, 0.0);
      for (niter_ = 0; niter_ < max_iter; ++niter_) {
        I = estimate_intensity(p, c, b, I0);
        DIALS_ASSERT(I[0] >= 0.0);
        if ((error_ = std::abs(I[0] - I0)) < eps) {
          break;
        }
        I0 = I[0];
      }

      // Set the intensity and variance
      intensity_ = I[0];
      variance_ = I[1];
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
      return variance_;
    }

    /**
     * @returns The number of iterations
     */
    std::size_t niter() const {
      return niter_;
    }

    /**
     * @returns The error in the fit
     */
    double error() const {
      return error_;
    }

  private:

    /**
     * Evaluate the next intensity iteration.
     * @ returns The estimate of the intensity
     */
    vec2<double>
    estimate_intensity(const af::const_ref<FloatType, af::c_grid<3> > &p,
                       const af::const_ref<FloatType, af::c_grid<3> > &c,
                       const af::const_ref<FloatType, af::c_grid<3> > &b,
                       double I) {
      double df = 0.0, d2f = 0.0, sum_v = 0.0;
      for (std::size_t i = 0; i < p.size(); ++i) {
        double v = b[i] + p[i] * I;
        double v2 = v*v;
        double v3 = v2*v;
        if (v2 > 0) {
          df  += p[i]*(v - c[i])*(v + c[i])/v2;
          d2f += 2.0*c[i]*c[i]*p[i]*p[i]/v3;
          sum_v += v;
        }
      }
      return vec2<double>(I - (d2f != 0 ? df / d2f : 0.0), sum_v);
    }

    double intensity_;
    double variance_;
    std::size_t niter_;
    double error_;
  };

  /**
   * Class to fix the observed with the reference profile
   */
  template <typename FloatType = double>
  class ProfileFitting2 {
  public:

    typedef FloatType float_type;

    /**
     * Instantiate the fitting algorithm with the reflection profile
     * @param p The profile to fit to
     * @param c The contents of the pixels
     * @param b The background of the pixels
     */
    ProfileFitting2(const af::const_ref<FloatType, af::c_grid<3> > &p,
                   const af::const_ref<FloatType, af::c_grid<3> > &c,
                   const af::const_ref<FloatType, af::c_grid<3> > &b,
                   double eps = 1e-3,
                   std::size_t max_iter = 10)
    {
      // Check the input
      DIALS_ASSERT(p.size() == c.size());
      DIALS_ASSERT(p.size() == b.size());
      DIALS_ASSERT(eps > 0.0);
      DIALS_ASSERT(max_iter >= 1);

      // Iterate to calculate the intensity. Exit if intensity goes less
      // than zero or if the tolerance or number of iteration is reached.
      double I = 0.0, I0 = sum(c);
      for (niter_ = 0; niter_ < max_iter; ++niter_) {
        I = estimate_intensity(p, c, b, I0);
        std::cout << I0 << ", " << I << std::endl;
        DIALS_ASSERT(I >= 0.0);
        if ((error_ = std::abs(I - I0)) < eps) {
          break;
        }
        I0 = I;
      }

      // Set the intensity and variance
      intensity_ = I;
      variance_ = estimate_variance(p, b, I);
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
      return variance_;
    }

    /**
     * @returns The number of iterations
     */
    std::size_t niter() const {
      return niter_;
    }

    /**
     * @returns The error in the fit
     */
    double error() const {
      return error_;
    }

  private:

    /**
     * Evaluate the next intensity iteration.
     * @ returns The estimate of the intensity
     */
    double estimate_intensity(const af::const_ref<FloatType, af::c_grid<3> > &p,
                              const af::const_ref<FloatType, af::c_grid<3> > &c,
                              const af::const_ref<FloatType, af::c_grid<3> > &b,
                              double I0) {
      double s1 = 0.0, s2 = 0.0;
      for (std::size_t i = 0; i < p.size(); ++i) {
        double v = b[i] + I0 * p[i];
        if (v == 0.0) continue;
        double pv = p[i] / v;
        s1 += (c[i] - b[i]) * pv;
        s2 += p[i] * pv;
      }
      return s2 == 0.0 ? 0.0 : s1 / s2;
    }

    /**
     * Calculate the total variance in the profile.
     * @returns The total variance
     */
    double estimate_variance(const af::const_ref<FloatType, af::c_grid<3> > &p,
                             const af::const_ref<FloatType, af::c_grid<3> > &b,
                             double I) {
      double V = 0.0;
      for (std::size_t i = 0; i < p.size(); ++i) {
        V += b[i] + I * p[i];
      }
      return V;
    }

    double intensity_;
    double variance_;
    std::size_t niter_;
    double error_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_FITTING_H */
