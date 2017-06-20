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
#ifndef DIALS_ALGORITHMS_INTEGRATION_MAXIMUM_LIKELIHOOD_FITTING_H
#define DIALS_ALGORITHMS_INTEGRATION_MAXIMUM_LIKELIHOOD_FITTING_H

#include <algorithm>
#include <vector>
#include <scitbx/vec2.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <dials/model/data/mask_code.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::af::sum;
  using scitbx::vec2;

  /**
   * Class to fit the observed with the reference profile
   */
  template <typename FloatType = double>
  class MLProfileFitting {
  public:

    typedef FloatType float_type;

    /**
     * Instantiate the fitting algorithm with the reflection profile
     * @param p The profile to fit to
     * @param c The contents of the pixels
     * @param b The background of the pixels
     */
    MLProfileFitting(const af::const_ref<FloatType, af::c_grid<3> > &s,
                   const af::const_ref<bool, af::c_grid<3> > &m,
                   const af::const_ref<FloatType, af::c_grid<3> > &c,
                   const af::const_ref<FloatType, af::c_grid<3> > &b,
                   double eps = 1e-3,
                   std::size_t max_iter = 10)
    {
      // Check the input
      DIALS_ASSERT(s.size() == m.size());
      DIALS_ASSERT(s.size() == c.size());
      DIALS_ASSERT(s.size() == b.size());
      DIALS_ASSERT(eps > 0.0);
      DIALS_ASSERT(max_iter >= 1);

      // Normalize input array
      std::vector <FloatType> bb(b.size());
      std::vector <FloatType> ss(s.size());
      double sum_b = af::sum(b);
      double sum_s = af::sum(s);
      DIALS_ASSERT(sum_b > 0 && sum_s > 0);
      for (std::size_t i = 0; i < s.size(); ++i) {
        bb[i] = b[i] / sum_b;
        if (sum_s > 1) {
          ss[i] = s[i] / sum_s;
        } else {
          ss[i] = s[i];
        }
      }
      background_ = sum_b;

      // Iterate to calculate the intensity. Exit if intensity goes less
      // than zero or if the tolerance or number of iteration is reached.
      estimate_intensity(
          af::const_ref<FloatType, af::c_grid<3> >(&ss[0], s.accessor()),
          m,
          c,
          af::const_ref<FloatType, af::c_grid<3> >(&bb[0], b.accessor()),
          eps,
          max_iter);

      // Set the intensity and variance
      correlation_ = compute_correlation(
          af::const_ref<FloatType, af::c_grid<3> >(&ss[0], s.accessor()),
          m,
          c,
          af::const_ref<FloatType, af::c_grid<3> >(&bb[0], b.accessor()));
    }

    /**
     * @returns The intensity
     */
    double intensity() const {
      return intensity_;
    }

    /**
     * @returns The total background counts
     */
    double background() const {
      return background_;
    }

    /**
     * @returns the variance
     */
    double variance() const {
      return variance_;
    }

    /**
     * @returns the correlation
     */
    double correlation() const {
      return correlation_;
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
    void estimate_intensity(const af::const_ref<FloatType, af::c_grid<3> > &s,
                       const af::const_ref<bool, af::c_grid<3> > &m,
                       const af::const_ref<FloatType, af::c_grid<3> > &c,
                       const af::const_ref<FloatType, af::c_grid<3> > &b,
                       double eps,
                       std::size_t max_iter) {
      double TOL = 1e-7;
      double sum_c = 0;
      double sum_s = 0;
      double sum_b = 0;
      for (std::size_t i = 0; i < c.size(); ++i) {
        if (m[i]) {
          sum_c += c[i];
          sum_s += s[i];
          sum_b += b[i];
        }
      }
      DIALS_ASSERT(sum_s > 0.1);
      DIALS_ASSERT(sum_b > 0.1);
      DIALS_ASSERT(sum_s <= (1.0+TOL));
      DIALS_ASSERT(sum_b <= (1.0+TOL));
      double S = std::max(sum_c - background_, 1.0) / sum_s; //sum_c / 2;
      double B = background_ / sum_b;//sum_c / 2;
      double V = 0;
      for (niter_ = 0; niter_ < max_iter; ++niter_) {
        double sum1 = 0.0;
        double sum2 = 0.0;
        double sumv = 0.0;
        for (std::size_t i = 0; i < s.size(); ++i) {
          if (m[i]) {
            double v = B*b[i] + S*s[i];
            if (v > 0) {
              sum1 += b[i] * c[i] / v;
              sum2 += s[i] * c[i] / v;
              if (s[i] > 0) {
                sumv += v;
              }
            }
          }
        }
        DIALS_ASSERT(sum1 >= 0 && sum2 >= 0);
        double Bold = B;
        double Sold = S;
        B = B * sum1 / sum_b;
        S = S * sum2 / sum_s;
        V = sumv;
        if ((Bold-B)*(Bold-B) + (Sold-S)*(Sold-S) < eps*eps) {
          break;
        }
      }
      intensity_ = S;
      background_ = B;
      variance_ = V;
    }

    /**
     * Compute the correlation coefficient between the profile and reference
     */
    double
    compute_correlation(const af::const_ref<FloatType, af::c_grid<3> > &s,
                        const af::const_ref<bool, af::c_grid<3> > &m,
                        const af::const_ref<FloatType, af::c_grid<3> > &c,
                        const af::const_ref<FloatType, af::c_grid<3> > &b) const {
      double xb = 0.0, yb = 0.0;
      std::size_t count = 0;
      for (std::size_t i = 0; i < s.size(); ++i) {
        if (m[i]) {
          xb += intensity_*s[i] + background_*b[i];
          yb += c[i];
          count++;
        }
      }
      DIALS_ASSERT(count > 0);
      xb /= count;
      yb /= count;
      double sdxdy = 0.0, sdx2 = 0.0, sdy2 = 0.0;
      for (std::size_t i = 0; i < s.size(); ++i) {
        if (m[i] && s[i] > 0) {
          double dx = (intensity_*s[i] + background_*b[i]) - xb;
          double dy = c[i] - yb;
          sdxdy += dx*dy;
          sdx2 += dx*dx;
          sdy2 += dy*dy;
        }
      }
      double result = 0.0;
      if (sdx2 > 0.0 && sdy2 > 0.0) {
        result = sdxdy / (std::sqrt(sdx2) * std::sqrt(sdy2));
      }
      return result;
    }

    double intensity_;
    double variance_;
    double background_;
    double correlation_;
    std::size_t niter_;
    double error_;
  };

}}

#endif /* DIALS_ALGORITHMS_INTEGRATION_MAXIMUM_LIKELIHOOD_FITTING_H */
