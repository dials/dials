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
#ifndef DIALS_ALGORITHMS_INTEGRATION_FIT_FITTING_H
#define DIALS_ALGORITHMS_INTEGRATION_FIT_FITTING_H

#include <algorithm>
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
                   const af::const_ref<bool, af::c_grid<3> > &m,
                   const af::const_ref<FloatType, af::c_grid<3> > &c,
                   const af::const_ref<FloatType, af::c_grid<3> > &b,
                   double eps = 1e-3,
                   std::size_t max_iter = 10)
    {
      // Check the input
      DIALS_ASSERT(p.size() == m.size());
      DIALS_ASSERT(p.size() == c.size());
      DIALS_ASSERT(p.size() == b.size());
      DIALS_ASSERT(eps > 0.0);
      DIALS_ASSERT(max_iter >= 1);

      double sumc = 0;
      double sumb = 0;
      double sums = 0;
      double minc = -1;
      double minI = 9999;
      for (std::size_t i = 0; i < m.size(); ++i) {
        DIALS_ASSERT(p[i] >= 0);
        DIALS_ASSERT(b[i] >= 0);
        DIALS_ASSERT(c[i] >= 0);
        if (m[i]) {
          sumc += c[i];
          sumb += b[i];
          sums += p[i];
          if (minc < 0 || c[i] < minc) minc = c[i];
          if (p[i] > 0 && b[i] / p[i] < minI) minI = b[i] / p[i];
        }
      }
      minI = -minI;
      DIALS_ASSERT(sumb >= 0);
      DIALS_ASSERT(sumc >= 0);
      DIALS_ASSERT(sums > 0);
      DIALS_ASSERT(minI <= 0);

      // Iterate to calculate the intensity. Exit if intensity goes less
      // than zero or if the tolerance or number of iteration is reached.
      double I0 = sumc - sumb;
      if (I0 < minI) {
        I0 = minI+1e-3;
      }
      vec2<double> I(0.0, 0.0);
      for (niter_ = 0; niter_ < max_iter; ++niter_) {
        double sum1 = 0.0;
        double sum2 = 0.0;
        double sumv = 0.0;
        for (std::size_t i = 0; i < p.size(); ++i) {
          if (m[i] && p[i] > 0) {
            double v = b[i] + I0 * p[i];
            DIALS_ASSERT(v > 0);
            sumv += v;
            if (v > 0) {
              sum1 += (c[i] - b[i]) * p[i] / v;
              sum2 += p[i] * p[i] / v;
            }
          }
        }
        DIALS_ASSERT(sum2 > 0);
        I[0] = sum1 / sum2;
        I[1] = sumv;
        if ((error_ = std::abs(I[0] - I0)) < eps) {
          break;
        }
        if (I[0] < minI+1e-3) {
          I[0] = minI+1e-3;
        }
        I0 = I[0];
      }
      DIALS_ASSERT(I[1] >= 0);

      if (niter_ >= max_iter-1) {
        niter_ = 0;
        error_ = 0;
        intensity_ = (sumc - sumb) / sums;
        variance_ = std::abs(intensity_ * sums) + sumb;
        correlation_ = 0;
        return;
      }

      // Set the intensity and variance
      intensity_ = I[0];
      variance_ = I[1];
      correlation_ = compute_correlation(p, m, c, b);
      //rmsd_ = compute_rmsd(I[0], p, m, c, b);
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

    /**
     * @returns The rmsd in the fit
     */
    /* double rmsd() const { */
    /*   return rmsd_; */
    /* } */

  private:

    /**
     * Evaluate the next intensity iteration.
     * @ returns The estimate of the intensity
     */
    vec2<double>
    estimate_intensity(const af::const_ref<FloatType, af::c_grid<3> > &p,
                       const af::const_ref<bool, af::c_grid<3> > &m,
                       const af::const_ref<FloatType, af::c_grid<3> > &c,
                       const af::const_ref<FloatType, af::c_grid<3> > &b,
                       double I) const {
      double sum1 = 0.0;
      double sum2 = 0.0;
      double sumv = 0.0;
      for (std::size_t i = 0; i < p.size(); ++i) {
        if (m[i]) {
          double v = std::abs(b[i]) + std::abs(I * p[i]);
          sumv += v;
          if (v > 0) {
            sum1 += (c[i] - b[i]) * p[i] / v;
            sum2 += p[i] * p[i] / v;
          }
        }
      }
      return vec2<double>(sum2 != 0 ? sum1 / sum2 : 0.0, sumv);
    }

    /**
     * Compute the correlation coefficient between the profile and reference
     */
    double
    compute_correlation(const af::const_ref<FloatType, af::c_grid<3> > &p,
                        const af::const_ref<bool, af::c_grid<3> > &m,
                        const af::const_ref<FloatType, af::c_grid<3> > &c,
                        const af::const_ref<FloatType, af::c_grid<3> > &b) const {
      double xb = 0.0, yb = 0.0;
      std::size_t count = 0;
      for (std::size_t i = 0; i < p.size(); ++i) {
        if (m[i]) {
          xb += intensity_*p[i] + b[i];
          yb += c[i];
          count++;
        }
      }
      DIALS_ASSERT(count > 0);
      xb /= count;
      yb /= count;
      double sdxdy = 0.0, sdx2 = 0.0, sdy2 = 0.0;
      for (std::size_t i = 0; i < p.size(); ++i) {
        if (m[i]) {
          double dx = (intensity_*p[i] + b[i]) - xb;
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

    /**
     * Compute the rmsd between the profile and reference
     */
    /* double */
    /* compute_rmsd(double I, */
    /*              const af::const_ref<FloatType, af::c_grid<3> > &p, */
    /*              const af::const_ref<bool, af::c_grid<3> > &m, */
    /*              const af::const_ref<FloatType, af::c_grid<3> > &c, */
    /*              const af::const_ref<FloatType, af::c_grid<3> > &b) const { */
    /*   double sum = 0.0; */
    /*   double ymax = 0.0; */
    /*   double ymin = 0.0; */
    /*   std::size_t count = 0; */
    /*   for (std::size_t i = 0; i < p.size(); ++i) { */
    /*     if (m[i]) { */
    /*       double y = (c[i] - b[i]); */
    /*       if (count == 0) { */
    /*         ymax = y; */
    /*         ymin = y; */
    /*       } else { */
    /*         if (ymax < y) ymax = y; */
    /*         if (ymin > y) ymin = y; */
    /*       } */
    /*       double x = I * p[i] - y; */
    /*       sum += x * x; */
    /*       count++; */
    /*     } */
    /*   } */
    /*   DIALS_ASSERT(count > 0); */
    /*   DIALS_ASSERT(ymax > ymin); */
    /*   sum /= count; */
    /*   return std::sqrt(sum) / (ymax - ymin); */
    /* } */

    double intensity_;
    double variance_;
    double correlation_;
    std::size_t niter_;
    double error_;
    /* double rmsd_; */
  };

}}

#endif /* DIALS_ALGORITHMS_INTEGRATION_FIT_FITTING_H */
