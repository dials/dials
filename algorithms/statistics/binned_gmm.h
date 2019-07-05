/*
 * binned_gmm.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_STATISTICS_BINNED_GMM_H
#define DIALS_ALGORITHMS_STATISTICS_BINNED_GMM_H

#include <cmath>
#include <boost/math/special_functions/erf.hpp>
#include <scitbx/constants.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using boost::math::erf;
  using scitbx::constants::pi;

  namespace detail {

    /**
     * Compute the expectations required for the 1d gaussian model
     */
    class Expectation {
    public:
      /**
       * Compute the expectations
       * @param a The lower bound of the bin
       * @param b The upper bound of the bin
       * @param mu The mean parameters
       * @param sigma The standard deviation parameter
       */
      Expectation(double a, double b, double mu, double sigma) {
        double e1 = erf((b - mu) / (std::sqrt(2.0) * sigma));
        double e2 = erf((a - mu) / (std::sqrt(2.0) * sigma));
        double e3 = std::exp(-(a - mu) * (a - mu) / (2 * sigma * sigma))
                    / (std::sqrt(2.0 * pi) * sigma);
        double e4 = std::exp(-(b - mu) * (b - mu) / (2 * sigma * sigma))
                    / (std::sqrt(2.0 * pi) * sigma);
        expectation0_ = 0.5 * (e1 - e2);
        expectation1_ =
          0.5 * mu * (e1 - e2) + (sigma * sigma / std::sqrt(2 * pi)) * (e3 - e4);
        expectation2_ = (sigma * sigma / 2.0) * (e1 - e2)
                        + sigma * sigma * ((a - mu) * e3 - (b - mu) * e4);
      }

      /**
       * Return the expectation
       */
      double expectation0() {
        return expectation0_;
      }

      /**
       * Return the expectation
       */
      double expectation1() {
        return expectation1_;
      }

      /**
       * Return the expectation
       */
      double expectation2() {
        return expectation2_;
      }

    protected:
      double expectation0_;
      double expectation1_;
      double expectation2_;
    };

  }  // namespace detail

  /**
   * Compute the binned gaussian model withe a fixed mean using expectation
   * maximization
   */
  class BinnedGMMSingle1DFixedMean {
  public:
    /**
     * Compute the parameters
     * @param a The lower bounds of the bin
     * @param b The upper bounds of the bin
     * @param n The number of counts in the bin
     * @param mu The mean parameter estimate
     * @param sigma The sigma parameter estimate
     * @param epsilon The convergence tolerance
     * @param max_iter The maximum number of iterations
     */
    BinnedGMMSingle1DFixedMean(const af::const_ref<double> &a,
                               const af::const_ref<double> &b,
                               const af::const_ref<double> &n,
                               double mu,
                               double sigma,
                               double epsilon,
                               std::size_t max_iter)
        : max_iter_(max_iter), num_iter_(0), epsilon_(epsilon), mu_(mu), sigma_(sigma) {
      // Check the input
      DIALS_ASSERT(epsilon > 0);
      DIALS_ASSERT(max_iter > 1);
      DIALS_ASSERT(sigma > 0);
      DIALS_ASSERT(a.size() == b.size());
      DIALS_ASSERT(a.size() == n.size());

      // The E step in this case is constant
      af::const_ref<double> m = n;
      double c = af::sum(m);
      double logL0 = 0.0;

      af::shared<bool> use(a.size(), true);
      for (std::size_t i = 0; i < m.size(); ++i) {
        detail::Expectation e(a[i], b[i], mu_, sigma_);
        if (m[i] < 1) {
          use[i] = false;
        }
      }

      // Do the iterations
      for (num_iter_ = 0; num_iter_ < max_iter; ++num_iter_) {
        double sum_n_log_p = 0;
        double sum_n = 0;
        double sum_log_p = 0;
        double va_new = 0;
        double min_P = 1;

        double sum_P = 0;
        double sum_E2 = 0;

        for (std::size_t i = 0; i < m.size(); ++i) {
          detail::Expectation e(a[i], b[i], mu_, sigma_);
          double P = e.expectation0();
          double E2 = e.expectation2();
          if (use[i]) {
            if (P > 1e-100) {
              va_new += m[i] * E2 / P;
            } else {
            }
            sum_n_log_p += n[i] * (P > 1e-100 ? std::log(P) : std::log(1e-100));
            sum_n += n[i];
            sum_log_p += (P > 1e-100 ? std::log(P) : std::log(1e-100));
            if (P < min_P) min_P = P;
            sum_P += P;
            sum_E2 += E2;
          }
        }

        /* double Pr = (1.0 - sum_P); */
        /* double mr = sum_n * Pr; */
        /* double E2r = (sigma_*sigma_ - sum_E2); */
        /* double XX = mr * E2r / Pr; */
        /* va_new += XX; */

        sigma_ = std::sqrt(va_new / c);

        // Check if we've converged
        double logL = sum_n_log_p - sum_n * sum_log_p;
        double error = std::abs((logL - logL0) / std::max(logL0, 1e-10));
        if (num_iter_ > 0 && error < epsilon) {
          break;
        }
        logL0 = logL;
      }
    }

    /**
     * @returns The maximum number of iterations
     */
    std::size_t max_iter() const {
      return max_iter_;
    }

    /**
     * @returns The number of iterations
     */
    std::size_t num_iter() const {
      return num_iter_;
    }

    /**
     * @returns The epsilon
     */
    double epsilon() const {
      return epsilon_;
    }

    /**
     * @returns The mu parameter estimate
     */
    double mu() const {
      return mu_;
    }

    /**
     * @returns The sigma parameter estimate
     */
    double sigma() const {
      return sigma_;
    }

  protected:
    std::size_t max_iter_;
    std::size_t num_iter_;
    double epsilon_;
    double mu_;
    double sigma_;
  };

  /**
   * Compute the binned gaussian model using expectation
   * maximization
   */
  class BinnedGMMSingle1D {
  public:
    /**
     * Compute the parameters
     * @param a The lower bounds of the bin
     * @param b The upper bounds of the bin
     * @param n The number of counts in the bin
     * @param mu The mean parameter estimate
     * @param sigma The sigma parameter estimate
     * @param epsilon The convergence tolerance
     * @param max_iter The maximum number of iterations
     */
    BinnedGMMSingle1D(const af::const_ref<double> &a,
                      const af::const_ref<double> &b,
                      const af::const_ref<double> &n,
                      double mu,
                      double sigma,
                      double epsilon,
                      std::size_t max_iter)
        : max_iter_(max_iter), num_iter_(0), epsilon_(epsilon), mu_(mu), sigma_(sigma) {
      // Check the input
      DIALS_ASSERT(epsilon > 0);
      DIALS_ASSERT(max_iter > 1);
      DIALS_ASSERT(sigma > 0);
      DIALS_ASSERT(a.size() == b.size());
      DIALS_ASSERT(a.size() == n.size());

      // The E step in this case is constant
      af::const_ref<double> m = n;
      double c = af::sum(m);
      double logL0 = 0.0;

      // Do the iterations
      for (num_iter_ = 0; num_iter_ < max_iter; ++num_iter_) {
        double sum_n_log_p = 0;
        double sum_n = 0;
        double sum_log_p = 0;
        double mu_new = 0;
        double va_new = 0;
        for (std::size_t i = 0; i < m.size(); ++i) {
          detail::Expectation e(a[i], b[i], mu_, sigma_);
          double P = e.expectation0();
          if (n[i] > 0 && P > 1e-10) {
            mu_new += m[i] * e.expectation1() / P;
            va_new += m[i] * e.expectation2() / P;
            sum_n_log_p += n[i] * std::log(P);
            sum_n += n[i];
            sum_log_p += std::log(P);
          }
        }
        mu_ = mu_new / c;
        sigma_ = std::sqrt(va_new / c);

        // Check the convergence
        double logL = sum_n_log_p - sum_n * sum_log_p;
        double error = std::abs((logL - logL0) / std::max(logL0, 1e-10));
        if (num_iter_ > 0 && error < epsilon) {
          break;
        }
        logL0 = logL;
      }
    }

    /**
     * @returns The maximum number of iterations
     */
    std::size_t max_iter() const {
      return max_iter_;
    }

    /**
     * @returns The number of iterations
     */
    std::size_t num_iter() const {
      return num_iter_;
    }

    /**
     * @returns The epsilon
     */
    double epsilon() const {
      return epsilon_;
    }

    /**
     * @returns The mu parameter estimate
     */
    double mu() const {
      return mu_;
    }

    /**
     * @returns The sigma parameter estimate
     */
    double sigma() const {
      return sigma_;
    }

  protected:
    std::size_t max_iter_;
    std::size_t num_iter_;
    double epsilon_;
    double mu_;
    double sigma_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_STATISTICS_BINNED_GMM_H
