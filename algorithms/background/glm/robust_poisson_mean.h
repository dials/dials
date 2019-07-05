/*
 * robust_poisson_mean.h
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef SCITBX_GLMTBX_ROBUST_POISSON_MEAN_H
#define SCITBX_GLMTBX_ROBUST_POISSON_MEAN_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <scitbx/matrix/inversion.h>
#include <scitbx/matrix/multiply.h>
#include <scitbx/glmtbx/family.h>
#include <scitbx/glmtbx/robust_glm.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  /**
   * An algorithm to do robust generalized linear model as described in
   * Cantoni and Rochetti (2001) "Robust Inference for Generalized Linear
   * Models"
   */
  class RobustPoissonMean {
    typedef scitbx::glmtbx::poisson family;

  public:
    /**
     * Compute the generalized linear model using iteratively reweighted least
     * squares. The input expects a design matrix of size (nobs, ncoef), a list
     * of observations of size (nobs) and a list of initial estimates of size
     * (ncoef).
     * @param X The design matrix
     * @param Y The observations
     * @param B The initial estimate
     * @param c The huber tuning constant
     * @param tolerance The stopping critera
     * @param max_iter The maximum number of iterations
     */
    RobustPoissonMean(const af::const_ref<double> &Y,
                      double mean0,
                      double c,
                      double tolerance,
                      std::size_t max_iter)
        : niter_(0), error_(0), c_(c), tolerance_(tolerance), max_iter_(max_iter) {
      SCITBX_ASSERT(Y.size() > 0);
      SCITBX_ASSERT(mean0 > 0);
      SCITBX_ASSERT(c > 0);
      SCITBX_ASSERT(tolerance > 0);
      SCITBX_ASSERT(max_iter > 0);
      beta_ = std::log(mean0);
      compute(Y);
    }

    /**
     * @returns The poisson mean
     */
    double mean() const {
      DIALS_ASSERT(beta_ > -300 && beta_ < 300);
      return std::exp(beta_);
    }

    /**
     * @returns The number of iterations
     */
    std::size_t niter() const {
      return niter_;
    }

    /**
     * @returns The reletive error at the last iteration
     */
    double error() const {
      return error_;
    }

    /**
     * @returns Did the algorithm converge
     */
    bool converged() const {
      return niter_ < max_iter_;
    }

  private:
    void compute(const af::const_ref<double> &Y) {
      // Number of observations and coefficients
      std::size_t n_obs = Y.size();

      // Loop until we reach the maximum number of iterations
      for (niter_ = 0; niter_ < max_iter_; ++niter_) {
        // Initialize the sum to zero
        double U = 0;
        double H = 0;
        double w = 1.0;
        double eta = beta_;
        double mu = family::linkinv(eta);
        double var = family::variance(mu);
        double dmu = family::dmu_deta(eta);
        double phi = family::dispersion();
        SCITBX_ASSERT(phi > 0);
        SCITBX_ASSERT(var > 0);
        double svar = std::sqrt(phi * var);

        // Compute expectation values
        scitbx::glmtbx::expectation<family> epsi(mu, svar, c_);

        // The value of the b diagonal parts
        double b = epsi.epsi2 * w * dmu * dmu / svar;

        // Build the matrices from the observations
        for (std::size_t i = 0; i < n_obs; ++i) {
          // Compute some required values
          double res = (Y[i] - mu) / svar;

          // Compute the difference between psi and its expected value
          double psi = scitbx::glmtbx::huber(res, c_);
          double psi_m_epsi = psi - epsi.epsi1;

          // Compute the value of Psi and B_diag for this observation
          double q = psi_m_epsi * w * dmu / svar;

          // Update the BX = B * X and U matrices
          U += q;
          H += b;
        }

        // Compute delta = H^-1 U
        U = U / H;

        // Compute the relative error in the parameters and update
        double sum_delta_sq = U * U;
        double sum_beta_sq = beta_ * beta_;
        beta_ += U;

        // If error is within tolerance then break
        error_ = std::sqrt(sum_delta_sq / std::max(1e-10, sum_beta_sq));
        if (error_ < tolerance_) {
          break;
        }
      }
    }

    double beta_;
    std::size_t niter_;
    double error_;
    double c_;
    double tolerance_;
    std::size_t max_iter_;
  };

}}  // namespace dials::algorithms

#endif  // SCITBX_GLMTBX_ROBUST_POISSON_MEAN_H
