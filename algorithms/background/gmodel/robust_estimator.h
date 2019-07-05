/*
 * robust_estimator.h
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_GMODEL_ROBUST_ESTIMATOR_H
#define DIALS_ALGORITHMS_BACKGROUND_GMODEL_ROBUST_ESTIMATOR_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <scitbx/matrix/multiply.h>
#include <scitbx/glmtbx/family.h>
#include <scitbx/glmtbx/robust_glm.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  class ExpectationTable {
  public:
    ExpectationTable()
        : c_(-1), max_(1000), div_(100), size_(max_ * div_), epsi_table_(size_) {}

    void set_c(double c) {
      DIALS_ASSERT(c > 0);
      if (c_ <= 0 || std::abs(c - c_) > 1e-3) {
        c_ = c;
        for (std::size_t i = 0; i < size_; ++i) {
          epsi_table_[i] = calculate((double)i / (double)div_);
        }
      }
    }

    vec2<double> get(double mu) {
      return mu < max_ ? interpolate(mu) : calculate(mu);
    }

  private:
    vec2<double> interpolate(double mu) const {
      DIALS_ASSERT(mu >= 0);
      std::size_t index = (std::size_t)std::floor(mu * div_);
      DIALS_ASSERT(index < epsi_table_.size());
      return epsi_table_[index];
    }

    vec2<double> calculate(double mu) const {
      DIALS_ASSERT(mu > 0);
      DIALS_ASSERT(c_ > 0);
      scitbx::glmtbx::expectation<scitbx::glmtbx::poisson> e(mu, std::sqrt(mu), c_);
      return vec2<double>(e.epsi1, e.epsi2);
    }

    double c_;
    int max_;
    int div_;
    std::size_t size_;
    af::shared<vec2<double> > epsi_table_;
  };

  inline ExpectationTable get_expectation_table(double c) {
    static ExpectationTable table;
    table.set_c(c);
    return table;
  }

  /**
   * An algorithm to do robust generalized linear model as described in
   * Cantoni and Rochetti (2001) "Robust Inference for Generalized Linear
   * Models"
   */
  class robust_estimator {
  public:
    typedef scitbx::glmtbx::poisson family;

    /**
     * Compute the generalized linear model using iteratively reweighted least
     * squares. The input expects a design matrix of size (nobs, ncoef), a list
     * of observations of size (nobs) and a list of initial estimates of size
     * (ncoef).
     * @param X The background profile
     * @param Y The background observations
     * @param B The initial estimate
     * @param c The huber tuning constant
     * @param tolerance The stopping critera
     * @param max_iter The maximum number of iterations
     */
    robust_estimator(const af::const_ref<double> &X,
                     const af::const_ref<double> &Y,
                     double B,
                     double c,
                     double tolerance,
                     std::size_t max_iter)
        : beta_(B),
          niter_(0),
          error_(0),
          c_(c),
          tolerance_(tolerance),
          max_iter_(max_iter) {
      SCITBX_ASSERT(X.size() == Y.size());
      SCITBX_ASSERT(c > 0);
      SCITBX_ASSERT(tolerance > 0);
      SCITBX_ASSERT(max_iter > 0);
      compute(X, Y);
    }

    /**
     * @returns The parameters
     */
    double scale_parameter() const {
      return beta_;
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

    /**
     * Compute the values of mu, at X given the computed parameters
     * @param X The design matrix
     * @return The values
     */
    af::shared<double> mu(const af::const_ref<double> &X) const {
      af::shared<double> result(X.size());
      for (std::size_t i = 0; i < result.size(); ++i) {
        double eta = X[i] * beta_;
        result[i] = family::linkinv(eta);
      }
      return result;
    }

  private:
    void compute(const af::const_ref<double> &X, const af::const_ref<double> &Y) {
      // Number of observations and coefficients
      std::size_t n_obs = X.size();

      // Initialize the required matrices and vectors
      double U = 0;
      double H = 0;

      // Get the expectation table
      ExpectationTable expectation = get_expectation_table(c_);

      // Loop until we reach the maximum number of iterations
      for (niter_ = 0; niter_ < max_iter_; ++niter_) {
        // Initialize the sum to zero
        U = 0.0;

        // Build the matrices from the observations
        for (std::size_t i = 0; i < n_obs; ++i) {
          // Compute the values for eta
          double eta = X[i] * beta_;

          // Compute some required values
          double mu = family::linkinv(eta);
          double var = mu;
          double dmu = family::dmu_deta(eta);
          SCITBX_ASSERT(var > 0);
          double svar = std::sqrt(var);
          double res = (Y[i] - mu) / svar;

          // Compute expectation values
          vec2<double> epsi = expectation.get(mu);
          double epsi1 = epsi[0];
          double epsi2 = epsi[1];

          // Compute the difference between psi and its expected value
          double psi = scitbx::glmtbx::huber(res, c_);
          double psi_m_epsi = psi - epsi1;

          // Compute the value of Psi and B_diag for this observation
          double q = psi_m_epsi * dmu / svar;
          double b = epsi2 * dmu * dmu / svar;

          // Update the BX = B * X and U matrices
          U += q * X[i];
          H += X[i] * b * X[i];
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

#endif  // DIALS_ALGORITHMS_BACKGROUND_GMODEL_ROBUST_ESTIMATOR_H
