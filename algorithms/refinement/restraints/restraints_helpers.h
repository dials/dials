/*
 * restraints_helpers.h
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman.
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_REFINEMENT_RESTRAINTS_HELPERS_H
#define DIALS_REFINEMENT_RESTRAINTS_HELPERS_H

#ifndef RAD2DEG
#define RAD2DEG(x) ((x)*57.29577951308232087721)
#endif

#include <scitbx/mat3.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/math/angle_derivative.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace refinement {

  using scitbx::mat3;
  using scitbx::vec3;
  using scitbx::math::angle_derivative_wrt_vectors;

  /**
   * Given a reciprocal space orthogonalisation matrix B and the derivatives of
   * that matrix with respect to some parameters, calculate the derivatives of
   * the real space cell with respect to those parameters.
   *
   * That is, from [B] and d[B]/dp, calculate da/dp, db/dp, dc/dp, daa/dp,
   * dbb/dp and dcc/dp.
   */
  class CalculateCellGradients {
  public:
    CalculateCellGradients(const mat3<double> &B,
                           const af::const_ref<mat3<double> > &dB_dp) {
      // calculate the real space orthogonalisation matrix and its derivatives
      Omat_ = (B.transpose()).inverse();
      for (std::size_t i = 0; i < dB_dp.size(); ++i) {
        mat3<double> dBT = dB_dp[i].transpose();
        dO_dp_.push_back(-Omat_ * dBT * Omat_);
      }

      // extract the real space cell vectors and lengths
      avec_ = Omat_.get_column(0);
      a_ = avec_.length();
      bvec_ = Omat_.get_column(1);
      b_ = bvec_.length();
      cvec_ = Omat_.get_column(2);
      c_ = cvec_.length();

      // calculate the derivatives of real space angles with respect to the
      // vectors that form them
      af::tiny<vec3<double>, 2> dalpha = angle_derivative_wrt_vectors(bvec_, cvec_);
      af::tiny<vec3<double>, 2> dbeta = angle_derivative_wrt_vectors(avec_, cvec_);
      af::tiny<vec3<double>, 2> dgamma = angle_derivative_wrt_vectors(avec_, bvec_);
      dalpha_db_ = dalpha[0];
      dalpha_dc_ = dalpha[1];
      dbeta_da_ = dbeta[0];
      dbeta_dc_ = dbeta[1];
      dgamma_da_ = dgamma[0];
      dgamma_db_ = dgamma[1];
    }

    // gradients of parameter a
    af::shared<double> da_dp() {
      af::shared<double> result;
      for (std::size_t i = 0; i < dO_dp_.size(); ++i) {
        vec3<double> dav_dp = dO_dp_[i].get_column(0);
        // 1./a * avec.dot(dav_dp)
        result.push_back(1. / a_ * avec_ * dav_dp);
      }
      return result;
    }

    // gradients of parameter b
    af::shared<double> db_dp() {
      af::shared<double> result;
      for (std::size_t i = 0; i < dO_dp_.size(); ++i) {
        vec3<double> dbv_dp = dO_dp_[i].get_column(1);
        // 1./b * bvec.dot(dbv_dp)
        result.push_back(1. / b_ * bvec_ * dbv_dp);
      }
      return result;
    }

    // gradients of parameter c
    af::shared<double> dc_dp() {
      af::shared<double> result;
      for (std::size_t i = 0; i < dO_dp_.size(); ++i) {
        vec3<double> dcv_dp = dO_dp_[i].get_column(2);
        // 1./c * cvec.dot(dcv_dp)
        result.push_back(1. / c_ * cvec_ * dcv_dp);
      }
      return result;
    }

    // gradients of parameter alpha
    af::shared<double> daa_dp() {
      af::shared<double> result;
      for (std::size_t i = 0; i < dO_dp_.size(); ++i) {
        vec3<double> dbv_dp = dO_dp_[i].get_column(1);
        vec3<double> dcv_dp = dO_dp_[i].get_column(2);
        // dbv_dp.dot(dalpha_db) + dcv_dp.dot(dalpha_dc)
        result.push_back(RAD2DEG(dbv_dp * dalpha_db_ + dcv_dp * dalpha_dc_));
      }
      return result;
    }

    // gradients of parameter beta
    af::shared<double> dbb_dp() {
      af::shared<double> result;
      for (std::size_t i = 0; i < dO_dp_.size(); ++i) {
        vec3<double> dav_dp = dO_dp_[i].get_column(0);
        vec3<double> dcv_dp = dO_dp_[i].get_column(2);
        // dav_dp.dot(dbeta_da) + dcv_dp.dot(dbeta_dc)
        result.push_back(RAD2DEG(dav_dp * dbeta_da_ + dcv_dp * dbeta_dc_));
      }
      return result;
    }

    // gradients of parameter gamma
    af::shared<double> dcc_dp() {
      af::shared<double> result;
      for (std::size_t i = 0; i < dO_dp_.size(); ++i) {
        vec3<double> dav_dp = dO_dp_[i].get_column(0);
        vec3<double> dbv_dp = dO_dp_[i].get_column(1);
        // dav_dp.dot(dgamma_da) + dbv_dp.dot(dgamma_db)
        result.push_back(RAD2DEG(dav_dp * dgamma_da_ + dbv_dp * dgamma_db_));
      }
      return result;
    }

  private:
    mat3<double> Omat_;
    double a_, b_, c_;
    af::shared<mat3<double> > dO_dp_;
    vec3<double> avec_, bvec_, cvec_;
    vec3<double> dalpha_db_, dalpha_dc_;
    vec3<double> dbeta_da_, dbeta_dc_;
    vec3<double> dgamma_da_, dgamma_db_;
  };
}}  // namespace dials::refinement

#endif  // DIALS_REFINEMENT_RESTRAINTS_HELPERS_H
