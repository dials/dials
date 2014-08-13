/*
 * maximum_likelihood_fitting.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_PROFILE_MAXIMUM_LIKELIHOOD_FITTING_H
#define DIALS_ALGORITHMS_INTEGRATION_PROFILE_MAXIMUM_LIKELIHOOD_FITTING_H

#include <vector>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  /**
   * A class to do a maximum likelihood estimate of 2 poisson distributed
   * signals.
   */
  class MLPoisson2 {
  public:

    MLPoisson2(
          const af::const_ref<double> &c,
          const af::const_ref<double> &s1,
          const af::const_ref<double> &s2,
          double S1,
          double S2,
          double step)
        : c_(c.begin(), c.end()),
          s1_(s1.begin(), s1.end()),
          s2_(s2.begin(), s2.end()) {
      DIALS_ASSERT(step > 0);
      DIALS_ASSERT(c.size() == s1.size());
      DIALS_ASSERT(c.size() == s2.size());
      step_ = step;
      S1_ = S1;
      S2_ = S2;
      sum_s1_ = af::sum(s1);
      sum_s2_ = af::sum(s2);
      DIALS_ASSERT(c.all_ge(0));
      DIALS_ASSERT(s1.all_ge(0));
      DIALS_ASSERT(s2.all_ge(0));
      K11_ = 1.0;
      K12_ = 0.0;
      K21_ = 0.0;
      K22_ = 1.0;
      bool S1T = false;
      bool S2T = false;
      for (std::size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != 0) {
          double m = s2[i] / s1[i];
          if (!S1T || m < K12_) {
            K12_ = m;
          }
          S1T = true;
        }
        if (s2[i] != 0) {
          double m = s1[i] / s2[i];
          if (!S2T || m < K21_) {
            K21_ = m;
          }
          S2T = true;
        }
      }
      DIALS_ASSERT(S1T && S2T);
    }

    void step() {
      
      // Evaluate the gradient and do a step
      double sum_s1c = 0;
      double sum_s2c = 0;
      for (std::size_t i = 0; i < c_.size(); ++i) {
        double den = S1_*s1_[i]+S2_*s2_[i];
        DIALS_ASSERT(den > 0);
        double inv_den = 1.0 / den;
        sum_s1c = c_[i]*s1_[i]*inv_den;
        sum_s2c = c_[i]*s2_[i]*inv_den;
      }
      double dlds1 = sum_s1c - sum_s1_;
      double dlds2 = sum_s2c - sum_s2_;
      double dS1 = step_*dlds1;
      double dS2 = step_*dlds2;
      double S1_new = S1_ + dS1;
      double S2_new = S2_ + dS2;
      
      // Hack to make sure we obey the constraints
      if (S1_new*K11_ + S2_new*K12_ <= 0) {
        double num = S1_*S2_new-S2_*S1_new;
        double den = K12_*(S1_-S1_new)-K11_*(S2_-S2_new);
        DIALS_ASSERT(den != 0);
        double frac = num / den;
        double S1X = K11_*frac;
        double S2X = K12_*frac;
        S1_new = (S1X + S1_) / 2.0;
        S2_new = (S2X + S2_) / 2.0;
      }
      if (S1_new*K21_ + S2_new*K22_ <= 0) {
        double num = S1_*S2_new-S2_*S1_new;
        double den = K22_*(S1_-S1_new)-K21_*(S2_-S2_new);
        DIALS_ASSERT(den != 0);
        double frac = num / den;
        double S1X = K21_*frac;
        double S2X = K22_*frac;
        S1_new = (S1X + S1_) / 2.0;
        S2_new = (S2X + S2_) / 2.0;
      }
      S1_ = S1_new;
      S2_ = S1_new;
    }
    
    void solve(double eps) {
      for (;;) {
        double S1 = S1_;
        double S2 = S2_;
        step();
        if ((S1-S1_)*(S1-S1_)+(S2-S2_)*(S2-S2_) < eps) {
          break;
        }
      }
    }

    double S1() const {
      return S1_;
    }

    double S2() const {
      return S2_;
    }

  private:

    double step_;
    double S1_;
    double S2_;
    double sum_s1_;
    double sum_s2_;
    double K11_;
    double K12_;
    double K21_;
    double K22_;
    std::vector<double> c_;
    std::vector<double> s1_;
    std::vector<double> s2_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_MAXIMUM_LIKELIHOOD_FITTING_H */
