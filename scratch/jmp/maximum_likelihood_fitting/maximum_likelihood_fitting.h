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
#include <boost/math/special_functions/gamma.hpp>
#include <scitbx/vec2.h>
#include <scitbx/mat2.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using boost::math::lgamma;
  using scitbx::vec2;
  using scitbx::mat2;

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
      DIALS_ASSERT(sum_s1_ > 0);
      DIALS_ASSERT(sum_s2_ > 0);
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
        sum_s1c += c_[i]*s1_[i]*inv_den;
        sum_s2c += c_[i]*s2_[i]*inv_den;
      }
      double dlds1 = sum_s1c - sum_s1_;
      double dlds2 = sum_s2c - sum_s2_;
      double dS1 = step_*dlds1;
      double dS2 = step_*dlds2;
      double S1_new = S1_ + dS1;
      double S2_new = S2_ + dS2;

      // Hack to make sure we obey the constraints
      if (K11_ * S1_new <= -K12_ * S2_new) {
        S1_new = S1_ / 2.0;
      }
      if (K22_ * S2_new <= -K21_ * S1_new) {
        S2_new = S2_ / 2.0;
      }
      S1_ = S1_new;
      S2_ = S2_new;
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

  class MLPoisson2Objective {
  public:

    MLPoisson2Objective(
          const af::const_ref<double> &c,
          const af::const_ref<double> &s1,
          const af::const_ref<double> &s2)
        : c_(c.begin(), c.end()),
          s1_(s1.begin(), s1.end()),
          s2_(s2.begin(), s2.end()) {
      DIALS_ASSERT(c.size() == s1.size());
      DIALS_ASSERT(c.size() == s2.size());
      sums1_ = af::sum(s1);
      sums2_ = af::sum(s2);
      sumcfac_ = 0;
      for (std::size_t i = 0; i < c_.size(); ++i) {
        sumcfac_ += lgamma(c_[i] + 1);
      }
      C1_ = 0.0;
      C2_ = 0.0;
      bool S1T = false;
      bool S2T = false;
      for (std::size_t i = 0; i < s1.size(); ++i) {
        if (s1[i] != 0) {
          double m = s2[i] / s1[i];
          if (!S1T || m < C1_) {
            C1_ = m;
          }
          S1T = true;
        }
        if (s2[i] != 0) {
          double m = s1[i] / s2[i];
          if (!S2T || m < C2_) {
            C2_ = m;
          }
          S2T = true;
        }
      }
      DIALS_ASSERT(S1T && S2T);
    }

    double C1() const {
      return C1_;
    }

    double C2() const {
      return C2_;
    }

    double f(vec2<double> x) const {
      double sumclog = 0;
      for (std::size_t i = 0; i < c_.size(); ++i) {
        double den = x[0]*s1_[i] + x[1]*s2_[i];
        DIALS_ASSERT(den > 0);
        sumclog += c_[i] * std::log(den);
      }
      return sumcfac_ + x[0]*sums1_ + x[1]*sums2_ - sumclog;
    }

    vec2<double> df(vec2<double> x) const {
      double sumcs1 = 0;
      double sumcs2 = 0;
      for (std::size_t i = 0; i < c_.size(); ++i) {
        double den = x[0]*s1_[i] + x[1]*s2_[i];
        DIALS_ASSERT(den > 0);
        double denr = 1.0 / den;
        sumcs1 += c_[i]*s1_[i] * denr;
        sumcs2 += c_[i]*s2_[i] * denr;
      }
      return vec2<double>(sums1_ - sumcs1, sums2_ - sumcs2);
    }

    double lsf(vec2<double> x0, vec2<double> p, double alpha) const {
      return f(x0 + alpha * p);
    }

    double lsdf(vec2<double> x0, vec2<double> p, double alpha) const {
      double sumps1 = p[0]*sums1_;
      double sumps2 = p[1]*sums2_;
      double sumc = 0;
      for (std::size_t i = 0; i < c_.size(); ++i) {
        double pbs = p[0]*s1_[i] + p[1]*s2_[i];
        double bs0 = x0[0]*s1_[i] + x0[1]*s2_[i];
        double den = bs0 + alpha * pbs;
        DIALS_ASSERT(den > 0);
        double denr = 1.0 / den;
        sumc += c_[i]*pbs * denr;
      }
      return sumps1 + sumps2 - sumc;
    }

  private:

    std::vector<double> c_;
    std::vector<double> s1_;
    std::vector<double> s2_;
    double sums1_;
    double sums2_;
    double sumcfac_;
    double C1_;
    double C2_;
  };

  class MLPoisson2Stepper {
  public:

    MLPoisson2Stepper(
          const af::const_ref<double> &c,
          const af::const_ref<double> &s1,
          const af::const_ref<double> &s2,
          vec2<double> x)
      : obj_(c, s1, s2),
        X0_(x),
        B_(1, 0, 0, 1) {

    }

    void step() {
      const mat2<double> I(1.0, 0.0, 0.0, 1.0);
      vec2<double> DF0 = obj_.df(X0_);
      vec2<double> P = -B_*DF0;
      double cden1 = (P[0]+P[1]*obj_.C1());
      double cden2 = (obj_.C2()*P[0]+P[1]);
      double alpha_max1 = -(X0_[0] + X0_[1]*obj_.C1()) / cden1;
      double alpha_max2 = -(obj_.C2()*X0_[0] + X0_[1]) / cden2;
      double alpha_max = 100;
      if (alpha_max1 > 0) {
        alpha_max = std::min(alpha_max, alpha_max1);
      }
      if (alpha_max2 > 0) {
        alpha_max = std::min(alpha_max, alpha_max2);
      }
      std::cout << alpha_max << std::endl;
      double alpha = linesearch(P, alpha_max);
      vec2<double> S = alpha * P;
      vec2<double> X1 = X0_ + S;
      vec2<double> DF1 = obj_.df(X1);
      vec2<double> Y = DF1 - DF0;
      vec2<double> YT = Y;
      vec2<double> ST = S;
      double YTS = YT * S;
      mat2<double> SYT(
          S[0] * YT[0], S[0] * YT[1],
          S[1] * YT[0], S[1] * YT[1]);
      mat2<double> YST(
          Y[0] * ST[0], Y[0] * ST[1],
          Y[1] * ST[0], Y[1] * ST[1]);
      mat2<double> SST(
          S[0] * ST[0], S[0] * ST[1],
          S[1] * ST[0], S[1] * ST[1]);
      B_ = (I - SYT/YTS)*B_*(I - YST/YTS) + SST/YTS;
      /* std::cout << SYT[0] << ", " << SYT[1] << ", " << SYT[2] << ", " << SYT[3] << std::endl; */
      /* std::cout << YST[0] << ", " << YST[1] << ", " << YST[2] << ", " << YST[3] << std::endl; */
      /* std::cout << SST[0] << ", " << SST[1] << ", " << SST[2] << ", " << SST[3] << std::endl; */
      /* std::cout << B_[0] << ", " << B_[1] << ", " << B_[2] << ", " << B_[3] << std::endl; */
      X0_ = X1;
    }

    vec2<double> X() const {
      return X0_;
    }

    mat2<double> B() const {
      return B_;
    }

  private:

    double linesearch(vec2<double> P, double alpha_max) const {
      double c1 = 1e-4;
      double c2 = 0.9;
      std::size_t maxiter = 10;
      double alpha0 = 0;
      double alpha1 = alpha_max / pow(1.618, 10);
      double f00 = lsf(P, alpha0);
      double d00 = lsdf(P, alpha0);
      double f0 = f00;
      double d0 = d00;
      DIALS_ASSERT(d00 < 0);
      for (std::size_t i = 0; i < maxiter; ++i) {
        double f1 = lsf(P, alpha1);
        if (f1 > f00 + c1*alpha1*d00 || (f1 >= f0 && i > 0)) {
          alpha1 = zoom(P, alpha0, alpha1, f00, d00, f0);
          break;
        }
        double d1 = lsdf(P, alpha1);
        if (std::abs(d1) <= -c2*d00) {
          break;
        }
        if (d1 >= 0) {
          alpha1 = zoom(P, alpha1, alpha0, f00, d00, f1);
          break;
        }
        f0 = f1;
        d0 = d1;
        alpha0 = alpha1;
        alpha1 *= 1.618;
      }
      return alpha1;
    }

    double zoom(vec2<double> P, double alpha0, double alpha1,
        double f00, double d00, double f0) const {
      double c1 = 1e-4;
      double c2 = 0.9;
      for (std::size_t i = 0; ; ++i) {
        double alpha = 0.5 * (alpha0 + alpha1);
        double f1 = lsf(P, alpha);
        double d1 = lsdf(P, alpha);
        if (f1 > f00 + c1*alpha*d00 || (f1 >= f0 && i > 0)) {
          alpha1 = alpha;
        } else {
          if (std::abs(d1) <= -c2*d00) {
            return alpha;
          }
          if (d1 * (alpha1 - alpha0) >= 0) {
            alpha1 = alpha0;
          }
          alpha0 = alpha;
        }
      }
    }

    double lsf(vec2<double> P, double alpha) const {
      return obj_.lsf(X0_, P, alpha);
    }

    double lsdf(vec2<double> P, double alpha) const {
      return obj_.lsdf(X0_, P, alpha);
    }

    MLPoisson2Objective obj_;
    vec2<double> X0_;
    mat2<double> B_;
  };

}} // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_PROFILE_MAXIMUM_LIKELIHOOD_FITTING_H */
