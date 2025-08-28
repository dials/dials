#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE1D_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_TOFPROFILE1D_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <array>
#include <cassert>

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/NonLinearOptimization>

/*
1D profile fitting method based on
Yano, N., Yamada, T., Hosoya, T. et al.
Application of profile fitting method to neutron time-of-flight protein
single crystal diffraction data collected at the iBIX. Sci Rep 6,
36628 (2016). https://doi.org/10.1038/srep36628
*/

inline double exp_safe(double x) {
  if (x > 700.0) return std::exp(700.0);
  if (x < -700.0) return std::exp(-700.0);
  return std::exp(x);
}
inline double erfc_safe(double x) {
  if (x > 10.0) return std::erfc(10.0);
  if (x < -10.0) return std::erfc(-10.0);
  return std::erfc(x);
}
inline bool is_finite_double(double x) {
  return std::isfinite(x);
}

inline double simpson_integrate(const std::vector<double>& y,
                                const std::vector<double>& x) {
  /*
   * Based on https://en.wikipedia.org/wiki/Simpson%27s_rule
   */

  const size_t n = y.size();
  if (n < 2) return 0.0;
  // If spacing is uniform we can use simple Simpson; handle general by composite
  // Simpson on subintervals
  if (n % 2 == 0) {
    // Points must be odd; for even n, drop last point
    // use trapezoid on final point
    double res = 0.0;
    size_t m = n - 1;
    for (size_t i = 0; i + 2 < m; i += 2) {
      double h0 = x[i + 1] - x[i];
      double h1 = x[i + 2] - x[i + 1];
      // map to local quadratic: use Simpson on three points i,i+1,i+2 with nonuniform
      // spacing via quadratic interpolation Simpler: res += (x[i+2] - x[i]) / 6.0 *
      // (y[i] + 4.0*y[i+1] + y[i+2]);
      res += (x[i + 2] - x[i]) / 6.0 * (y[i] + 4.0 * y[i + 1] + y[i + 2]);
    }
    // trapezoid for last interval (m-1,m)
    res += 0.5 * (x[n - 1] - x[n - 2]) * (y[n - 2] + y[n - 1]);
    return res;
  } else {
    // odd n -> standard composite Simpson across whole range
    double res = 0.0;
    for (size_t i = 0; i + 2 < n; i += 2) {
      res += (x[i + 2] - x[i]) / 6.0 * (y[i] + 4.0 * y[i + 1] + y[i + 2]);
    }
    return res;
  }
}

static std::vector<double> profile1d_func(const std::vector<double>& tof,
                                          double A,
                                          double alpha,
                                          double beta,
                                          double sigma,
                                          double T_ph) {
  const size_t m = tof.size();
  std::vector<double> out(m, 0.0);
  double sigma2 = sigma * sigma;
  double sigma_sqrt = std::sqrt(2.0 * sigma2);
  double N = (alpha * beta) / (2.0 * (alpha + beta + 1e-300));  // avoid /0

  for (size_t i = 0; i < m; ++i) {
    double dT = tof[i] - T_ph;
    double u = alpha * 0.5 * (alpha * sigma2 + 2.0 * dT);
    double v = beta * 0.5 * (beta * sigma2 - 2.0 * dT);
    double y = (alpha * sigma2 + dT) / (sigma_sqrt + 1e-300);
    double z = (beta * sigma2 - dT) / (sigma_sqrt + 1e-300);
    double exp_u = exp_safe(u);
    double exp_v = exp_safe(v);
    double erfc_y = erfc_safe(y);
    double erfc_z = erfc_safe(z);
    double val = A * N * (exp_u * erfc_y + exp_v * erfc_z);
    if (!is_finite_double(val)) val = 1e-12;
    out[i] = val;
  }
  return out;
}

// Eigen Functor wrapper
template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor {
  typedef _Scalar Scalar;
  enum { InputsAtCompileTime = NX, ValuesAtCompileTime = NY };
  typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

  int m_inputs, m_values;
  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}
  int inputs() const {
    return m_inputs;
  }
  int values() const {
    return m_values;
  }
};

struct TOFProfileFunctor : Functor<double> {
  const std::vector<double>& tof;
  const std::vector<double>& y_norm;  // assumed normalized
  std::array<double, 5> min_bounds;
  std::array<double, 5> max_bounds;

  TOFProfileFunctor(const std::vector<double>& tof_,
                    const std::vector<double>& y_norm_,
                    const std::array<double, 5>& minb,
                    const std::array<double, 5>& maxb)
      : Functor<double>(5, static_cast<int>(tof_.size())), tof(tof_), y_norm(y_norm_) {
    min_bounds = minb;
    max_bounds = maxb;
  }

  inline Eigen::VectorXd clamp_params(const Eigen::VectorXd& x) const {
    Eigen::VectorXd xc = x;
    for (int i = 0; i < x.size(); ++i) {
      xc[i] = std::min(std::max(x[i], min_bounds[i]), max_bounds[i]);
    }
    return xc;
  }

  int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
    Eigen::VectorXd xc = clamp_params(x);
    double A = xc[0];
    double alpha = xc[1];
    double beta = xc[2];
    double sigma = xc[3];
    double T_ph = xc[4];

    std::vector<double> model = profile1d_func(tof, A, alpha, beta, sigma, T_ph);
    assert(model.size() == static_cast<size_t>(values()));
    for (int i = 0; i < values(); ++i) {
      fvec[i] = y_norm[i] - model[i];
    }
    return 0;
  }

  // numerical Jacobian (central differences)
  int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) const {
    const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
    Eigen::VectorXd xc = clamp_params(x);
    int p = inputs();
    int m = values();

    Eigen::VectorXd f0(m);
    operator()(xc, f0);

    for (int j = 0; j < p; ++j) {
      Eigen::VectorXd xp = xc;
      Eigen::VectorXd xm = xc;
      double step = (std::abs(xc[j]) + 1.0) * eps;
      xp[j] += step;
      xm[j] -= step;
      Eigen::VectorXd fp(m), fm(m);
      operator()(xp, fp);
      operator()(xm, fm);
      for (int i = 0; i < m; ++i) {
        J(i, j) = (fp[i] - fm[i]) / (2.0 * step);
      }
    }
    return 0;
  }
};

class TOFProfile1D {
public:
  std::vector<double> tof;
  std::vector<double> intensities;  // raw intensities
  std::vector<double> y_norm;       // normalized intensities
  double intensity_max;

  // params
  double A, alpha, beta, sigma, T_ph;
  std::array<double, 5> min_bounds;
  std::array<double, 5> max_bounds;

  TOFProfile1D(const std::vector<double>& tof_,
               const std::vector<double>& intensities_,
               double A_,
               double alpha_,
               double beta_,
               double sigma_,
               double T_ph_,
               const std::array<double, 5>& minb = {1.0, 0.0, 0.0, 1.0, 0.0},
               const std::array<double, 5>& maxb = {1e6, 1.0, 1e5, 1e7, 1e6})
      : tof(tof_),
        intensities(intensities_),
        A(A_),
        alpha(alpha_),
        beta(beta_),
        sigma(sigma_),
        T_ph(T_ph_),
        min_bounds(minb),
        max_bounds(maxb) {
    assert(tof.size() == intensities.size());

    // Ensure no negative values
    for (auto& v : intensities)
      if (v < 0) v = 0.0;
    intensity_max = 1.0;
    if (!intensities.empty()) {
      intensity_max = *std::max_element(intensities.begin(), intensities.end());
      if (intensity_max <= 0.0) intensity_max = 1.0;
    }

    // Rescale max A based on max intensity
    max_bounds[0] = 4.0 * intensity_max;

    // build normalized y vector
    const size_t n = intensities.size();
    y_norm.resize(n);
    for (size_t i = 0; i < n; ++i) {
      double v = intensities[i];
      if (!is_finite_double(v)) v = 0.0;
      y_norm[i] = v / intensity_max;
    }
  }

  // compute model (unnormalized)
  std::vector<double> model(double A_,
                            double alpha_,
                            double beta_,
                            double sigma_,
                            double T_ph_) const {
    std::vector<double> m = profile1d_func(tof, A_, alpha_, beta_, sigma_, T_ph_);
    for (auto& v : m)
      v *= intensity_max;
    return m;
  }

  // return fitted (unnormalized) model using current params
  std::vector<double> result() const {
    return model(A, alpha, beta, sigma, T_ph);
  }

  double calc_intensity() const {
    std::vector<double> r = result();
    return simpson_integrate(r, tof);
  }

  double calc_variance(const std::vector<double>& projected_variance) const {
    std::vector<double> profile = result();  // current best-fit model (includes A)
    DIALS_ASSERT(profile.size() == projected_variance.size());

    // Compute normalized profile shape (divide out amplitude A)
    std::vector<double> p(profile.size());
    double A_local = A;
    if (A_local == 0.0) return 0.0;
    for (size_t i = 0; i < p.size(); ++i) {
      p[i] = profile[i] / A_local;
    }

    // Variance of amplitude A
    double denom = 0.0;
    for (size_t i = 0; i < p.size(); ++i) {
      double var = projected_variance[i];
      if (var <= 1e-7) continue;
      denom += (p[i] * p[i]) / var;
    }
    if (denom <= 0.0) return 0.0;
    double varA = 1.0 / denom;

    // Integral of profile shape
    double integral_p = simpson_integrate(p, tof);

    // Propagate to integrated intensity
    return varA * (integral_p * integral_p);
  }

  // Fit method: does LM and updates A,alpha,beta,sigma,T_ph
  bool fit(int maxfev = 2000, double xtol = 1e-8, double ftol = 1e-8) {
    const int ndata = static_cast<int>(tof.size());
    if (ndata < 5) {
      std::cerr << "Too few data points for fitting\n";
      return false;
    }

    // Build functor
    TOFProfileFunctor functor(tof, y_norm, min_bounds, max_bounds);

    typedef Eigen::LevenbergMarquardt<TOFProfileFunctor, double> LM;
    LM lm(functor);

    // initial parameter vector
    Eigen::VectorXd x(5);
    x[0] = A;
    x[1] = alpha;
    x[2] = beta;
    x[3] = sigma;
    x[4] = T_ph;

    // set LM parameters
    lm.parameters.maxfev = maxfev;
    lm.parameters.xtol = xtol;
    lm.parameters.ftol = ftol;

    // run
    int ret = lm.minimize(x);

    // clamp final
    for (int i = 0; i < 5; ++i) {
      x[i] = std::min(std::max(x[i], min_bounds[i]), max_bounds[i]);
    }

    // store back
    A = x[0];
    alpha = x[1];
    beta = x[2];
    sigma = x[3];
    T_ph = x[4];

    // check convergence
    if (ret == Eigen::LevenbergMarquardtSpace::ImproperInputParameters) {
      std::cerr << "LM improper input parameters\n";
      return false;
    }
    if (ret < 0) {
      std::cerr << "LM failed with code " << ret << "\n";
    }

    return true;
  }
};

#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE1D */