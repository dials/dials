#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE_1D_PROFILELIB_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE_1D_PROFILELIB_H

#include <dials/algorithms/integration/tof/tof_profile_1d_ibix.h>

/*
profile_1d_profilelib fits the same back-to-back exponential peak shape as
profile_1d_ibix, and shares profile1d_func with it, but measures that shape only
where the counts determine it and reuses it elsewhere, after Mantid's
IntegratePeaksProfileFitting.

The fitting machinery below is deliberately NOT shared with profile_1d_ibix. It
is a corrected fork: parameter bounds are imposed by reparameterisation rather
than by clamping the vector inside the residual and the Jacobian, the index the
peak-position test reads is initialised before it is used, and the acceptance
thresholds that a fixed bin count and the tallest single channel made bin-width
dependent are expressed as a time and measured against a smoothed curve.
profile_1d_ibix is left exactly as it was so that its results do not move; if
those fixes are applied to it later the two can be merged back together.
*/

namespace dials { namespace algorithms {

  /*
   * A peak shape, as fitted on one reflection and reusable on another.  Held
   * separately from TOFProfile1DProfilelibParams so that a library can be built
   * from many reflections without mutating shared fitting parameters.
   */
  struct IBIXShape {
    double alpha = 0.0;
    double beta = 0.0;
    double sigma = 0.0;
    bool valid = false;
  };

  /*
   * Everything profile_1d_profilelib needs: the shape fit, which reflections may
   * contribute a shape to the library, and the acceptance thresholds. Held
   * separately from TOFProfile1DProfilelibParams, whose fitter this one forks, so
   * that profile_1d_ibix keeps exactly the parameters it had.
   */
  struct TOFProfile1DProfilelibParams {
    double A;
    double A_min;
    double A_max;
    double alpha;
    double alpha_min;
    double alpha_max;
    double beta;
    double beta_min;
    double beta_max;
    int n_restarts;
    bool optimize_profile;
    bool show_profile_failures;
    double peak_height_smoothing;  // Time to smooth over before measuring the data peak
    double trust_min_corr;         // Correlation with the data a shape fit must reach
    double
      trust_peak_tolerance;  // How far in us the fitted peak may sit from the data's
    double trust_peak_height_fraction;  // Allowed relative error on the peak height
    double library_min_i_sigma;  // Summation I/sigma a shape's reflection must reach
    double library_min_corr;     // Correlation its free fit must reach

    TOFProfile1DProfilelibParams(double A_min,
                                 double A_max,
                                 double alpha,
                                 double alpha_min,
                                 double alpha_max,
                                 double beta,
                                 double beta_min,
                                 double beta_max,
                                 int n_restarts,
                                 bool optimize_profile,
                                 bool show_profile_failures,
                                 double peak_height_smoothing,
                                 double trust_min_corr,
                                 double trust_peak_tolerance,
                                 double trust_peak_height_fraction,
                                 double library_min_i_sigma,
                                 double library_min_corr)
        : A(A_min),
          A_min(A_min),
          A_max(A_max),
          alpha(alpha),
          alpha_min(alpha_min),
          alpha_max(alpha_max),
          beta(beta),
          beta_min(beta_min),
          beta_max(beta_max),
          n_restarts(n_restarts),
          optimize_profile(optimize_profile),
          show_profile_failures(show_profile_failures),
          peak_height_smoothing(peak_height_smoothing),
          trust_min_corr(trust_min_corr),
          trust_peak_tolerance(trust_peak_tolerance),
          trust_peak_height_fraction(trust_peak_height_fraction),
          library_min_i_sigma(library_min_i_sigma),
          library_min_corr(library_min_corr) {}
  };

  /*
   * The amplitude enters the model linearly, so for a given shape and position
   * the best A is a one-line least squares rather than something to search for:
   * A* = sum(y g) / sum(g g), with g the same profile at unit amplitude.  Used
   * to seed the fit, because until A is roughly right the model is flat
   * compared with the data and the gradient with respect to the peak position
   * all but vanishes, which leaves the optimiser unable to find the peak.
   */
  static double analytic_amplitude(scitbx::af::const_ref<double> tof,
                                   scitbx::af::const_ref<double> y_norm,
                                   double alpha,
                                   double beta,
                                   double sigma,
                                   double T_ph) {
    scitbx::af::shared<double> g = profile1d_func(tof, 1.0, alpha, beta, sigma, T_ph);
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < g.size(); ++i) {
      if (!is_finite_double(g[i])) {
        continue;
      }
      num += y_norm[i] * g[i];
      den += g[i] * g[i];
    }
    if (!(den > 0.0) || !is_finite_double(num)) {
      return 1.0;
    }
    double a = num / den;
    return is_finite_double(a) && a > 0.0 ? a : 1.0;
  }

  /*
   * Bounds by reparameterisation rather than by clamping.
   *
   * Eigen's LevenbergMarquardt is unconstrained.  Imposing bounds by clamping
   * the parameter vector inside the residual and the Jacobian, as this file
   * used to, makes the model constant outside the box: the optimiser then sees
   * no change from a step that leaves it, a parameter sitting on a bound gets a
   * zero Jacobian column, and the trust region has no direction to move in.
   * Instead every bounded parameter is carried internally as an unbounded
   * variable u, with p = lo + (hi - lo) / (1 + exp(-u)).  The map is smooth and
   * strictly monotone, every u is feasible, and the optimisation is genuinely
   * unconstrained, so the bounds cost nothing in conditioning.
   */
  inline double bounded_to_internal(double p, double lo, double hi) {
    double span = hi - lo;
    if (!(span > 0.0)) {
      return 0.0;
    }
    double frac = (p - lo) / span;
    // Keep strictly inside, or the logit is infinite
    const double eps = 1e-9;
    frac = std::min(std::max(frac, eps), 1.0 - eps);
    return std::log(frac / (1.0 - frac));
  }

  inline double internal_to_bounded(double u, double lo, double hi) {
    double span = hi - lo;
    if (!(span > 0.0)) {
      return lo;
    }
    // Logistic, written to avoid overflow for either sign of u
    double frac;
    if (u >= 0.0) {
      frac = 1.0 / (1.0 + std::exp(-u));
    } else {
      double e = std::exp(u);
      frac = e / (1.0 + e);
    }
    return lo + span * frac;
  }

  struct ProfilelibShapeFunctor {
    scitbx::af::const_ref<double> tof;
    scitbx::af::const_ref<double> y_norm;  // Assumed normalized
    std::array<double, 5> min_bounds;      // parameter bounds
    std::array<double, 5> max_bounds;      // parameter bounds
    int num_data_points, num_params;

    ProfilelibShapeFunctor(scitbx::af::const_ref<double> tof_,
                           scitbx::af::const_ref<double> y_norm_,
                           const std::array<double, 5>& minb,
                           const std::array<double, 5>& maxb)
        : tof(tof_), y_norm(y_norm_) {
      min_bounds = minb;
      max_bounds = maxb;
      num_data_points = tof.size();
      num_params = 5;
    }

    int values() const {
      return num_data_points;
    }

    int inputs() const {
      return num_params;
    }

    // The optimiser's vector is internal; these convert to and from parameters
    Eigen::VectorXd to_internal(const Eigen::VectorXd& p) const {
      Eigen::VectorXd u(num_params);
      for (int i = 0; i < num_params; ++i) {
        u[i] = bounded_to_internal(p[i], min_bounds[i], max_bounds[i]);
      }
      return u;
    }

    Eigen::VectorXd to_params(const Eigen::VectorXd& u) const {
      Eigen::VectorXd p(num_params);
      for (int i = 0; i < num_params; ++i) {
        p[i] = internal_to_bounded(u[i], min_bounds[i], max_bounds[i]);
      }
      return p;
    }

    int operator()(const Eigen::VectorXd& u, Eigen::VectorXd& fvec) const {
      Eigen::VectorXd x = to_params(u);
      scitbx::af::shared<double> model =
        profile1d_func(tof, x[0], x[1], x[2], x[3], x[4]);
      assert(model.size() == num_data_points);
      for (int i = 0; i < num_data_points; ++i) {
        fvec[i] = y_norm[i] - model[i];
      }
      return 0;
    }

    int df(const Eigen::VectorXd& u, Eigen::MatrixXd& J) const {
      // Central differences in the internal variable, where every step is
      // feasible and no clamping can flatten a column
      const double eps = 1e-5;
      J.resize(num_data_points, num_params);

      for (int j = 0; j < num_params; ++j) {
        double delta = eps * std::max(1.0, std::abs(u[j]));
        Eigen::VectorXd up = u, um = u;
        up[j] += delta;
        um[j] -= delta;

        Eigen::VectorXd fp(num_data_points), fm(num_data_points);
        operator()(up, fp);
        operator()(um, fm);
        J.col(j) = (fp - fm) / (2.0 * delta);
      }

      return 0;
    }
  };

  /*
   * As ProfilelibShapeFunctor, but with the peak shape held fixed and only the
   * amplitude and the peak position free.  Used when the shape comes from a
   * profile library built on strong reflections, as in Mantid's
   * BVGFitTools.doBVGFit with forceParams: a weak reflection has too few
   * counts to determine its own shape, but two parameters against the same
   * number of points are well constrained.
   */
  struct ProfilelibForcedFunctor {
    scitbx::af::const_ref<double> tof;
    scitbx::af::const_ref<double> y_norm;  // Assumed normalized
    double alpha, beta, sigma;             // fixed shape
    std::array<double, 2> min_bounds;      // bounds on A and T_ph
    std::array<double, 2> max_bounds;
    int num_data_points, num_params;

    ProfilelibForcedFunctor(scitbx::af::const_ref<double> tof_,
                            scitbx::af::const_ref<double> y_norm_,
                            double alpha_,
                            double beta_,
                            double sigma_,
                            const std::array<double, 2>& minb,
                            const std::array<double, 2>& maxb)
        : tof(tof_),
          y_norm(y_norm_),
          alpha(alpha_),
          beta(beta_),
          sigma(sigma_),
          min_bounds(minb),
          max_bounds(maxb) {
      num_data_points = tof.size();
      num_params = 2;
    }

    int values() const {
      return num_data_points;
    }

    int inputs() const {
      return num_params;
    }

    Eigen::VectorXd to_internal(const Eigen::VectorXd& p) const {
      Eigen::VectorXd u(num_params);
      for (int i = 0; i < num_params; ++i) {
        u[i] = bounded_to_internal(p[i], min_bounds[i], max_bounds[i]);
      }
      return u;
    }

    Eigen::VectorXd to_params(const Eigen::VectorXd& u) const {
      Eigen::VectorXd p(num_params);
      for (int i = 0; i < num_params; ++i) {
        p[i] = internal_to_bounded(u[i], min_bounds[i], max_bounds[i]);
      }
      return p;
    }

    int operator()(const Eigen::VectorXd& u, Eigen::VectorXd& fvec) const {
      Eigen::VectorXd x = to_params(u);
      scitbx::af::shared<double> model =
        profile1d_func(tof, x[0], alpha, beta, sigma, x[1]);
      for (int i = 0; i < num_data_points; ++i) {
        fvec[i] = y_norm[i] - model[i];
      }
      return 0;
    }

    int df(const Eigen::VectorXd& u, Eigen::MatrixXd& J) const {
      const double eps = 1e-5;
      J.resize(num_data_points, num_params);

      for (int j = 0; j < num_params; ++j) {
        double delta = eps * std::max(1.0, std::abs(u[j]));
        Eigen::VectorXd up = u, um = u;
        up[j] += delta;
        um[j] -= delta;

        Eigen::VectorXd fp(num_data_points), fm(num_data_points);
        operator()(up, fp);
        operator()(um, fm);
        J.col(j) = (fp - fm) / (2.0 * delta);
      }

      return 0;
    }
  };

  class TOFProfile1DProfilelibShape {
  public:
    scitbx::af::const_ref<double> tof;
    scitbx::af::const_ref<double> intensities;  // raw intensities
    scitbx::af::shared<double> y_norm;          // normalized intensities
    double intensity_max;
    int n_restarts;

    // params
    double A, alpha, beta, sigma, T_ph;
    double peak_height_smoothing;
    double trust_min_corr;
    double trust_peak_tolerance;
    double trust_peak_height_fraction;
    std::array<double, 5> min_bounds;
    std::array<double, 5> max_bounds;

    TOFProfile1DProfilelibShape(scitbx::af::const_ref<double> tof_,
                                scitbx::af::const_ref<double> intensities_,
                                double A_,
                                double alpha_,
                                double beta_,
                                double T_ph_,
                                const std::array<double, 2> A_bounds,
                                const std::array<double, 2> alpha_bounds,
                                const std::array<double, 2> beta_bounds,
                                int n_restarts_,
                                double peak_height_smoothing_,
                                double trust_min_corr_,
                                double trust_peak_tolerance_,
                                double trust_peak_height_fraction_)
        : tof(tof_),
          intensities(intensities_),
          A(A_),
          alpha(alpha_),
          beta(beta_),
          sigma(1.0),
          T_ph(T_ph_),
          peak_height_smoothing(peak_height_smoothing_),
          trust_min_corr(trust_min_corr_),
          trust_peak_tolerance(trust_peak_tolerance_),
          trust_peak_height_fraction(trust_peak_height_fraction_),
          n_restarts(n_restarts_) {
      DIALS_ASSERT(tof.size() > 0);
      DIALS_ASSERT(tof.size() == intensities.size());

      // Get max intensity
      intensity_max = 1.0;
      if (!(intensities.size() == 0)) {
        intensity_max = *std::max_element(intensities.begin(), intensities.end());
        if (intensity_max <= 0.0) intensity_max = 1.0;
      }

      // Get normalized y vector
      const size_t n = intensities.size();
      y_norm.resize(n);
      for (size_t i = 0; i < n; ++i) {
        double v = intensities[i];
        if (!is_finite_double(v)) v = 0.0;
        if (v < 0) v = 0.0;
        y_norm[i] = v / intensity_max;
      }

      // Set initial sigma from estimating peak width
      sigma = estimate_sigma_from_fwhm(tof, y_norm.const_ref());

      // Seed the amplitude at its analytic optimum for that shape
      double A_seed =
        analytic_amplitude(tof, y_norm.const_ref(), alpha, beta, sigma, T_ph);
      A = std::min(std::max(A_seed, A_bounds[0]), A_bounds[1]);

      // Param bounds (A, alpha, beta, sigma, T_ph)
      min_bounds = {
        A_bounds[0], alpha_bounds[0], beta_bounds[0], sigma / 4.0, tof.front()};

      max_bounds = {A_bounds[1],
                    alpha_bounds[1],
                    beta_bounds[1],
                    std::max(100., sigma * 4.0),
                    tof.back()};

      // Sanity check params
      DIALS_ASSERT(A >= min_bounds[0] && A <= max_bounds[0]);
      DIALS_ASSERT(alpha >= min_bounds[1] && alpha <= max_bounds[1]);
      DIALS_ASSERT(beta >= min_bounds[2] && beta <= max_bounds[2]);
      DIALS_ASSERT(sigma >= min_bounds[3] && sigma <= max_bounds[3]);
      DIALS_ASSERT(T_ph >= min_bounds[4] && T_ph <= max_bounds[4]);
    }

    scitbx::af::shared<double> result() const {
      /*
       * Generates (unnormalized) profile for all positions in tof
       */

      scitbx::af::shared<double> m = profile1d_func(tof, A, alpha, beta, sigma, T_ph);
      for (auto& v : m)
        v *= intensity_max;
      return m;
    }

    double estimate_sigma_from_fwhm(scitbx::af::const_ref<double> tof,
                                    scitbx::af::const_ref<double> y) {
      /*
       * Estimates sigma param using full width at half maximum of peak in y
       */

      // Not enough data
      if (tof.size() <= 3) {
        return 1.0;
      }

      // locate peak
      size_t imax = std::distance(y.begin(), std::max_element(y.begin(), y.end()));
      double ymax = y[imax];

      // Negative peak
      if (ymax <= 0.0) {
        return 1.0;
      }

      double half_max = 0.5 * ymax;

      // Search left crossing
      double tL = tof.front();
      for (size_t i = imax; i-- > 0;) {
        if (y[i] <= half_max && y[i + 1] > half_max) {
          double t0 = tof[i], t1 = tof[i + 1];
          double y0 = y[i], y1 = y[i + 1];
          double frac = (half_max - y0) / (y1 - y0);
          tL = t0 + frac * (t1 - t0);
          break;
        }
      }

      // Search right crossing
      double tR = tof.back();
      for (size_t i = imax; i + 1 < y.size(); ++i) {
        if (y[i] > half_max && y[i + 1] <= half_max) {
          double t0 = tof[i], t1 = tof[i + 1];
          double y0 = y[i], y1 = y[i + 1];
          double frac = (half_max - y0) / (y1 - y0);
          tR = t0 + frac * (t1 - t0);
          break;
        }
      }

      // Full width at half maximum
      double fwhm = std::max(tR - tL, 0.0);
      DIALS_ASSERT(fwhm > 0.0);

      // 2.354520045 = approx(sqrt(2ln2))
      double sigma0 = fwhm / 2.354820045;

      // Unphysical sigma (check large as at least one sample spacing)
      double mean_dt = (tof.back() - tof.front()) / std::max<size_t>(tof.size() - 1, 1);
      sigma0 = std::max(sigma0, mean_dt);
      return sigma0;
    }

    double smoothed_data_peak() const {
      const int n = static_cast<int>(y_norm.size());
      if (n == 0) {
        return 0.0;
      }
      double mean_dt = (tof.back() - tof.front()) / std::max(n - 1, 1);
      int half_window = static_cast<int>(0.5 * peak_height_smoothing / mean_dt);
      if (half_window < 1) {
        return *std::max_element(y_norm.begin(), y_norm.end());
      }
      double best = 0.0;
      for (int i = 0; i < n; ++i) {
        int lo = std::max(0, i - half_window);
        int hi = std::min(n - 1, i + half_window);
        double total = 0.0;
        for (int j = lo; j <= hi; ++j) {
          total += y_norm[j];
        }
        best = std::max(best, total / (hi - lo + 1));
      }
      return best;
    }

    double calc_intensity() const {
      /**
       * Get overall intensity with Simpsons rule then divide by mean_dt to
       * approximate summation scale
       */

      scitbx::af::shared<double> r = result();
      double mean_dt = (tof[tof.size() - 1] - tof[0]) / (tof.size() - 1);
      return simpson_integrate(r.const_ref(), tof) / mean_dt;
    }

    std::size_t get_max_profile_index() {
      /*
       * Returns the index of the max of the profile
       */

      auto profile_result = this->result();
      auto max_profile_it =
        std::max_element(profile_result.begin(), profile_result.end());
      std::size_t max_profile_index =
        std::distance(profile_result.begin(), max_profile_it);
      return max_profile_index;
    }

    bool fit(std::size_t max_sum_index,  // Peak index of the projected intensity
             bool show_profile_failures,
             int maxfev = 200,
             double xtol = 1e-8,
             double ftol = 1e-8) {
      /*
       * Least-squares minimization
       * Updates A, alpha, beta, sigma, T_ph
       * If fitting fails, params are perturbed n_restarts to find a solution
       */

      // Check enough data for fitting
      const int ndata = static_cast<int>(tof.size());
      if (ndata < 5) return false;

      ProfilelibShapeFunctor functor(tof, y_norm.const_ref(), min_bounds, max_bounds);
      typedef Eigen::LevenbergMarquardt<ProfilelibShapeFunctor, double> LM;

      auto run_single_fit = [&](const Eigen::VectorXd& x_init,
                                double& final_error) -> bool {
        LM lm(functor);
        lm.parameters.maxfev = maxfev;
        lm.parameters.xtol = xtol;
        lm.parameters.ftol = ftol;

        Eigen::VectorXd u = functor.to_internal(x_init);
        int result = lm.minimize(u);
        if (result < 0) return false;

        Eigen::VectorXd x = functor.to_params(u);

        // Compute residual norm (the functor takes internal variables)
        Eigen::VectorXd fvec(functor.num_data_points);
        functor(u, fvec);
        final_error = fvec.squaredNorm();

        // Update fitted parameters
        A = x[0];
        alpha = x[1];
        beta = x[2];
        sigma = x[3];
        T_ph = x[4];

        return true;
      };

      // First fit attempt
      Eigen::VectorXd x0(5);
      x0 << A, alpha, beta, sigma, T_ph;
      double fit_resid = std::numeric_limits<double>::infinity();
      bool success = run_single_fit(x0, fit_resid);
      std::size_t max_profile_index;
      double I_prf, I_var;

      if (success) {
        I_prf = this->calc_intensity();
        max_profile_index = this->get_max_profile_index();
        if (this->trust_result(fit_resid,
                               I_prf,
                               max_sum_index,
                               max_profile_index,
                               show_profile_failures)) {
          return true;
        }
      }

      // Initial fit failed, perturb initial params
      std::mt19937 rng(std::random_device{}());
      std::uniform_real_distribution<double> unit_dist(0.0, 1.0);

      for (int i = 0; i < n_restarts; ++i) {
        Eigen::VectorXd x_try(5);
        for (int j = 0; j < 5; ++j) {
          double span = max_bounds[j] - min_bounds[j];
          double rand_frac = (unit_dist(rng) - 0.5) * 0.4;
          double perturbed = x0[j] + rand_frac * span;
          x_try[j] = std::max(min_bounds[j], std::min(perturbed, max_bounds[j]));
        }

        // Attempt fit
        success = run_single_fit(x_try, fit_resid);
        if (!success) continue;

        I_prf = this->calc_intensity();
        max_profile_index = this->get_max_profile_index();
        if (this->trust_result(fit_resid,
                               I_prf,
                               max_sum_index,
                               max_profile_index,
                               show_profile_failures)) {
          return true;
        }
      }

      return false;
    }

    bool trust_result(double error,
                      double I_prf,
                      std::size_t max_sum_index,
                      std::size_t max_profile_index,
                      bool show_error = false) {
      /*
       * Tests to check the fit is reasonable
       */
      if (!std::isfinite(error) || error <= 0.0) {
        if (show_error) {
          std::cerr << "profile1d fitting failure: invalid error value (error=" << error
                    << ")\n";
        }
        return false;
      }

      // Check reasonable intensity
      if (I_prf < 1e-7) {
        if (show_error) {
          std::cerr << "profile1d fitting failure: profile intensity too small (I_prf="
                    << I_prf << ")\n";
        }
        return false;
      }

      /*
       * Check the peak position is close to the data peak.  With a tolerance in
       * microseconds the comparison is made as a time, which does not shrink
       * with the slice width; with none it falls back to the original fixed
       * three bins.
       */
      if (trust_peak_tolerance > 0.0) {
        double peak_delta = std::abs(tof[max_profile_index] - tof[max_sum_index]);
        if (peak_delta > trust_peak_tolerance) {
          if (show_error) {
            std::cerr << "profile1d fitting failure: peak position mismatch (delta="
                      << peak_delta << " > " << trust_peak_tolerance << ")\n";
          }
          return false;
        }
      } else {
        int peak_delta = std::abs(static_cast<int>(max_sum_index)
                                  - static_cast<int>(max_profile_index));
        if (peak_delta > 3) {
          if (show_error) {
            std::cerr << "profile1d fitting failure: peak index mismatch (delta="
                      << peak_delta << ")\n";
          }
          return false;
        }
      }

      // Check peak isn't very flat
      auto m = result();
      double max_val = *std::max_element(m.begin(), m.end());
      double mean_val = std::accumulate(m.begin(), m.end(), 0.0) / m.size();
      double contrast = (max_val - mean_val) / (max_val + 1e-12);
      if (contrast < 0.1) {
        if (show_error) {
          std::cerr << "profile1d fitting failure: insufficient peak contrast "
                    << "(contrast=" << contrast << ", max_val=" << max_val
                    << ", mean_val=" << mean_val << ")\n";
        }
        return false;
      }

      // Check correlation with data
      double profile_peak = 0.0, data_peak = 0.0;
      double num = 0.0, denom_y = 0.0, denom_m = 0.0;

      for (std::size_t i = 0; i < tof.size(); ++i) {
        double y = y_norm[i];
        double p = m[i] / intensity_max;

        if (i == 0 || y > data_peak) {
          data_peak = y;
        }
        if (i == 0 || p > profile_peak) {
          profile_peak = p;
        }

        num += y * p;
        denom_y += y * y;
        denom_m += p * p;
      }

      double corr = num / std::sqrt(denom_y * denom_m + 1e-12);
      if (corr < trust_min_corr) {
        if (show_error) {
          std::cerr << "profile1d fitting failure: low correlation (corr=" << corr
                    << ")\n";
        }
        return false;
      }

      /*
       * Check the peak height against the data, measured on a smoothed copy so
       * that the comparison is with the curve rather than with its noisiest
       * channel, and so that it does not tighten as the slices get finer.
       */
      double reference_peak = data_peak;
      if (peak_height_smoothing > 0.0) {
        reference_peak = std::max(smoothed_data_peak(), 1e-12);
      }
      double peak_diff = std::abs(profile_peak - reference_peak);
      if (peak_diff > reference_peak * trust_peak_height_fraction) {
        if (show_error) {
          std::cerr << "profile1d fitting failure: peak height mismatch "
                    << "(profile_peak=" << profile_peak
                    << ", data_peak=" << reference_peak
                    << ", raw_data_peak=" << data_peak << ", diff=" << peak_diff
                    << ")\n";
        }
        return false;
      }

      return true;
    }
  };

  bool fit_profile_1d_profilelib_shape(
    scitbx::af::const_ref<double> projected_intensity,
    scitbx::af::const_ref<double> tof_z,
    TOFProfile1DProfilelibParams& profile_params,
    double& I_prf_out,
    boost::optional<scitbx::af::shared<double>> line_profile_out = boost::none,
    bool update_params = false,
    IBIXShape* shape_out = nullptr) {
    /**
     * Wrapper for fitting a given reflection
     * If line_profile_out is provided the profile is returned at every
     * position in tof_z
     */

    // Get T_ph (peak position)
    auto max_it =
      std::max_element(projected_intensity.begin(), projected_intensity.end());
    size_t max_index = std::distance(projected_intensity.begin(), max_it);
    double T_ph = tof_z[max_index];

    // Fit profile
    const std::array<double, 2> A_bounds = {profile_params.A_min, profile_params.A_max};
    const std::array<double, 2> alpha_bounds = {profile_params.alpha_min,
                                                profile_params.alpha_max};
    const std::array<double, 2> beta_bounds = {profile_params.beta_min,
                                               profile_params.beta_max};

    TOFProfile1DProfilelibShape profile(tof_z,
                                        projected_intensity,
                                        profile_params.A,
                                        profile_params.alpha,
                                        profile_params.beta,
                                        T_ph,
                                        A_bounds,
                                        alpha_bounds,
                                        beta_bounds,
                                        profile_params.n_restarts,
                                        profile_params.peak_height_smoothing,
                                        profile_params.trust_min_corr,
                                        profile_params.trust_peak_tolerance,
                                        profile_params.trust_peak_height_fraction);

    bool profile_success = true;
    if (profile_params.optimize_profile) {
      profile_success = profile.fit(max_index, profile_params.show_profile_failures);
    }

    if (profile_success) {
      if (update_params) {
        profile_params.alpha = profile.alpha;
        profile_params.beta = profile.beta;
        profile_params.A = profile.A;
      }
      if (shape_out != nullptr) {
        shape_out->alpha = profile.alpha;
        shape_out->beta = profile.beta;
        shape_out->sigma = profile.sigma;
        shape_out->valid = true;
      }
      double I_prf = profile.calc_intensity();
      auto profile_result = profile.result();
      DIALS_ASSERT(projected_intensity.size() == profile_result.size());

      I_prf_out = I_prf;

      if (!line_profile_out) {
        return profile_success;
      }

      scitbx::af::shared<double> line_profile = *line_profile_out;
      DIALS_ASSERT(line_profile.size() == profile_result.size());
      for (std::size_t i = 0; i < profile_result.size(); ++i) {
        line_profile[i] = profile_result[i];
      }
      return profile_success;
    }
    return false;
  }

  /*
   * Fit a reflection with its peak shape taken from elsewhere, leaving only
   * the amplitude and the peak position free.
   *
   * Deliberately not gated on the correlation between model and data, unlike
   * the free fit: a forced fit has fewer parameters and so follows the noise
   * less closely, which lowers that correlation while raising the accuracy of
   * the intensity.  Gating on it would reject the better answer.
   */
  bool fit_profile_1d_profilelib_forced(
    scitbx::af::const_ref<double> projected_intensity,
    scitbx::af::const_ref<double> tof_z,
    const IBIXShape& shape,
    double A_min,
    double A_max,
    double& I_prf_out,
    boost::optional<scitbx::af::shared<double>> line_profile_out = boost::none) {
    const int ndata = static_cast<int>(tof_z.size());
    if (ndata < 5 || !shape.valid) {
      return false;
    }
    if (!(shape.sigma > 0.0) || !(shape.alpha > 0.0) || !(shape.beta > 0.0)) {
      return false;
    }

    // Normalize, as the free fitter does, so that A is on a comparable scale
    double intensity_max =
      *std::max_element(projected_intensity.begin(), projected_intensity.end());
    if (!(intensity_max > 0.0)) {
      return false;
    }
    scitbx::af::shared<double> y_norm(ndata);
    for (int i = 0; i < ndata; ++i) {
      double v = projected_intensity[i];
      if (!is_finite_double(v) || v < 0.0) v = 0.0;
      y_norm[i] = v / intensity_max;
    }

    // Seed the peak position at the tallest channel
    std::size_t max_index = std::distance(
      projected_intensity.begin(),
      std::max_element(projected_intensity.begin(), projected_intensity.end()));

    std::array<double, 2> min_bounds = {A_min, tof_z.front()};
    std::array<double, 2> max_bounds = {A_max, tof_z.back()};

    ProfilelibForcedFunctor functor(tof_z,
                                    y_norm.const_ref(),
                                    shape.alpha,
                                    shape.beta,
                                    shape.sigma,
                                    min_bounds,
                                    max_bounds);

    double A_seed = analytic_amplitude(tof_z,
                                       y_norm.const_ref(),
                                       shape.alpha,
                                       shape.beta,
                                       shape.sigma,
                                       tof_z[max_index]);
    Eigen::VectorXd x(2);
    x << std::min(std::max(A_seed, A_min), A_max), tof_z[max_index];

    Eigen::LevenbergMarquardt<ProfilelibForcedFunctor, double> lm(functor);
    lm.parameters.maxfev = 200;
    lm.parameters.xtol = 1e-8;
    lm.parameters.ftol = 1e-8;
    Eigen::VectorXd u = functor.to_internal(x);
    if (lm.minimize(u) < 0) {
      return false;
    }
    x = functor.to_params(u);

    scitbx::af::shared<double> model =
      profile1d_func(tof_z, x[0], shape.alpha, shape.beta, shape.sigma, x[1]);
    for (auto& v : model) {
      v *= intensity_max;
      if (!is_finite_double(v)) {
        return false;
      }
    }

    double mean_dt = (tof_z[ndata - 1] - tof_z[0]) / (ndata - 1);
    if (!(mean_dt > 0.0)) {
      return false;
    }
    double I_prf = simpson_integrate(model.const_ref(), tof_z) / mean_dt;
    if (!is_finite_double(I_prf) || I_prf < 1e-7) {
      return false;
    }

    I_prf_out = I_prf;
    if (line_profile_out) {
      scitbx::af::shared<double> line_profile = *line_profile_out;
      if (line_profile.size() != model.size()) {
        return false;
      }
      for (std::size_t i = 0; i < model.size(); ++i) {
        line_profile[i] = model[i];
      }
    }
    return true;
  }

}}  // namespace dials::algorithms
#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE_1D_PROFILELIB_H */
