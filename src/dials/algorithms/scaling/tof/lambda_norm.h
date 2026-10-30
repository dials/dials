#ifndef DIALS_ALGORITHMS_SCALING_TOF_LAMBDA_NORM_H
#define DIALS_ALGORITHMS_SCALING_TOF_LAMBDA_NORM_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/constants.h>
#include <dials/error.h>

#include <eigen3/Eigen/Dense>

/*
Gaussian-process wavelength normalisation for Laue/ToF data, after the
LAUE NORMGPR method of lawless (the Laue branch of aimless).
*/
namespace dials_scaling {

using scitbx::af::const_ref;
using scitbx::af::shared;

inline shared<double> to_shared(const std::vector<double>& v) {
  return shared<double>(v.data(), v.data() + v.size());
}

/**
 * sigma_exp^2 = max(sigma^2 - b*I, 0) + b*max(I_pred, 0), with b(lambda)
 * piecewise constant over equal-count wavelength bins.
 */
class LambdaNormVarianceModel {
public:
  LambdaNormVarianceModel() {}

  void fit(const const_ref<double>& lam,
           const const_ref<double>& intensity,
           const const_ref<double>& variance) {
    DIALS_ASSERT(lam.size() == intensity.size());
    DIALS_ASSERT(lam.size() == variance.size());
    const std::size_t n_lambda_bins = 8;
    const std::size_t n_intensity_bins = 10;
    struct IV {
      double lam, I, v;
    };
    std::vector<IV> data;
    for (std::size_t i = 0; i < lam.size(); ++i) {
      if (intensity[i] > 0.0 && variance[i] > 0.0) {
        IV d = {lam[i], intensity[i], variance[i]};
        data.push_back(d);
      }
    }
    lambda_edges_.clear();
    b_.clear();
    share_.clear();
    const std::size_t N = data.size();
    if (N < n_lambda_bins * n_intensity_bins * 5) return;
    std::sort(
      data.begin(), data.end(), [](const IV& a, const IV& b) { return a.lam < b.lam; });
    for (std::size_t k = 0; k < n_lambda_bins; ++k) {
      std::size_t i0 = (k * N) / n_lambda_bins;
      std::size_t i1 = ((k + 1) * N) / n_lambda_bins;
      std::vector<IV> bin(data.begin() + i0, data.begin() + i1);
      std::sort(
        bin.begin(), bin.end(), [](const IV& a, const IV& b) { return a.I < b.I; });
      const std::size_t n = bin.size();
      // straight line sigma^2 = a + b*I through per-quantile medians
      double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
      for (std::size_t q = 0; q < n_intensity_bins; ++q) {
        std::size_t j0 = (q * n) / n_intensity_bins;
        std::size_t j1 = ((q + 1) * n) / n_intensity_bins;
        std::vector<double> vv;
        for (std::size_t j = j0; j < j1; ++j)
          vv.push_back(bin[j].v);
        std::nth_element(vv.begin(), vv.begin() + vv.size() / 2, vv.end());
        double mI = bin[(j0 + j1) / 2].I;
        double mV = vv[vv.size() / 2];
        sx += mI;
        sy += mV;
        sxx += mI * mI;
        sxy += mI * mV;
      }
      double nq = double(n_intensity_bins);
      double den = nq * sxx - sx * sx;
      double b = (den > 0.0) ? (nq * sxy - sx * sy) / den : 0.0;
      if (!(b > 0.0)) b = 0.0;
      std::vector<double> share;
      for (std::size_t j = 0; j < n; ++j) {
        share.push_back(std::min(1.0, b * bin[j].I / bin[j].v));
      }
      std::nth_element(share.begin(), share.begin() + n / 2, share.end());
      lambda_edges_.push_back(data[i0].lam);
      b_.push_back(b);
      share_.push_back(share[n / 2]);
    }
    lambda_edges_.push_back(data[N - 1].lam);
  }

  double b(double lambda) const {
    if (b_.empty()) return 0.0;
    std::ptrdiff_t k =
      std::upper_bound(lambda_edges_.begin(), lambda_edges_.end() - 1, lambda)
      - lambda_edges_.begin() - 1;
    k = std::max<std::ptrdiff_t>(0, std::min<std::ptrdiff_t>(k, b_.size() - 1));
    return b_[k];
  }

  double variance(double lambda, double I, double var, double I_pred) const {
    double b_lam = b(lambda);
    double v = std::max(var - b_lam * I, 0.0) + b_lam * std::max(I_pred, 0.0);
    return std::max(v, 0.01 * var);
  }

  shared<double> lambda_edges() const {
    return to_shared(lambda_edges_);
  }
  shared<double> b_values() const {
    return to_shared(b_);
  }
  /// Median fraction b*I/sigma^2 of the variance that is counting noise.
  shared<double> counting_share() const {
    return to_shared(share_);
  }

private:
  std::vector<double> lambda_edges_;  // n_bins + 1
  std::vector<double> b_;
  std::vector<double> share_;
};

/// One leave-one-out ratio per usable observation.
struct LambdaNormSamples {
  shared<std::size_t> index;
  shared<double> rho;  // I_i / M, M the mates' mean with correction divided out
  shared<double> u;    // ratio-estimator weight, 1 / (M var(rho))
};

/**
 * Leave-one-out samples over groups of symmetry equivalents.
 *
 * Observations are sorted by group; group g is
 * [group_offsets[g], group_offsets[g+1]). correction[i] is the factor
 * (w(lambda) times any per-experiment factor) taken out of observation i
 * when it serves as a mate.
 */
inline LambdaNormSamples leave_one_out_samples(
  const const_ref<double>& intensity,
  const const_ref<double>& variance,
  const const_ref<double>& lam,
  const const_ref<double>& correction,
  const const_ref<std::size_t>& group_offsets,
  const LambdaNormVarianceModel& vmodel,
  bool expected_sigma) {
  const std::size_t N = intensity.size();
  DIALS_ASSERT(variance.size() == N);
  DIALS_ASSERT(lam.size() == N);
  DIALS_ASSERT(correction.size() == N);
  DIALS_ASSERT(group_offsets.size() >= 1);
  DIALS_ASSERT(group_offsets[group_offsets.size() - 1] == N);
  LambdaNormSamples out;
  for (std::size_t g = 0; g + 1 < group_offsets.size(); ++g) {
    const std::size_t first = group_offsets[g];
    const std::size_t last = group_offsets[g + 1];
    DIALS_ASSERT(first <= last);
    if (last - first < 2) continue;
    for (std::size_t i = first; i < last; ++i) {
      double sw = 0.0, swi = 0.0;
      for (std::size_t j = first; j < last; ++j) {
        if (j == i) continue;
        double c = correction[j];
        double w = c * c / variance[j];
        sw += w;
        swi += w * intensity[j] / c;
      }
      if (!(sw > 0.0)) continue;
      double M = swi / sw;
      if (!(M > 0.0)) continue;
      double ve = variance[i];
      if (expected_sigma) {
        // reweight the mates by their variance at the predicted intensity
        sw = 0.0;
        swi = 0.0;
        for (std::size_t j = first; j < last; ++j) {
          if (j == i) continue;
          double c = correction[j];
          double w = c * c / vmodel.variance(lam[j], intensity[j], variance[j], c * M);
          sw += w;
          swi += w * intensity[j] / c;
        }
        if (!(sw > 0.0)) continue;
        M = swi / sw;
        if (!(M > 0.0)) continue;
        ve = vmodel.variance(lam[i], intensity[i], variance[i], correction[i] * M);
      }
      out.index.push_back(i);
      out.rho.push_back(intensity[i] / M);
      out.u.push_back(M / ve);
    }
  }
  return out;
}

/**
 * Smooth w(lambda) from leave-one-out ratios: equal-count bins, ratio
 * estimator per bin, GP in log(intensity) over x = log(lambda) with
 * hyperparameters by grid-searched marginal likelihood.
 */
class GPRWavelengthScale {
public:
  /// lengthscale in Angstrom at lambda_ref (<= 0: optimise); lam_min,
  /// lam_max <= 0 take the range from the data.
  GPRWavelengthScale(std::size_t nbins,
                     double lengthscale,
                     bool matern32,
                     double lam_min,
                     double lam_max)
      : nbins_(nbins),
        fixed_lengthscale_(lengthscale),
        matern32_(matern32),
        lam_min_req_(lam_min),
        lam_max_req_(lam_max),
        active_(false),
        lam_min_(0),
        lam_max_(0),
        lambda_ref_(0),
        g_ref_(0),
        gmean_(0),
        ell_(0),
        sigf_(0),
        noise_inflation_(1),
        mean_noise_(0),
        lml_(0),
        n_samples_(0),
        ell_at_lower_(false),
        ell_at_upper_(false),
        noise_at_upper_(false),
        grid_step_(0) {}

  /// Returns false, with message() set, if there is too little data.
  bool fit(const const_ref<double>& lam,
           const const_ref<double>& rho,
           const const_ref<double>& u,
           double lambda_ref) {
    DIALS_ASSERT(lam.size() == rho.size());
    DIALS_ASSERT(lam.size() == u.size());
    reset();
    const std::size_t min_samples = 50;
    const std::size_t min_count = 20;

    lam_min_ = lam_min_req_;
    lam_max_ = lam_max_req_;
    if (lam_min_ <= 0.0 || lam_max_ <= 0.0 || lam_max_ <= lam_min_) {
      lam_min_ = std::numeric_limits<double>::max();
      lam_max_ = -std::numeric_limits<double>::max();
      for (std::size_t n = 0; n < lam.size(); ++n) {
        lam_min_ = std::min(lam_min_, lam[n]);
        lam_max_ = std::max(lam_max_, lam[n]);
      }
    }
    if (!(lam_max_ > lam_min_) || !(lam_min_ > 0.0)) {
      message_ = "degenerate wavelength range";
      return false;
    }

    std::vector<Sample> smp;
    for (std::size_t n = 0; n < lam.size(); ++n) {
      if (lam[n] < lam_min_ || lam[n] > lam_max_) continue;
      if (!(u[n] > 0.0) || !std::isfinite(u[n]) || !std::isfinite(rho[n])) continue;
      Sample s = {lam[n], rho[n], u[n]};
      smp.push_back(s);
    }
    const std::size_t N = smp.size();
    n_samples_ = N;
    if (N < min_samples) {
      message_ = "too few usable samples (" + std::to_string(N) + ")";
      return false;
    }
    std::sort(smp.begin(), smp.end(), [](const Sample& a, const Sample& b) {
      return a.lam < b.lam;
    });

    // Equal-count bins, closed early at a maximum width in log(lambda) so
    // the sparse tails still get training points near the ends of the range.
    std::size_t nbins = nbins_ > 0 ? nbins_ : 50;
    nbins = std::min(nbins, N / min_count);
    if (nbins < 4) {
      message_ = "too few samples for binning (" + std::to_string(N) + ")";
      return false;
    }
    const double xspan = std::log(smp[N - 1].lam) - std::log(smp[0].lam);
    const double target = double(N) / nbins;
    const double wcap = 4.0 * xspan / nbins;
    const std::size_t mincnt =
      std::max(min_count, static_cast<std::size_t>(0.25 * target));
    std::vector<double> tl, ty, tn;
    std::size_t i0 = 0;
    while (i0 < N) {
      double x0 = std::log(smp[i0].lam);
      std::size_t i1 = i0 + 1;
      while (i1 < N) {
        std::size_t cnt = i1 - i0;
        if (cnt >= mincnt
            && (cnt >= target || std::log(smp[i1 - 1].lam) - x0 >= wcap)) {
          break;
        }
        ++i1;
      }
      if (N - i1 < mincnt) i1 = N;
      add_bin(smp, i0, i1, min_count, tl, ty, tn);
      i0 = i1;
    }
    const std::size_t B = tl.size();
    if (B < 4) {
      message_ = "too few usable bins (" + std::to_string(B) + ")";
      return false;
    }

    std::vector<double> tnsort(tn);
    std::sort(tnsort.begin(), tnsort.end());
    double nfloor = std::max(1.0e-6, 0.01 * tnsort[B / 2]);
    for (std::size_t i = 0; i < B; ++i)
      tn[i] = std::max(tn[i], nfloor);

    std::vector<double> tx(B);
    for (std::size_t i = 0; i < B; ++i)
      tx[i] = std::log(tl[i]);
    const double xrange = tx[B - 1] - tx[0];
    if (!(xrange > 0.0)) {
      message_ = "degenerate training range";
      return false;
    }
    std::vector<double> dxs;
    for (std::size_t i = 1; i < B; ++i)
      dxs.push_back(std::fabs(tx[i] - tx[i - 1]));
    std::sort(dxs.begin(), dxs.end());
    const double dxmed = dxs[dxs.size() / 2];

    // constant mean: away from the data the posterior reverts to it
    double swm = 0.0, swmy = 0.0;
    for (std::size_t i = 0; i < B; ++i) {
      swm += 1.0 / tn[i];
      swmy += ty[i] / tn[i];
    }
    gmean_ = swmy / swm;
    Eigen::VectorXd Y(B);
    double yvar = 0.0, nmean = 0.0;
    for (std::size_t i = 0; i < B; ++i) {
      Y(i) = ty[i] - gmean_;
      yvar += Y(i) * Y(i);
      nmean += tn[i];
    }
    yvar /= B;
    nmean /= B;
    const double sig0 = std::sqrt(std::max(yvar - nmean, 1.0e-4));

    lambda_ref_ = lambda_ref;
    std::vector<double> ellcand;
    const double elllo = std::max(2.0 * dxmed, 0.02 * xrange);
    const double ellhi = 1.5 * xrange;
    const bool ell_fixed = fixed_lengthscale_ > 0.0;
    if (ell_fixed) {
      double lscale = lambda_ref_ > 0.0 ? lambda_ref_ : std::sqrt(lam_min_ * lam_max_);
      ellcand.push_back(fixed_lengthscale_ / lscale);
    } else {
      const int n_ell = 16;
      for (int k = 0; k < n_ell; ++k) {
        double f = double(k) / (n_ell - 1);
        ellcand.push_back(elllo * std::pow(ellhi / elllo, f));
      }
    }
    std::vector<double> sigcand, alfcand;
    const int n_sigf = 9;
    for (int k = 0; k < n_sigf; ++k) {
      sigcand.push_back(sig0 * std::pow(10.0, -0.7 + 1.4 * double(k) / (n_sigf - 1)));
    }
    const int n_alpha = 8;
    for (int k = 0; k < n_alpha; ++k) {
      alfcand.push_back(std::pow(10.0, -0.3 + 1.5 * double(k) / (n_alpha - 1)));
    }

    const double jitter = 1.0e-10;
    const double log2pi = std::log(2.0 * scitbx::constants::pi);
    double best_lml = -std::numeric_limits<double>::max();
    double best_ell = ellcand[0], best_sig = sigcand[0], best_alf = alfcand[0];
    Eigen::VectorXd best_alpha;
    Eigen::MatrixXd K0(B, B);
    for (std::size_t e = 0; e < ellcand.size(); ++e) {
      for (std::size_t i = 0; i < B; ++i) {
        for (std::size_t j = 0; j < B; ++j) {
          K0(i, j) = kernel(tx[i], tx[j], ellcand[e], 1.0);
        }
      }
      for (std::size_t s = 0; s < sigcand.size(); ++s) {
        const double sf2 = sigcand[s] * sigcand[s];
        for (std::size_t a = 0; a < alfcand.size(); ++a) {
          Eigen::MatrixXd A = sf2 * K0;
          for (std::size_t i = 0; i < B; ++i)
            A(i, i) += alfcand[a] * tn[i] + jitter;
          Eigen::LLT<Eigen::MatrixXd> llt(A);
          if (llt.info() != Eigen::Success) continue;
          Eigen::VectorXd alpha = llt.solve(Y);
          double logdet = 0.0;
          Eigen::MatrixXd L = llt.matrixL();
          for (std::size_t i = 0; i < B; ++i)
            logdet += std::log(L(i, i));
          double lml = -0.5 * Y.dot(alpha) - logdet - 0.5 * B * log2pi;
          if (lml > best_lml) {
            best_lml = lml;
            best_ell = ellcand[e];
            best_sig = sigcand[s];
            best_alf = alfcand[a];
            best_alpha = alpha;
          }
        }
      }
    }
    if (best_alpha.size() == 0) {
      message_ = "Cholesky factorisation failed";
      return false;
    }
    ell_ = best_ell;
    sigf_ = best_sig;
    noise_inflation_ = best_alf;
    lml_ = best_lml;
    ell_at_lower_ = !ell_fixed && ell_ < 1.01 * elllo;
    ell_at_upper_ = !ell_fixed && ell_ > 0.99 * ellhi;
    noise_at_upper_ = noise_inflation_ > 0.99 * alfcand.back();

    Eigen::MatrixXd A(B, B);
    for (std::size_t i = 0; i < B; ++i) {
      for (std::size_t j = 0; j < B; ++j)
        A(i, j) = kernel(tx[i], tx[j]);
      A(i, i) += noise_inflation_ * tn[i] + jitter;
    }
    Eigen::LLT<Eigen::MatrixXd> llt(A);

    // Posterior mean held at its edge value beyond the outermost bins, so
    // the fit cannot run away where there is no data; the SD is not held.
    const double xlo = tx[0], xhi = tx[B - 1];
    const std::size_t ngrid = std::max<std::size_t>(200, 4 * nbins);
    grid_step_ = (lam_max_ - lam_min_) / (ngrid - 1);
    grid_g_.assign(ngrid, 0.0);
    grid_sd_.assign(ngrid, 0.0);
    Eigen::VectorXd kstar(B);
    for (std::size_t m = 0; m < ngrid; ++m) {
      double x = std::log(lam_min_ + m * grid_step_);
      grid_g_[m] = posterior_mean(std::min(std::max(x, xlo), xhi), tx, best_alpha);
      for (std::size_t i = 0; i < B; ++i)
        kstar(i) = kernel(x, tx[i]);
      Eigen::VectorXd v = llt.matrixL().solve(kstar);
      grid_sd_[m] = std::sqrt(std::max(kernel(x, x) - v.dot(v), 0.0));
    }
    double xref =
      std::log(lambda_ref_ > 0.0 ? lambda_ref_ : std::sqrt(lam_min_ * lam_max_));
    g_ref_ = posterior_mean(std::min(std::max(xref, xlo), xhi), tx, best_alpha);

    double nsum = 0.0;
    train_lam_ = tl;
    train_h_ = ty;
    train_sd_.resize(B);
    for (std::size_t i = 0; i < B; ++i) {
      train_sd_[i] = std::sqrt(noise_inflation_ * tn[i]);
      nsum += train_sd_[i];
    }
    mean_noise_ = nsum / B;
    active_ = true;
    return true;
  }

  /// w(lambda) = exp(g(lambda) - g(lambda_ref)); held at the range ends.
  shared<double> scale(const const_ref<double>& lam) const {
    shared<double> out(lam.size(), 1.0);
    if (!active_) return out;
    for (std::size_t i = 0; i < lam.size(); ++i) {
      out[i] = std::exp(interpolate(grid_g_, lam[i]) - g_ref_);
    }
    return out;
  }

  /// Posterior SD of log w, i.e. the relative uncertainty of w.
  shared<double> uncertainty(const const_ref<double>& lam) const {
    shared<double> out(lam.size(), 0.0);
    if (!active_) return out;
    for (std::size_t i = 0; i < lam.size(); ++i) {
      out[i] = interpolate(grid_sd_, lam[i]);
    }
    return out;
  }

  bool active() const {
    return active_;
  }
  std::string message() const {
    return message_;
  }
  double lambda_min() const {
    return lam_min_;
  }
  double lambda_max() const {
    return lam_max_;
  }
  double lambda_ref() const {
    return lambda_ref_;
  }
  /// In log(lambda) units; multiply by lambda for Angstrom.
  double lengthscale() const {
    return ell_;
  }
  double sigma_f() const {
    return sigf_;
  }
  double noise_inflation() const {
    return noise_inflation_;
  }
  double mean_noise() const {
    return mean_noise_;
  }
  double log_marginal_likelihood() const {
    return lml_;
  }
  std::size_t n_samples() const {
    return n_samples_;
  }
  bool lengthscale_at_lower_limit() const {
    return ell_at_lower_;
  }
  bool lengthscale_at_upper_limit() const {
    return ell_at_upper_;
  }
  bool noise_inflation_at_upper_limit() const {
    return noise_at_upper_;
  }
  shared<double> training_lambda() const {
    return to_shared(train_lam_);
  }
  /// Training bins as w, relative to g(lambda_ref) like scale().
  shared<double> training_scale() const {
    shared<double> out;
    for (std::size_t i = 0; i < train_h_.size(); ++i) {
      out.push_back(std::exp(train_h_[i] - g_ref_));
    }
    return out;
  }
  /// Relative SD of each training bin, noise inflation included.
  shared<double> training_sigma() const {
    return to_shared(train_sd_);
  }

private:
  struct Sample {
    double lam, rho, u;
  };

  void reset() {
    active_ = false;
    message_.clear();
    grid_g_.clear();
    grid_sd_.clear();
    train_lam_.clear();
    train_h_.clear();
    train_sd_.clear();
    n_samples_ = 0;
    ell_at_lower_ = ell_at_upper_ = noise_at_upper_ = false;
  }

  // Ratio estimator r = sum(u rho)/sum(u) with sandwich variance; the
  // training target is log r. A mean of per-observation log ratios is
  // biased and heavy-tailed for weak data.
  static void add_bin(const std::vector<Sample>& smp,
                      std::size_t i0,
                      std::size_t i1,
                      std::size_t min_count,
                      std::vector<double>& tl,
                      std::vector<double>& ty,
                      std::vector<double>& tn) {
    if (i1 - i0 < min_count) return;
    double su = 0.0, sur = 0.0, sul = 0.0;
    for (std::size_t i = i0; i < i1; ++i) {
      su += smp[i].u;
      sur += smp[i].u * smp[i].rho;
      sul += smp[i].u * smp[i].lam;
    }
    if (!(su > 0.0)) return;
    double r = sur / su;
    if (!(r > 0.0)) return;
    double sv = 0.0;
    for (std::size_t i = i0; i < i1; ++i) {
      double e = smp[i].u * (smp[i].rho - r);
      sv += e * e;
    }
    double varh = sv / (su * su) / (r * r);
    if (!std::isfinite(varh)) return;
    tl.push_back(sul / su);
    ty.push_back(std::log(r));
    tn.push_back(varh);
  }

  double kernel(double xa, double xb, double ell, double sigf) const {
    double r = (xa - xb) / ell;
    if (matern32_) {
      double d = std::sqrt(3.0) * std::fabs(r);
      return sigf * sigf * (1.0 + d) * std::exp(-d);
    }
    return sigf * sigf * std::exp(-0.5 * r * r);
  }

  double kernel(double xa, double xb) const {
    return kernel(xa, xb, ell_, sigf_);
  }

  double posterior_mean(double x,
                        const std::vector<double>& tx,
                        const Eigen::VectorXd& alpha) const {
    double g = gmean_;
    for (std::size_t i = 0; i < tx.size(); ++i)
      g += kernel(x, tx[i]) * alpha(i);
    return g;
  }

  double interpolate(const std::vector<double>& grid, double lambda) const {
    const std::size_t n = grid.size();
    double t = (std::min(std::max(lambda, lam_min_), lam_max_) - lam_min_) / grid_step_;
    std::size_t i = std::min(static_cast<std::size_t>(t), n - 2);
    double f = t - i;
    return grid[i] * (1.0 - f) + grid[i + 1] * f;
  }

  std::size_t nbins_;
  double fixed_lengthscale_;
  bool matern32_;
  double lam_min_req_, lam_max_req_;

  bool active_;
  std::string message_;
  double lam_min_, lam_max_, lambda_ref_, g_ref_, gmean_;
  double ell_, sigf_, noise_inflation_, mean_noise_, lml_;
  std::size_t n_samples_;
  bool ell_at_lower_, ell_at_upper_, noise_at_upper_;
  double grid_step_;
  std::vector<double> grid_g_, grid_sd_;
  std::vector<double> train_lam_, train_h_, train_sd_;
};

}  // namespace dials_scaling

#endif  // DIALS_ALGORITHMS_SCALING_TOF_LAMBDA_NORM_H
