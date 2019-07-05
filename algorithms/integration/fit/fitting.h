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
#include <vector>
#include <scitbx/vec2.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <dials/model/data/mask_code.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <scitbx/matrix/inversion.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::vec2;
  using scitbx::af::sum;
  using scitbx::matrix::inversion_in_place;

  namespace detail {

    /**
     * Return multi dimensional array as 2D (i.e. flatten dimensons 1-N)
     */
    template <typename T, std::size_t N>
    af::const_ref<T, af::c_grid<2> > as_2d(
      const af::const_ref<T, af::c_grid<N> > &src) {
      DIALS_ASSERT(N > 2);
      std::size_t size = 1;
      for (std::size_t i = 1; i < src.accessor().size(); ++i) {
        size *= src.accessor()[i];
      }
      af::c_grid<2> accessor(src.accessor()[0], size);
      DIALS_ASSERT(accessor[1] * accessor[0] == src.size());
      return af::const_ref<T, af::c_grid<2> >(&src[0], accessor);
    }

  }  // namespace detail

  /**
   * A class to do profile fitting.
   *
   * The class is able to do profile fitting of a single reflection or
   * deconvolution of multiple overlapping reflections.
   */
  template <typename T = double>
  class ProfileFitter {
  public:
    typedef T float_type;

    /**
     * Profile fit a single reflection
     */
    ProfileFitter(const af::const_ref<T> &d,
                  const af::const_ref<T> &b,
                  const af::const_ref<bool> &m,
                  const af::const_ref<T> &p,
                  double eps = 1e-3,
                  std::size_t maxiter = 10) {
      fit(d, b, m, p, eps, maxiter);
    }

    /**
     * Profile fit and deconvolve multiple reflections
     */
    ProfileFitter(const af::const_ref<T> &d,
                  const af::const_ref<T> &b,
                  const af::const_ref<bool> &m,
                  const af::const_ref<T, af::c_grid<2> > &p,
                  double eps = 1e-3,
                  std::size_t maxiter = 10) {
      fit(d, b, m, p, eps, maxiter);
    }

    /**
     * Profile fit a single 2D reflection
     */
    ProfileFitter(const af::const_ref<T, af::c_grid<2> > &d,
                  const af::const_ref<T, af::c_grid<2> > &b,
                  const af::const_ref<bool, af::c_grid<2> > &m,
                  const af::const_ref<T, af::c_grid<2> > &p,
                  double eps = 1e-3,
                  std::size_t maxiter = 10) {
      fit(d.as_1d(), b.as_1d(), m.as_1d(), p.as_1d(), eps, maxiter);
    }

    /**
     * Profile fit and deconvolve multiple 2D reflections
     */
    ProfileFitter(const af::const_ref<T, af::c_grid<3> > &d,
                  const af::const_ref<T, af::c_grid<3> > &b,
                  const af::const_ref<bool, af::c_grid<3> > &m,
                  const af::const_ref<T, af::c_grid<3> > &p,
                  double eps = 1e-3,
                  std::size_t maxiter = 10) {
      fit(d.as_1d(), b.as_1d(), m.as_1d(), p.as_1d(), eps, maxiter);
    }

    /**
     * Profile fit a single 3D reflection
     */
    ProfileFitter(const af::const_ref<T, af::c_grid<2> > &d,
                  const af::const_ref<T, af::c_grid<2> > &b,
                  const af::const_ref<bool, af::c_grid<2> > &m,
                  const af::const_ref<T, af::c_grid<3> > &p,
                  double eps = 1e-3,
                  std::size_t maxiter = 10) {
      fit(d.as_1d(), b.as_1d(), m.as_1d(), detail::as_2d(p), eps, maxiter);
    }

    /**
     * Profile fit and deconvolve multiple 3D reflections
     */
    ProfileFitter(const af::const_ref<T, af::c_grid<3> > &d,
                  const af::const_ref<T, af::c_grid<3> > &b,
                  const af::const_ref<bool, af::c_grid<3> > &m,
                  const af::const_ref<T, af::c_grid<4> > &p,
                  double eps = 1e-3,
                  std::size_t maxiter = 10) {
      fit(d.as_1d(), b.as_1d(), m.as_1d(), detail::as_2d(p), eps, maxiter);
    }

    /**
     * @returns The intensity
     */
    af::small<double, 10> intensity() const {
      DIALS_ASSERT(intensity_.size() > 0);
      return intensity_;
    }

    /**
     * @returns the variance
     */
    af::small<double, 10> variance() const {
      DIALS_ASSERT(variance_.size() > 0);
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
     * @returns The maximum number of iterations
     */
    std::size_t maxiter() const {
      return maxiter_;
    }

    /**
     * @returns The error in the fit
     */
    double error() const {
      return error_;
    }

  protected:
    /**
     * Profile fit a single reflection
     *
     * @param d The data array
     * @param b The background array
     * @param m The mask array
     * @param p The profile array
     * @param eps The tolerance
     * @param maxiter The maximum number of iterations
     */
    void fit(const af::const_ref<T> &d,
             const af::const_ref<T> &b,
             const af::const_ref<bool> &m,
             const af::const_ref<T> &p,
             double eps,
             std::size_t maxiter) {
      // Save the max iter
      maxiter_ = maxiter;

      // Check the input
      DIALS_ASSERT(d.size() == b.size());
      DIALS_ASSERT(d.size() == m.size());
      DIALS_ASSERT(d.size() == p.size());
      DIALS_ASSERT(eps > 0.0);
      DIALS_ASSERT(maxiter >= 1);

      // Compute the sums of the background and foreground
      double sumd = 0;
      double sumb = 0;
      double sump = 0;
      for (std::size_t i = 0; i < m.size(); ++i) {
        if (m[i]) {
          DIALS_ASSERT(p[i] >= 0);
          sumd += d[i];
          sumb += b[i];
          sump += p[i];
        }
      }
      DIALS_ASSERT(sumb >= 0);
      DIALS_ASSERT(sumd >= 0);
      DIALS_ASSERT(sump > 0);

      // Iterate to calculate the intensity. Exit if intensity goes less
      // than zero or if the tolerance or number of iteration is reached.
      double I0 = sumd - sumb;
      double I = 0.0;
      double V = 0.0;
      for (niter_ = 0; niter_ < maxiter; ++niter_) {
        double sum1 = 0.0;
        double sum2 = 0.0;
        for (std::size_t i = 0; i < p.size(); ++i) {
          if (m[i] && p[i] > 0) {
            double v = 1e-10 + std::abs(b[i]) + std::abs(I0 * p[i]);
            DIALS_ASSERT(v > 0);
            sum1 += (d[i] - b[i]) * p[i] / v;
            sum2 += p[i] * p[i] / v;
          }
        }
        DIALS_ASSERT(sum2 > 0);
        I = sum1 / sum2;
        V = std::abs(I) + std::abs(sumb);
        if ((error_ = std::abs(I - I0)) < eps) {
          break;
        }
        I0 = I;
      }

      // If niter is too large replace with summation results
      if (niter_ >= maxiter) {
        I = sumd - sumb;
        V = std::abs(I) + std::abs(sumb);
      }
      DIALS_ASSERT(V >= 0);
      DIALS_ASSERT(V >= I);

      // Set the intensity and variance
      intensity_.push_back(I);
      variance_.push_back(V);

      // Compute the correlation
      correlation_ = compute_correlation(d, b, m, p);
    }

    /**
     * Profile fit and deconvolve multiple reflections
     *
     * The profile array must be a 2D array with dimensions (M, N) where M is
     * the number of profiles and N is the size of the data array. Each slice of
     * the profile array is then a profile for a separate reflection which must
     * be zero outside the reflection foreground.
     *
     * FIXME The variance structure for I < 0 needs fixing
     *
     * @param d The data array
     * @param b The background array
     * @param m The mask array
     * @param p The profile array
     * @param eps The tolerance
     * @param maxiter The maximum number of iterations
     */
    void fit(const af::const_ref<T> &d,
             const af::const_ref<T> &b,
             const af::const_ref<bool> &m,
             const af::const_ref<T, af::c_grid<2> > &p,
             double eps,
             std::size_t maxiter) {
      // Save the max iter
      maxiter_ = maxiter;

      // Check some input
      DIALS_ASSERT(d.size() == b.size());
      DIALS_ASSERT(d.size() == m.size());
      DIALS_ASSERT(p.accessor()[0] > 0);
      DIALS_ASSERT(p.accessor()[0] < intensity_.max_size());
      DIALS_ASSERT(p.accessor()[1] == d.size());
      DIALS_ASSERT(eps > 0.0);
      DIALS_ASSERT(maxiter >= 1);

      // The dimensions of the problem
      std::size_t N = d.size();
      std::size_t M = p.accessor()[0];

      // Allocate some arrays
      std::vector<double> v(N);
      std::vector<double> I(M, 1);
      std::vector<double> I0(M, 1);
      std::vector<double> A(M * M);

      // Iterate a number of times
      for (niter_ = 0; niter_ < maxiter; ++niter_) {
        //
        // Compute the variance for the given estimate
        std::copy(b.begin(), b.end(), v.begin());
        for (std::size_t j = 0; j < M; ++j) {
          for (std::size_t i = 0; i < N; ++i) {
            v[i] += std::max(1.0 / N, std::abs(I[j])) * p(j, i);
          }
        }

        // Compute the matrices to do the profile fitting
        std::fill(I.begin(), I.end(), 0);
        for (std::size_t k = 0; k < M; ++k) {
          for (std::size_t i = 0; i < N; ++i) {
            if (m[i]) {
              DIALS_ASSERT(v[i] > 0);
              I[k] += p(k, i) * (d[i] - b[i]) / v[i];
            }
          }
        }

        std::fill(A.begin(), A.end(), 0);
        for (std::size_t k = 0; k < M; ++k) {
          for (std::size_t i = 0; i < N; ++i) {
            if (m[i]) {
              A[k + k * M] += p(k, i) * p(k, i) / v[i];
            }
          }
          for (std::size_t j = k + 1; j < M; ++j) {
            for (std::size_t i = 0; i < N; ++i) {
              if (m[i]) {
                A[j + k * M] += p(k, i) * p(j, i) / v[i];
                A[k + j * M] = A[j + k * M];
              }
            }
          }
        }

        // Do the inversion to compute the profile fits
        inversion_in_place(&A[0], M, &I[0], 1);

        // Compute error
        error_ = 0;
        for (std::size_t j = 0; j < M; ++j) {
          // DIALS_ASSERT(I[j] > 0);
          error_ += (I[j] - I0[j]) * (I[j] - I0[j]);
        }
        if (error_ < eps * eps) {
          break;
        }

        // Set the old intensity
        I0 = I;
      }

      // Set the return values
      for (std::size_t j = 0; j < M; ++j) {
        double V = std::abs(I[j]);
        for (std::size_t i = 0; i < N; ++i) {
          if (m[i]) {
            V += b[i];
          }
        }

        DIALS_ASSERT(V >= 0);
        DIALS_ASSERT(V >= I[j]);

        intensity_.push_back(I[j]);
        variance_.push_back(V);
      }

      // Compute the correlation
      correlation_ = compute_correlation(d, b, m, p);
    }

    /**
     * Compute the correlation for a single reflection
     */
    double compute_correlation(const af::const_ref<T> &d,
                               const af::const_ref<T> &b,
                               const af::const_ref<bool> &m,
                               const af::const_ref<T> &p) const {
      // Compute the mean observed and predicted
      double xb = 0.0, yb = 0.0;
      std::size_t count = 0;
      for (std::size_t i = 0; i < p.size(); ++i) {
        if (m[i]) {
          xb += intensity_[0] * p[i] + b[i];
          yb += d[i];
          count++;
        }
      }
      DIALS_ASSERT(count > 0);
      xb /= count;
      yb /= count;

      // Compute the variance
      double sdxdy = 0.0, sdx2 = 0.0, sdy2 = 0.0;
      for (std::size_t i = 0; i < p.size(); ++i) {
        if (m[i]) {
          double dx = (intensity_[0] * p[i] + b[i]) - xb;
          double dy = d[i] - yb;
          sdxdy += dx * dy;
          sdx2 += dx * dx;
          sdy2 += dy * dy;
        }
      }

      // Compute the correlation
      double result = 0.0;
      if (sdx2 > 0.0 && sdy2 > 0.0) {
        result = sdxdy / (std::sqrt(sdx2) * std::sqrt(sdy2));
      }
      return result;
    }

    /**
     * Compute the correlation for multiple profiles
     */
    double compute_correlation(const af::const_ref<T> &d,
                               const af::const_ref<T> &b,
                               const af::const_ref<bool> &m,
                               const af::const_ref<T, af::c_grid<2> > &p) const {
      // Compute the mean observed and predicted
      double xb = 0.0, yb = 0.0;
      std::size_t count = 0;
      for (std::size_t i = 0; i < p.accessor()[1]; ++i) {
        if (m[i]) {
          count++;
          yb += d[i];
          xb += b[i];
          for (std::size_t j = 0; j < p.accessor()[0]; ++j) {
            xb += intensity_[j] * p(j, i);
          }
        }
      }
      DIALS_ASSERT(count > 0);
      xb /= count;
      yb /= count;

      // Compute the variance
      double sdxdy = 0.0, sdx2 = 0.0, sdy2 = 0.0;
      for (std::size_t i = 0; i < p.accessor()[1]; ++i) {
        if (m[i]) {
          double dy = d[i] - yb;
          double dx = b[i] - xb;
          for (std::size_t j = 0; j < p.accessor()[0]; ++j) {
            dx += intensity_[j] * p(j, i);
          }
          sdxdy += dx * dy;
          sdx2 += dx * dx;
          sdy2 += dy * dy;
        }
      }

      // Compute the correlation
      double result = 0.0;
      if (sdx2 > 0.0 && sdy2 > 0.0) {
        result = sdxdy / (std::sqrt(sdx2) * std::sqrt(sdy2));
      }
      return result;
    }

    af::small<double, 10> intensity_;
    af::small<double, 10> variance_;
    double correlation_;
    std::size_t niter_;
    std::size_t maxiter_;
    double error_;
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_FIT_FITTING_H */
