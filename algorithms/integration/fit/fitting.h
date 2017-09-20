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

  using scitbx::matrix::inversion_in_place;
  using scitbx::af::sum;
  using scitbx::vec2;

  /* /** */
  /*  * Class to fit the observed with the reference profile */
  /*  *1/ */
  /* template <typename FloatType = double> */
  /* class ProfileFitting { */
  /* public: */

  /*   typedef FloatType float_type; */

  /*   /** */
  /*    * Instantiate the fitting algorithm with the reflection profile */
  /*    * @param p The profile to fit to */
  /*    * @param c The contents of the pixels */
  /*    * @param b The background of the pixels */
  /*    *1/ */
  /*   ProfileFitting(const af::const_ref<FloatType, af::c_grid<3> > &p, */
  /*                  const af::const_ref<bool, af::c_grid<3> > &m, */
  /*                  const af::const_ref<FloatType, af::c_grid<3> > &c, */
  /*                  const af::const_ref<FloatType, af::c_grid<3> > &b, */
  /*                  double eps = 1e-3, */
  /*                  std::size_t max_iter = 10) */
  /*   { */
  /*     // Check the input */
  /*     DIALS_ASSERT(p.size() == m.size()); */
  /*     DIALS_ASSERT(p.size() == c.size()); */
  /*     DIALS_ASSERT(p.size() == b.size()); */
  /*     DIALS_ASSERT(eps > 0.0); */
  /*     DIALS_ASSERT(max_iter >= 1); */

  /*     double sumc = 0; */
  /*     double sumb = 0; */
  /*     double sums = 0; */
  /*     double minc = -1; */
  /*     double minI = 9999; */
  /*     for (std::size_t i = 0; i < m.size(); ++i) { */
  /*       if (m[i]) { */
  /*         DIALS_ASSERT(p[i] >= 0); */
  /*         DIALS_ASSERT(b[i] >= 0); */
  /*         DIALS_ASSERT(c[i] >= 0); */
  /*         sumc += c[i]; */
  /*         sumb += b[i]; */
  /*         sums += p[i]; */
  /*         if (minc < 0 || c[i] < minc) minc = c[i]; */
  /*         if (p[i] > 0 && b[i] / p[i] < minI) minI = b[i] / p[i]; */
  /*       } */
  /*     } */
  /*     minI = -minI; */
  /*     DIALS_ASSERT(sumb >= 0); */
  /*     DIALS_ASSERT(sumc >= 0); */
  /*     DIALS_ASSERT(sums > 0); */
  /*     DIALS_ASSERT(minI <= 0); */

  /*     // Iterate to calculate the intensity. Exit if intensity goes less */
  /*     // than zero or if the tolerance or number of iteration is reached. */
  /*     double I0 = sumc - sumb; */
  /*     if (I0 < minI) { */
  /*       I0 = minI+1e-3; */
  /*     } */
  /*     vec2<double> I(0.0, 0.0); */
  /*     for (niter_ = 0; niter_ < max_iter; ++niter_) { */
  /*       double sum1 = 0.0; */
  /*       double sum2 = 0.0; */
  /*       double sumv = 0.0; */
  /*       for (std::size_t i = 0; i < p.size(); ++i) { */
  /*         if (m[i] && p[i] > 0) { */
  /*           double v = b[i] + I0 * p[i]; */
  /*           DIALS_ASSERT(v > 0); */
  /*           sumv += v; */
  /*           if (v > 0) { */
  /*             sum1 += (c[i] - b[i]) * p[i] / v; */
  /*             sum2 += p[i] * p[i] / v; */
  /*           } */
  /*         } */
  /*       } */
  /*       DIALS_ASSERT(sum2 > 0); */
  /*       I[0] = sum1 / sum2; */
  /*       I[1] = sumv; */
  /*       if ((error_ = std::abs(I[0] - I0)) < eps) { */
  /*         break; */
  /*       } */
  /*       if (I[0] < minI+1e-5) { */
  /*         intensity_ = (sumc - sumb) / sums; */
  /*         variance_ = std::abs(intensity_ * sums) + sumb; */
  /*         correlation_ = 0; */
  /*         break; */
  /*       } */
  /*       I0 = I[0]; */
  /*     } */
  /*     DIALS_ASSERT(I[1] >= 0); */

  /*     // Set the intensity and variance */
  /*     intensity_ = I[0]; */
  /*     variance_ = I[1]; */
  /*     correlation_ = compute_correlation(p, m, c, b); */
  /*     //rmsd_ = compute_rmsd(I[0], p, m, c, b); */
  /*   } */

  /*   /** */
  /*    * @returns The intensity */
  /*    *1/ */
  /*   double intensity() const { */
  /*     return intensity_; */
  /*   } */

  /*   /** */
  /*    * @returns the variance */
  /*    *1/ */
  /*   double variance() const { */
  /*     return variance_; */
  /*   } */

  /*   /** */
  /*    * @returns the correlation */
  /*    *1/ */
  /*   double correlation() const { */
  /*     return correlation_; */
  /*   } */

  /*   /** */
  /*    * @returns The number of iterations */
  /*    *1/ */
  /*   std::size_t niter() const { */
  /*     return niter_; */
  /*   } */

  /*   /** */
  /*    * @returns The error in the fit */
  /*    *1/ */
  /*   double error() const { */
  /*     return error_; */
  /*   } */

  /*   /** */
  /*    * @returns The rmsd in the fit */
  /*    *1/ */
  /*   /1* double rmsd() const { *1/ */
  /*   /1*   return rmsd_; *1/ */
  /*   /1* } *1/ */

  /* private: */

  /*   /** */
  /*    * Evaluate the next intensity iteration. */
  /*    * @ returns The estimate of the intensity */
  /*    *1/ */
  /*   vec2<double> */
  /*   estimate_intensity(const af::const_ref<FloatType, af::c_grid<3> > &p, */
  /*                      const af::const_ref<bool, af::c_grid<3> > &m, */
  /*                      const af::const_ref<FloatType, af::c_grid<3> > &c, */
  /*                      const af::const_ref<FloatType, af::c_grid<3> > &b, */
  /*                      double I) const { */
  /*     double sum1 = 0.0; */
  /*     double sum2 = 0.0; */
  /*     double sumv = 0.0; */
  /*     for (std::size_t i = 0; i < p.size(); ++i) { */
  /*       if (m[i]) { */
  /*         double v = std::abs(b[i]) + std::abs(I * p[i]); */
  /*         sumv += v; */
  /*         if (v > 0) { */
  /*           sum1 += (c[i] - b[i]) * p[i] / v; */
  /*           sum2 += p[i] * p[i] / v; */
  /*         } */
  /*       } */
  /*     } */
  /*     return vec2<double>(sum2 != 0 ? sum1 / sum2 : 0.0, sumv); */
  /*   } */

  /*   /** */
  /*    * Compute the correlation coefficient between the profile and reference */
  /*    *1/ */
  /*   double */
  /*   compute_correlation(const af::const_ref<FloatType, af::c_grid<3> > &p, */
  /*                       const af::const_ref<bool, af::c_grid<3> > &m, */
  /*                       const af::const_ref<FloatType, af::c_grid<3> > &c, */
  /*                       const af::const_ref<FloatType, af::c_grid<3> > &b) const { */
  /*     double xb = 0.0, yb = 0.0; */
  /*     std::size_t count = 0; */
  /*     for (std::size_t i = 0; i < p.size(); ++i) { */
  /*       if (m[i]) { */
  /*         xb += intensity_*p[i] + b[i]; */
  /*         yb += c[i]; */
  /*         count++; */
  /*       } */
  /*     } */
  /*     DIALS_ASSERT(count > 0); */
  /*     xb /= count; */
  /*     yb /= count; */
  /*     double sdxdy = 0.0, sdx2 = 0.0, sdy2 = 0.0; */
  /*     for (std::size_t i = 0; i < p.size(); ++i) { */
  /*       if (m[i]) { */
  /*         double dx = (intensity_*p[i] + b[i]) - xb; */
  /*         double dy = c[i] - yb; */
  /*         sdxdy += dx*dy; */
  /*         sdx2 += dx*dx; */
  /*         sdy2 += dy*dy; */
  /*       } */
  /*     } */
  /*     double result = 0.0; */
  /*     if (sdx2 > 0.0 && sdy2 > 0.0) { */
  /*       result = sdxdy / (std::sqrt(sdx2) * std::sqrt(sdy2)); */
  /*     } */
  /*     return result; */
  /*   } */

  /*   /** */
  /*    * Compute the rmsd between the profile and reference */
  /*    *1/ */
  /*   /1* double *1/ */
  /*   /1* compute_rmsd(double I, *1/ */
  /*   /1*              const af::const_ref<FloatType, af::c_grid<3> > &p, *1/ */
  /*   /1*              const af::const_ref<bool, af::c_grid<3> > &m, *1/ */
  /*   /1*              const af::const_ref<FloatType, af::c_grid<3> > &c, *1/ */
  /*   /1*              const af::const_ref<FloatType, af::c_grid<3> > &b) const { *1/ */
  /*   /1*   double sum = 0.0; *1/ */
  /*   /1*   double ymax = 0.0; *1/ */
  /*   /1*   double ymin = 0.0; *1/ */
  /*   /1*   std::size_t count = 0; *1/ */
  /*   /1*   for (std::size_t i = 0; i < p.size(); ++i) { *1/ */
  /*   /1*     if (m[i]) { *1/ */
  /*   /1*       double y = (c[i] - b[i]); *1/ */
  /*   /1*       if (count == 0) { *1/ */
  /*   /1*         ymax = y; *1/ */
  /*   /1*         ymin = y; *1/ */
  /*   /1*       } else { *1/ */
  /*   /1*         if (ymax < y) ymax = y; *1/ */
  /*   /1*         if (ymin > y) ymin = y; *1/ */
  /*   /1*       } *1/ */
  /*   /1*       double x = I * p[i] - y; *1/ */
  /*   /1*       sum += x * x; *1/ */
  /*   /1*       count++; *1/ */
  /*   /1*     } *1/ */
  /*   /1*   } *1/ */
  /*   /1*   DIALS_ASSERT(count > 0); *1/ */
  /*   /1*   DIALS_ASSERT(ymax > ymin); *1/ */
  /*   /1*   sum /= count; *1/ */
  /*   /1*   return std::sqrt(sum) / (ymax - ymin); *1/ */
  /*   /1* } *1/ */

  /*   double intensity_; */
  /*   double variance_; */
  /*   double correlation_; */
  /*   std::size_t niter_; */
  /*   double error_; */
  /*   /1* double rmsd_; *1/ */
  /* }; */

  namespace detail {

    template <typename T, std::size_t N>
    af::const_ref<T, af::c_grid<2> > as_2d(const af::const_ref< T, af::c_grid<N> > &src) {
      DIALS_ASSERT(N > 2);
      std::size_t size = 1;
      for (std::size_t i = 1; i < src.accessor().size(); ++i) {
        size *= src.accessor()[i];
      }
      af::c_grid<2> accessor(src.accessor()[0], size);
      DIALS_ASSERT(accessor[1] * accessor[0] == src.size());
      return af::const_ref<T, af::c_grid<2> >(&src[0], accessor);
    }

  }

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
    ProfileFitter(
        const af::const_ref<T> &d,
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
    ProfileFitter(
        const af::const_ref<T> &d,
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
    ProfileFitter(
        const af::const_ref<T, af::c_grid<2> > &d,
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
    ProfileFitter(
        const af::const_ref<T, af::c_grid<3> > &d,
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
    ProfileFitter(
        const af::const_ref<T, af::c_grid<2> > &d,
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
    ProfileFitter(
        const af::const_ref<T, af::c_grid<3> > &d,
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
    af::small<double,10> intensity() const {
      return intensity_;
    }

    /**
     * @returns the variance
     */
    af::small<double, 10> variance() const {
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
    void fit(
        const af::const_ref<T> &d,
        const af::const_ref<T> &b,
        const af::const_ref<bool> &m,
        const af::const_ref<T> &p,
        double eps,
        std::size_t maxiter) {

      const double TINY = 1e-7;
      const double MASSIVE = 1e10;

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
      double mind = -1;
      double minI = MASSIVE;
      for (std::size_t i = 0; i < m.size(); ++i) {
        if (m[i]) {
          DIALS_ASSERT(p[i] >= 0);
          DIALS_ASSERT(b[i] >= 0);
          DIALS_ASSERT(d[i] >= 0);
          sumd += d[i];
          sumb += b[i];
          sump += p[i];
          if (mind < 0 || d[i] < mind) mind = d[i];
          if (p[i] > 0 && b[i] / p[i] < minI) minI = b[i] / p[i];
        }
      }
      minI = -minI;
      DIALS_ASSERT(sumb >= 0);
      DIALS_ASSERT(sumd >= 0);
      DIALS_ASSERT(sump > 0);
      DIALS_ASSERT(minI <= 0);

      // Iterate to calculate the intensity. Exit if intensity goes less
      // than zero or if the tolerance or number of iteration is reached.
      double I0 = sumd - sumb;
      if (I0 < minI + TINY) {
        I0 = minI + TINY;
      }
      double I = 0.0;
      double V = 0.0;
      for (niter_ = 0; niter_ < maxiter; ++niter_) {
        double sum1 = 0.0;
        double sum2 = 0.0;
        double sumv = 0.0;
        for (std::size_t i = 0; i < p.size(); ++i) {
          if (m[i] && p[i] > 0) {
            double v = b[i] + I0 * p[i];
            DIALS_ASSERT(v > 0);
            sumv += v;
            if (v > 0) {
              sum1 += (d[i] - b[i]) * p[i] / v;
              sum2 += p[i] * p[i] / v;
            }
          }
        }
        DIALS_ASSERT(sum2 > 0);
        I = sum1 / sum2;
        V = sumv;
        if ((error_ = std::abs(I - I0)) < eps) {
          break;
        }
        if (I < minI + TINY) {
          I = (sumd - sumb) / sump;
          V = std::abs(I * sump) + sumb;
          correlation_ = 0;
          break;
        }
        I0 = I;
      }
      DIALS_ASSERT(V >= 0);

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
     * @param d The data array
     * @param b The background array
     * @param m The mask array
     * @param p The profile array
     * @param eps The tolerance
     * @param maxiter The maximum number of iterations
     */
    void fit(
        const af::const_ref<T> &d,
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

        // Compute the variance for the given estimate
        std::copy(b.begin(), b.end(), v.begin());
        for (std::size_t j = 0; j < M; ++j) {
          for (std::size_t i = 0; i < N; ++i) {
            v[i] += I[j] * p(j,i);
          }
        }

        // Compute the matrices to do the profile fitting
        std::fill(I.begin(), I.end(), 0);
        for (std::size_t k = 0; k < M; ++k) {
          for (std::size_t i = 0; i < N; ++i) {
            if (m[i]) {
              DIALS_ASSERT(v[i] > 0);
              I[k] += p(k,i) * (d[i] - b[i]) / v[i];
            }
          }
        }

        std::fill(A.begin(), A.end(), 0);
        for (std::size_t k = 0; k < M; ++k) {
          for (std::size_t i = 0; i < N; ++i) {
            if (m[i]) {
              A[k+k*M] += p(k,i)*p(k,i) / v[i];
            }
          }
          for (std::size_t j = k+1; j < M; ++j) {
            for (std::size_t i = 0; i < N; ++i) {
              if (m[i]) {
                A[j+k*M] += p(k,i)*p(j,i) / v[i];
                A[k+j*M] = A[j+k*M];
              }
            }
          }
        }

        // Do the inversion to compute the profile fits
        inversion_in_place(&A[0], M, &I[0], 1);

        // Compute error
        error_ = 0;
        for (std::size_t j = 0; j < M; ++j) {
          DIALS_ASSERT(I[j] > 0);
          error_ += (I[j] - I0[j])*(I[j] - I0[j]);
        }
        if (error_ < eps*eps) {
          break;
        }

        // Set the old intensity
        I0 = I;
      }

      // Set the return values
      for (std::size_t j = 0; j < M; ++j) {

        double V = 0;
        for (std::size_t i = 0; i < N; ++i) {
          if (m[i] && p(j,i) > 0) {
            V += b[i] + I[j] * p(j,i);
          }
        }

        intensity_.push_back(I[j]);
        variance_.push_back(V);
      }

      // Compute the correlation
      correlation_ = compute_correlation(d, b, m, p);

    }

    /**
     * Compute the correlation for a single reflection
     */
    double compute_correlation(
        const af::const_ref<T> &d,
        const af::const_ref<T> &b,
        const af::const_ref<bool> &m,
        const af::const_ref<T> &p) const {

      // Compute the mean observed and predicted
      double xb = 0.0, yb = 0.0;
      std::size_t count = 0;
      for (std::size_t i = 0; i < p.size(); ++i) {
        if (m[i]) {
          xb += intensity_[0]*p[i] + b[i];
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
          double dx = (intensity_[0]*p[i] + b[i]) - xb;
          double dy = d[i] - yb;
          sdxdy += dx*dy;
          sdx2 += dx*dx;
          sdy2 += dy*dy;
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
    double compute_correlation(
        const af::const_ref<T> &d,
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
            xb += intensity_[j]*p(j,i);
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
            dx += intensity_[j]*p(j,i);
          }
          sdxdy += dx*dy;
          sdx2 += dx*dx;
          sdy2 += dy*dy;
        }
      }

      // Compute the correlation
      double result = 0.0;
      if (sdx2 > 0.0 && sdy2 > 0.0) {
        result = sdxdy / (std::sqrt(sdx2) * std::sqrt(sdy2));
      }
      return result;
    }

    af::small<double,10> intensity_;
    af::small<double,10> variance_;
    double correlation_;
    std::size_t niter_;
    std::size_t maxiter_;
    double error_;
  };


}}

#endif /* DIALS_ALGORITHMS_INTEGRATION_FIT_FITTING_H */
