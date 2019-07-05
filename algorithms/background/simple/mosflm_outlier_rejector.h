/*
 * mosflm_outlier_rejector.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_BACKGROUND_MOSFLM_OUTLIER_REJECTOR_H
#define DIALS_ALGORITHMS_BACKGROUND_MOSFLM_OUTLIER_REJECTOR_H

#include <algorithm>
#include <scitbx/array_family/ref_reductions.h>
#include <boost/math/special_functions/erf.hpp>
#include <scitbx/math/mean_and_variance.h>
#include <scitbx/matrix/inversion.h>
#include <dials/array_family/sort_index.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace background {

  using model::Background;
  using model::BackgroundUsed;
  using model::Overlapped;
  using model::Valid;
  using scitbx::matrix::inversion_in_place;

  /**
   * Compute the background plane from a subset and select pixels within a few
   * standard deviations of the subset plane for use in calculating the final
   * background.
   */
  class MosflmOutlierRejector : public OutlierRejector {
  public:
    /**
     * @param fraction The fraction to use in initial estimate
     * @param nsigma The threshold for outliers
     */
    MosflmOutlierRejector(double fraction, double n_sigma)
        : fraction_(fraction), n_sigma_(n_sigma) {
      DIALS_ASSERT(fraction_ > 0 && fraction_ <= 1.0);
      DIALS_ASSERT(n_sigma_ > 0);
    }

    /**
     * @params shoebox The shoebox profile
     * @params mask The shoebox mask
     */
    virtual void mark(const af::const_ref<double, af::c_grid<3> > &data,
                      af::ref<int, af::c_grid<3> > mask) const {
      // Check the input
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));

      // Mark for all slices
      af::c_grid<2> accessor(data.accessor()[1], data.accessor()[2]);
      std::size_t xysize = accessor[0] * accessor[1];
      for (std::size_t i = 0; i < data.accessor()[0]; ++i) {
        // Get the 2D slices
        af::const_ref<double, af::c_grid<2> > data_2d(&data[i * xysize], accessor);
        af::ref<int, af::c_grid<2> > mask_2d(&mask[i * xysize], accessor);

        // Compute the initial mask using a subset of the available pixels
        compute_initial_mask(data_2d, mask_2d);

        // Compute the background plane using the subset of pixels
        double a = 0;
        double b = 0;
        double c = 0;
        compute_background(data_2d, mask_2d, a, b, c);

        // Apply a correction to the background based on the fraction of pixels
        // chosen for use in the background plane calculation.
        compute_correction(a, b, c);

        // Check all pixels for outliers against the initial plane
        compute_final_mask(data_2d, mask_2d, a, b, c);
      }
    }

  private:
    /**
     * Compare index array by pixel value
     */
    struct compare_pixel_value {
      af::const_ref<double> data_;
      compare_pixel_value(const af::const_ref<double> &data) : data_(data) {}
      bool operator()(std::size_t a, std::size_t b) {
        return data_[a] < data_[b];
      }
    };

    /**
     * Calculate the initial mask. Select the fraction of pixels with the lowest
     * intensity and then update the mask for those pixels.
     */
    void compute_initial_mask(const af::const_ref<double, af::c_grid<2> > &data,
                              af::ref<int, af::c_grid<2> > mask) const {
      int code = Valid | Background;
      std::vector<std::size_t> index;
      index.reserve(data.size());
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & code) == code && (mask[i] & Overlapped) == 0) {
          index.push_back(i);
        } else {
          mask[i] &= ~BackgroundUsed;
        }
      }
      DIALS_ASSERT(index.size() > 0);
      std::sort(index.begin(), index.end(), compare_pixel_value(data.as_1d()));
      std::size_t nactive = (std::size_t)std::floor(fraction_ * index.size() + 0.5);
      DIALS_ASSERT(nactive > 0 && nactive <= index.size());
      for (std::size_t i = 0; i < nactive; ++i) {
        mask[index[i]] |= BackgroundUsed;
      }
    }

    /**
     * Compute the background plane given shoebox data and mask
     */
    void compute_background(const af::const_ref<double, af::c_grid<2> > &data,
                            const af::const_ref<int, af::c_grid<2> > &mask,
                            double &a,
                            double &b,
                            double &c) const {
      std::vector<double> A(3 * 3, 0);
      std::vector<double> B(3, 0);
      int hy = data.accessor()[0] / 2;
      int hx = data.accessor()[1] / 2;
      for (std::size_t j = 0; j < mask.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < mask.accessor()[1]; ++i) {
          if (mask(j, i) & BackgroundUsed) {
            double x = ((int)i - hx);
            double y = ((int)j - hy);
            double p = data(j, i);
            A[0] += 1;
            A[1] += x;
            A[2] += y;
            A[4] += x * x;
            A[5] += x * y;
            A[8] += y * y;
            B[0] += p;
            B[1] += x * p;
            B[2] += y * p;
          }
        }
      }
      A[3] = A[1];
      A[6] = A[2];
      A[7] = A[5];
      inversion_in_place(&A[0], 3, &B[0], 1);
      c = B[0];
      a = B[1];
      b = B[2];
    }

    /**
     * Compute the correction to the background plane for the initial estimate
     * where only a fraction of the input pixels are given.
     */
    void compute_correction(double &a, double &b, double &c) const {
      double factor = 0.8;
      if (fraction_ >= 0.999) {
        factor = 0.0;
      } else if (fraction_ >= 0.9) {
        factor = 0.05;
      } else if (fraction_ >= 0.7) {
        factor = 0.2;
      } else if (fraction_ >= 0.5) {
        factor = 0.4;
      } else if (fraction_ >= 0.3) {
        factor = 0.6;
      }
      c = c + factor * std::sqrt(std::abs(c));
    }

    /**
     * Compute the final mask. Check the data for outliers by looking at the
     * deviation of each point from the plane.
     */
    void compute_final_mask(const af::const_ref<double, af::c_grid<2> > &data,
                            af::ref<int, af::c_grid<2> > mask,
                            double a,
                            double b,
                            double c) const {
      std::size_t noutlier = 0;
      std::size_t nbackground = 0;
      double max_background = n_sigma_ * std::sqrt(std::abs(c));
      int code = Valid | Background;
      int hy = data.accessor()[0] / 2;
      int hx = data.accessor()[1] / 2;
      for (std::size_t j = 0; j < mask.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < mask.accessor()[1]; ++i) {
          if ((mask(j, i) & code) == code && (mask(j, i) & Overlapped) == 0) {
            double x = ((int)i - hx);
            double y = ((int)j - hy);
            double p1 = data(j, i);
            double p2 = a * x + b * y + c;
            double d = std::abs(p1 - p2);
            if (d > max_background) {
              mask(j, i) &= ~BackgroundUsed;
              noutlier++;
            } else {
              mask(j, i) |= BackgroundUsed;
              nbackground++;
            }
          }
        }
      }
    }

    double fraction_, n_sigma_;
  };

}}}  // namespace dials::algorithms::background

#endif  // DIALS_ALGORITHMS_BACKGROUND_MOSFLM_OUTLIER_REJECTOR_H
