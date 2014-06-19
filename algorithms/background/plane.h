/*
 * plane.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *  Author: Luis Fuentes-Montero
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_PLANE_H
#define DIALS_ALGORITHMS_BACKGROUND_PLANE_H

#include <vector>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <dials/error.h>
#include <scitbx/matrix/inversion.h>

namespace dials { namespace algorithms { namespace background {

  using dials::model::Valid;
  using dials::model::Background;
  using dials::model::BackgroundUsed;
  using scitbx::matrix::inversion_in_place;

  /**
   * A class to do 2D background subtraction.
   */
  class PlaneModel {
  public:

    /**
     * Calculate the background plane
     * @param data The shoebox data
     * @param mask The shoebox mask
     * @param fraction The fraction to use in initial estimate
     * @param nsigma The threshold for outliers
     */
    PlaneModel(const af::const_ref<int, af::c_grid<2> > &data,
               af::ref<int, af::c_grid<2> > mask,
               double fraction,
               double nsigma)
        : a_(0.0),
          b_(0.0),
          c_(0.0),
          rmsd_(0.0),
          maxdiff_(0.0),
          noutlier_(0),
          nbackground_(0) {
      compute(data, mask, fraction, nsigma);
    }

    /**
     * @returns The plane constant a*x + b*y +c
     */
    double a() const {
      return a_;
    }

    /**
     * @returns The plane constant a*x + b*y +c
     */
    double b() const {
      return b_;
    }

    /**
     * @returns The plane constant a*x + b*y +c
     */
    double c() const {
      return c_;
    }

    /**
     * @returns The maximum difference from plane
     */
    double maxdiff() const {
      return maxdiff_;
    }

    /**
     * @returns The rmsd from the plane
     */
    double rmsd() const {
      return rmsd_;
    }

    /**
     * @returns The number of outliers
     */
    std::size_t noutlier() const {
      return noutlier_;
    }

    /**
     * @returns The number of pixels used as background
     */
    std::size_t nbackground() const {
      return nbackground_;
    }

  private:

    /**
     * Compare index array by pixel value
     */
    struct compare_pixel_value {
      af::const_ref<int> data_;
      compare_pixel_value(const af::const_ref<int> &data)
        : data_(data) {}
      bool operator()(std::size_t a, std::size_t b) {
        return data_[a] < data_[b];
      }
    };

    /**
     * Compute the background.
     */
    void compute(
        const af::const_ref<int, af::c_grid<2> > &data,
        af::ref<int, af::c_grid<2> > mask,
        double fraction, double nsigma) {

      // Ensure data is consistent
      DIALS_ASSERT(fraction > 0 && fraction <= 1);
      DIALS_ASSERT(nsigma > 0);
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));

      // Compute the initial mask using a subset of the available pixels
      compute_initial_mask(data, mask, fraction);

      // Compute the background plane using the subset of pixels
      compute_background(data, mask);

      // Apply a correction to the background based on the fraction of pixels
      // chosen for use in the background plane calculation.
      compute_correction(fraction);

      // Check all pixels for outliers against the initial plane
      compute_final_mask(data, mask, nsigma);

      // Compute the final background plane
      compute_background(data, mask);

      // Compute statistics about deviations from the plane
      compute_deviation(data, mask);
    }

    /**
     * Calculate the initial mask. Select the fraction of pixels with the lowest
     * intensity and then update the mask for those pixels.
     */
    void compute_initial_mask(
        const af::const_ref< int, af::c_grid<2> > &data,
        af::ref< int, af::c_grid<2> > mask,
        double fraction) {
      int code = Valid | Background;
      std::vector<std::size_t> index;
      index.reserve(data.size());
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & code) == code) {
          index.push_back(i);
        } else {
          mask[i] &= ~BackgroundUsed;
        }
      }
      DIALS_ASSERT(index.size() > 0);
      std::sort(index.begin(), index.end(), compare_pixel_value(data.as_1d()));
      std::size_t nactive = (std::size_t)std::floor(fraction * index.size() + 0.5);
      DIALS_ASSERT(nactive > 0 && nactive <= index.size());
      for (std::size_t i = 0; i < nactive; ++i) {
        mask[index[i]] |= BackgroundUsed;
      }
    }

    /**
     * Compute the background plane given shoebox data and mask
     */
    void compute_background(
        const af::const_ref< int, af::c_grid<2> > &data,
        const af::const_ref< int, af::c_grid<2> > &mask) {
      std::vector<double> A(3*3, 0);
      std::vector<double> B(3, 0);
      int hy = data.accessor()[0] / 2;
      int hx = data.accessor()[1] / 2;
      for (std::size_t j = 0; j < mask.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < mask.accessor()[1]; ++i) {
          if (mask(j,i) & BackgroundUsed) {
            double x = ((int)i - hx);
            double y = ((int)j - hy);
            double p = data(j,i);
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
      std::cout << A[0] << std::endl;
      inversion_in_place(&A[0], 3, &B[0], 1);
      c_ = B[0];
      a_ = B[1];
      b_ = B[2];
      std::cout << a_ << ", " << b_ << ", " << c_ << std::endl;
    }

    /**
     * Compute the correction to the background plane for the initial estimate
     * where only a fraction of the input pixels are given.
     */
    void compute_correction(double fraction) {
      double factor = 0.8;
      if (fraction >= 0.999) {
        factor = 0.0;
      } else if (fraction >= 0.9) {
        factor = 0.05;
      } else if (fraction >= 0.7) {
        factor = 0.2;
      } else if (fraction >= 0.5) {
        factor = 0.4;
      } else if (fraction >= 0.3) {
        factor = 0.6;
      }
      c_ = c_ + factor * std::sqrt(std::abs(c_));
    }

    /**
     * Compute the final mask. Check the data for outliers by looking at the
     * deviation of each point from the plane.
     */
    void compute_final_mask(
        const af::const_ref< int, af::c_grid<2> > &data,
        af::ref< int, af::c_grid<2> > mask,
        double nsigma) {
      noutlier_ = 0;
      nbackground_ = 0;
      double max_background = nsigma * std::sqrt(std::abs(c_));
      int code = Valid | Background;
      int hy = data.accessor()[0] / 2;
      int hx = data.accessor()[1] / 2;
      for (std::size_t j = 0; j < mask.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < mask.accessor()[1]; ++i) {
          if ((mask(j,i) & code) == code) {
            double x = ((int)i - hx);
            double y = ((int)j - hy);
            double p1 = data(j,i);
            double p2 = a_ * x + b_ * y + c_;
            double d = std::abs((p1 - p2)*(p1 - p2));
            if (d > max_background) {
              mask(j,i) &= ~BackgroundUsed;
              noutlier_++;
            } else {
              mask(j,i) |= BackgroundUsed;
              nbackground_++;
            }
          }
        }
      }
    }

    /**
     * Compute the deviation and rmsd
     */
    void compute_deviation(
        const af::const_ref< int, af::c_grid<2> > &data,
        const af::const_ref< int, af::c_grid<2> > &mask) {
      rmsd_ = 0;
      maxdiff_ = 0;
      int hy = data.accessor()[0] / 2;
      int hx = data.accessor()[1] / 2;
      for (std::size_t j = 0; j < mask.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < mask.accessor()[1]; ++i) {
          int mask_value = mask(j,i);
          if (mask_value & BackgroundUsed) {
            double x = ((int)i - hx);
            double y = ((int)j - hy);
            double p1 = data(j,i);
            double p2 = a_ * x + b_ * y + c_;
            double d = (p1 - p2)*(p1 - p2);
            double sd = d * d;
            rmsd_ += sd;
            maxdiff_ = std::max(std::abs(d), maxdiff_);
          }
        }
      }
      rmsd_ = std::sqrt(rmsd_ / nbackground_);
    }

    double a_, b_, c_;
    double rmsd_, maxdiff_;
    std::size_t noutlier_;
    std::size_t nbackground_;
  };

}}} // namespace dials::algorithms::background

#endif // DIALS_ALGORITHMS_BACKGROUND_PLANE_H
