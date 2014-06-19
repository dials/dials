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

  class PlaneModel {
  public:

    struct compare_pixel_value {
      af::const_ref<int> data_;
      compare_pixel_value(const af::const_ref<int> &data)
        : data_(data) {}
      bool operator()(std::size_t a, std::size_t b) {
        return data_[a] < data_[b];
      }
    };

    PlaneModel(const af::const_ref<int, af::c_grid<2> > &data,
               af::ref<int, af::c_grid<2> > mask,
               double fraction,
               double nsigma) {

      // Ensure data is consistent
      DIALS_ASSERT(fraction > 0 && fraction <= 1);
      DIALS_ASSERT(nsigma > 0);
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));

      // Add indices of good background pixels to array
      int code = Valid | Background;
      std::vector<std::size_t> index;
      index.reserve(data.size());
      for (std::size_t i = 0; i < index.size(); ++i) {
        if ((mask[i] & code) == code) {
          index.push_back(i);
        } else {
          mask[i] *= ~BackgroundUsed;   
        }
      }

      // Sort indices by ascending pixel value
      std::sort(index.begin(), index.end(), compare_pixel_value(data.as_1d()));
      std::size_t nactive = (std::size_t)(fraction * index.size());
      DIALS_ASSERT(nactive > 0);
      for (std::size_t i = 0; i < nactive; ++i) {
        mask[index[i]] |= BackgroundUsed;
      }

      // Compute the background plane
      std::vector<double> A(3*3, 0);
      std::vector<double> B(3, 0);
      for (std::size_t j = 0; j < mask.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < mask.accessor()[1]; ++i) {
          if (mask(j,i) & BackgroundUsed) {
            double x = (i + 0.5);
            double y = (j + 0.5);
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
      inversion_in_place(&A[0], 3, &B[0], 1);
      a_ = B[0];
      b_ = B[1];
      c_ = B[2];

      // The correction factor for different fractions used assuming gaussian
      // distribution of background counts
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
      c_ = c_ + factor * std::sqrt(c_);

      // Compute the rmsd of all background and active pixels.
      // Compute the maximum difference from the background plane
      double rmsd_all = 0.0;
      double rmsd_sub = 0.0;
      double max_diff = 0.0;
      for (std::size_t j = 0; j < mask.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < mask.accessor()[1]; ++i) {
          int mask_value = mask(j,i);
          if ((mask_value & code) == code) {
            double x = (i + 0.5);
            double y = (j + 0.5);
            double p1 = data(j,i);
            double p2 = a_ * x + b_ * y + c_;
            double d = (p1 - p2)*(p1 - p2); 
            double sd = d * d;
            rmsd_all += sd; 
            if (mask_value & BackgroundUsed) {
              rmsd_sub += sd;
            }
            max_diff = std::max(std::abs(d), max_diff);
          }
        }
      }
      rmsd_all = std::sqrt(rmsd_all / index.size());
      rmsd_sub = std::sqrt(rmsd_sub / nactive);

      // Compute the maximum expected background
      double max_background = nsigma * std::sqrt(c_);
      std::size_t noutlier = 0;
      if (max_diff > max_background) {
        for (std::size_t j = 0; j < mask.accessor()[0]; ++j) {
          for (std::size_t i = 0; i < mask.accessor()[1]; ++i) {
            int mask_value = mask(j,i);
            if ((mask_value & code) == code) {
              double x = (i + 0.5);
              double y = (j + 0.5);
              double p1 = data(j,i);
              double p2 = a_ * x + b_ * y + c_;
              double d = std::abs((p1 - p2)*(p1 - p2)); 
              if (d > max_background) {
                mask(j,i) &= ~BackgroundUsed;
                noutlier++;
              }
            }
          }
        }
      }
    }

    double a() const {
      return a_;
    }

    double b() const {
      return b_;
    }
    
    double c() const {
      return c_;
    }

    double rmsd() const {
      return rmsd_;
    }

  private:

    double a_, b_, c_, rmsd_;
  };

}}} // namespace dials::algorithms::background

#endif // DIALS_ALGORITHMS_BACKGROUND_PLANE_H
