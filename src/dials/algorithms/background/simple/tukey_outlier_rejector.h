/*
 * tukey_outlier_rejector.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_TUKEY_OUTLIER_REJECTOR_H
#define DIALS_ALGORITHMS_BACKGROUND_TUKEY_OUTLIER_REJECTOR_H

#include <algorithm>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/math/mean_and_variance.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace background {

  using scitbx::math::mean_and_variance;

  /**
   * Remove pixels with intensity <|> bounds
   */
  class TukeyOutlierRejector : public OutlierRejector {
  public:
    /**
     * Set the parameters 1.5 for both is Tukeys rule
     * @param lower The lower IQR multiplier
     * @param upper The upper IQR multiplier
     */
    TukeyOutlierRejector(double lower, double upper) : lower_(lower), upper_(upper) {
      DIALS_ASSERT(0 <= lower);
      DIALS_ASSERT(0 <= upper);
    }

    /**
     * @params shoebox The shoebox profile
     * @params mask The shoebox mask
     */
    virtual void mark(const af::const_ref<double, af::c_grid<3> > &shoebox,
                      af::ref<int, af::c_grid<3> > mask) const {
      const int mask_code = shoebox::Valid | shoebox::Background;

      // Ensure data is correctly sized.
      DIALS_ASSERT(shoebox.size() == mask.size());

      // Copy valid pixels and indices into list
      af::shared<double> data;
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code
            && (mask[i] & shoebox::Overlapped) == 0) {
          data.push_back(shoebox[i]);
        }
      }

      // Compute interquartile range
      std::sort(data.begin(), data.end());
      DIALS_ASSERT(data.size() > 2);
      std::size_t mid = data.size() / 2;
      std::size_t q1i = mid / 2;
      std::size_t q3i = mid + (data.size() - mid) / 2;
      DIALS_ASSERT(q1i < mid && mid < q3i && q3i < data.size());
      double q1 = data[q1i];
      double q3 = data[q3i];
      DIALS_ASSERT(q3 >= q1);
      double iqr = q3 - q1;
      double lower_bound = q1 - lower_ * iqr;
      double upper_bound = q3 + upper_ * iqr;

      // Set rejected pixels as 'not background'
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code
            && (mask[i] & shoebox::Overlapped) == 0) {
          if (lower_bound <= shoebox[i] && shoebox[i] <= upper_bound) {
            mask[i] |= shoebox::BackgroundUsed;
          }
        }
      }
    }

  private:
    double lower_, upper_;
  };

}}}  // namespace dials::algorithms::background

#endif /* DIALS_ALGORITHMS_BACKGROUND_TUKEY_OUTLIER_REJECTOR_H */
