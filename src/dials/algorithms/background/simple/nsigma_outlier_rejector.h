/*
 * nsigma_outlier_rejector.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_NSIGMA_OUTLIER_REJECTOR_H
#define DIALS_ALGORITHMS_BACKGROUND_NSIGMA_OUTLIER_REJECTOR_H

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
  class NSigmaOutlierRejector : public OutlierRejector {
  public:
    /**
     * @param lower The lower n sigma
     * @param upper The upper n sigma
     */
    NSigmaOutlierRejector(double lower, double upper) : lower_(lower), upper_(upper) {
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

      // Compute the mean and sigma
      DIALS_ASSERT(data.size() > 1);

      mean_and_variance<double> mv(data.const_ref());
      double mean = mv.mean();
      double sigma = mv.unweighted_sample_standard_deviation();
      double p0 = mean - lower_ * sigma;
      double p1 = mean + upper_ * sigma;

      // Set rejected pixels as 'not background'
      for (std::size_t i = 0; i < mask.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code
            && (mask[i] & shoebox::Overlapped) == 0) {
          if (p0 <= shoebox[i] && shoebox[i] <= p1) {
            mask[i] |= shoebox::BackgroundUsed;
          }
        }
      }
    }

  private:
    double lower_, upper_;
  };

}}}  // namespace dials::algorithms::background

#endif /* DIALS_ALGORITHMS_BACKGROUND_NSIGMA_OUTLIER_REJECTOR_H */
