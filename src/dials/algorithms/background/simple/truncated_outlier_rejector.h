/*
 * truncated_outlier_rejector.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_BACKGROUND_TRUNCATED_OUTLIER_REJECTOR_H
#define DIALS_ALGORITHMS_BACKGROUND_TRUNCATED_OUTLIER_REJECTOR_H

#include <algorithm>
#include <scitbx/array_family/ref_reductions.h>
#include <boost/math/special_functions/erf.hpp>
#include <scitbx/math/mean_and_variance.h>
#include <dials/array_family/sort_index.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace background {

  using dials::af::sort_index;

  /**
   * Remove top and bottom n% of pixels to use in background
   */
  class TruncatedOutlierRejector : public OutlierRejector {
  public:
    /**
     * @param lower The lower percent to reject
     * @param upper The upper percent to reject
     */
    TruncatedOutlierRejector(double lower, double upper)
        : lower_(lower), upper_(upper) {
      DIALS_ASSERT(0 <= lower && lower <= 1.0);
      DIALS_ASSERT(0 <= upper && upper <= 1.0);
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
      af::shared<int> indices;
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        if ((mask[i] & mask_code) == mask_code
            && (mask[i] & shoebox::Overlapped) == 0) {
          indices.push_back(i);
        }
      }

      // Sort the pixels into ascending intensity order
      sort_index(indices.begin(), indices.end(), shoebox.begin());
      af::shared<double> pixels(indices.size(), af::init_functor_null<double>());
      for (std::size_t i = 0; i < indices.size(); ++i) {
        pixels[i] = (double)shoebox[indices[i]];
      }

      // Set rejected pixels as 'not background'
      std::size_t num_data = indices.size();
      std::size_t i0 = (std::size_t)(lower_ * num_data / 2.0);
      std::size_t i1 = num_data - (std::size_t)(upper_ * num_data / 2.0);
      for (std::size_t i = i0; i < i1; ++i) {
        mask[indices[i]] |= shoebox::BackgroundUsed;
      }
    }

  private:
    double lower_, upper_;
  };

}}}  // namespace dials::algorithms::background

#endif /* DIALS_ALGORITHMS_BACKGROUND_TRUNCATED_OUTLIER_REJECTOR_H */
