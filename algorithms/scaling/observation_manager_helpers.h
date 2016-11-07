/*
 * grouped_obs.h
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_SCALING_GROUPED_OBS_H
#define DIALS_ALGORITHMS_SCALING_GROUPED_OBS_H

#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>
#include <cctbx/miller.h>

// For debugging
// #include <scitbx/array_family/simple_io.h>
// then can do things like this:
// std::cout << something.const_ref() << "\n";

namespace dials { namespace scaling {

  typedef cctbx::miller::index<> miller_index;

  /**
   * Create a selection that keeps only Miller indices with at least the
   * specified number of occurrences. This function is designed to be simple
   * and fast. It assumes that the input array of Miller indices is sorted,
   * so that indices are already grouped in runs of adjacent equal values.
   */
  af::shared<bool>
  minimum_multiplicity_selection(af::const_ref< miller_index > group_index,
                                 const std::size_t multiplicity) {

    std::size_t nobs = group_index.size();
    af::shared<bool> sel(nobs, true);

    af::ref<bool> r = sel.ref();

    miller_index gp_idx = group_index[0];
    std::size_t gp_start = 0;
    std::size_t gp_end = 0;
    for (std::size_t i = 1; i < nobs; ++i) {
      miller_index h = group_index[i];
      if (h != gp_idx)
      {
        // the group ended at the previous element
        gp_end = i;

        // if too few members of this group, set the selection accordingly
        if ((gp_end - gp_start) < multiplicity)
        {
          for (std::size_t j = gp_start; j < gp_end; ++j) {
             r[j] = false;
          }
        }

        // reset current group
        gp_start = i;
        gp_idx = h;
      }
    }
    // final group
    if ((nobs - gp_start) < multiplicity)
    {
      for (std::size_t j = gp_start; j < nobs; ++j) {
         r[j] = false;
      }
    }

    return sel;
  }

  /**
   * Maintain information about groups of intensities and their scales so that
   * the weighted average intensity of the observations can be calculated
   * according to the formula of Hamilton, Rollett and Sparks (1965). Groups of
   * observations are determined by consecutive runs with equal Miller indices.
   * Therefore, it is expected that the input data are sorted by Miller index.
   * No filtering by multiplicity is done here - the minimum group size is 1.
   * As single-member groups are unlikely to be what is wanted, the input
   * should be filtered by some minimum multiplicity level first. The function
   * minimum_multiplicity_selection can be used for this purpose.
   */
  class GroupedObservations {
  public:

    // constructor
    GroupedObservations(af::const_ref< miller_index > group_index,
                        af::shared<double> intensity,
                        af::shared<double> weight,
                        af::shared<double> phi,
                        af::shared<double> scale)
      : group_index_(group_index), // initialisation lists
        intensity_(intensity),
        weight_(weight),
        phi_(phi),
        scale_(scale){
      nobs_ = group_index_.size();
      DIALS_ASSERT(nobs_ > 1);
      DIALS_ASSERT(intensity_.size() == nobs_);
      DIALS_ASSERT(weight_.size() == nobs_);
      DIALS_ASSERT(phi_.size() == nobs_);
      DIALS_ASSERT(scale_.size() == nobs_);

      //calculate the group lengths here
      miller_index gp_idx = group_index_[0];
      std::size_t gp_start = 0;
      for (std::size_t i = 1; i < nobs_; ++i) {
        miller_index h = group_index_[i];
        if (h != gp_idx)
        {
          // record previous group size
          group_size_.push_back(i - gp_start);

          // reset current group
          gp_start = i;
          gp_idx = h;
        }
      }
      // record final group size
      group_size_.push_back(nobs_ - gp_start);

      ngroups_ = group_size_.size();
    }

    // intensity accessors
    af::shared<double> get_intensity() const {
      return intensity_;
    }

    void set_intensity(af::shared<double> intensity) {
      DIALS_ASSERT(intensity.size() == nobs_);
      intensity_ = intensity;
    }

    // weight accessors
    af::shared<double> get_weight() const {
      return weight_;
    }

    void set_weight(af::shared<double> weight) {
      DIALS_ASSERT(weight.size() == nobs_);
      weight_ = weight;
    }

    // scale accessors
    af::shared<double> get_scale() const {
      return scale_;
    }

    void set_scale(af::shared<double> scale) {
      DIALS_ASSERT(scale.size() == nobs_);
      scale_ = scale;
    }

    // phi getter
    af::shared<double> get_phi() const {
      return phi_;
    }

    // group size getter
    af::shared<std::size_t> get_group_size() const {
      return group_size_;
    }

    // calculate and return average intensity in groups
    af::shared<double> get_average_intensity() {

      af::shared<double> avI(nobs_);
      af::ref<double> r = avI.ref();

      // loop over the groups
      std::size_t gp_start = 0;
      for (std::size_t i = 0; i < ngroups_; ++i){
        std::size_t sz = group_size_[i];

        double u = 0;
        double v = 0;

        // loop over obs in the group
        for (std::size_t j = gp_start; j < gp_start + sz; ++j){
          u += weight_[j] * scale_[j] * intensity_[j];
          v += weight_[j] * scale_[j] * scale_[j];
        }
        DIALS_ASSERT(v > 0);
        double gp_avI = u / v;

        // assign elements
        for (std::size_t j = gp_start; j < gp_start + sz; ++j){
          r[j] = gp_avI;
        }

        // update group start index
        gp_start += sz;
      }

      return avI;
    }

  private:
    af::const_ref<miller_index> group_index_;
    af::shared<double> intensity_;
    af::shared<double> weight_;
    af::shared<double> phi_;
    af::shared<double> scale_;

    af::shared<std::size_t> group_size_;

    std::size_t nobs_;
    std::size_t ngroups_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SCALING_GROUPED_OBS_H
