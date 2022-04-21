/*
 * overload_checker.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SHOEBOX_OVERLOAD_CHECKER_H
#define DIALS_ALGORITHMS_SHOEBOX_OVERLOAD_CHECKER_H

#include <vector>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace shoebox {

  using dials::model::Shoebox;

  /**
   * A class to check for and mark overloaded pixels
   */
  class OverloadChecker {
    /**
     * The internal class used to do the checking
     */
    class Checker {
    public:
      typedef Shoebox<>::float_type float_type;

      /**
       * @param overload - The overload values for each pixel
       */
      Checker(const af::const_ref<double> &overload)
          : overload_(overload.begin(), overload.end()) {}

      /**
       * Mark all pixels that are overloaded as invalid
       * @param panel The panel number
       * @param data The data array
       * @param mask The mask array
       */
      bool operator()(std::size_t panel,
                      const af::const_ref<float_type, af::c_grid<3> > &data,
                      af::ref<int, af::c_grid<3> > mask) const {
        DIALS_ASSERT(panel < overload_.size());
        DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
        DIALS_ASSERT(data.size() == mask.size());
        bool result = false;
        double overload_value = overload_[panel];
        for (std::size_t i = 0; i < data.size(); ++i) {
          if (data[i] >= overload_value) {
            mask[i] &= ~Valid;
            result = true;
          }
        }
        return result;
      }

    private:
      af::shared<double> overload_;
    };

  public:
    /**
     * Add the overload values for this detector
     * @param overload The overloads for each panel
     */
    void add(const af::const_ref<double> &overload) {
      checker_.push_back(Checker(overload));
    }

    /**
     * Check each shoebox to see if it contains overloads
     * @param id The experiment id
     * @param shoebox The shoebox data
     * @returns flex.bool True contains outliers
     */
    af::shared<bool> operator()(const af::const_ref<int> id,
                                af::ref<Shoebox<> > shoebox) const {
      DIALS_ASSERT(id.size() == shoebox.size());
      af::shared<bool> result(id.size(), false);
      for (std::size_t i = 0; i < id.size(); ++i) {
        DIALS_ASSERT(id[i] >= 0);
        DIALS_ASSERT(id[i] < checker_.size());
        result[i] = checker_[id[i]](
          shoebox[i].panel, shoebox[i].data.const_ref(), shoebox[i].mask.ref());
      }
      return result;
    }

  private:
    std::vector<Checker> checker_;
  };

}}}  // namespace dials::algorithms::shoebox

#endif  // DIALS_ALGORITHMS_SHOEBOX_OVERLOAD_CHECKER_H
