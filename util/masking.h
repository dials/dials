/*
 * masking.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_UTIL_MASKING_H
#define DIALS_UTIL_MASKING_H

#include <algorithm>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/panel.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace util {

  using dxtbx::model::BeamBase;
  using dxtbx::model::Panel;
  using scitbx::vec2;
  using scitbx::vec3;

  /**
   * A class to mask multiple resolution ranges
   */
  class ResolutionMaskGenerator {
  public:
    /**
     * Initialise the resolution at each pixel
     * @param beam The beam model
     * @param panel The panel model
     */
    ResolutionMaskGenerator(const BeamBase &beam, const Panel &panel)
        : resolution_(
            af::c_grid<2>(panel.get_image_size()[1], panel.get_image_size()[0])) {
      vec3<double> s0 = beam.get_s0();
      for (std::size_t j = 0; j < resolution_.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < resolution_.accessor()[1]; ++i) {
          vec2<double> px(i + 0.5, j + 0.5);
          try {
            resolution_(j, i) = panel.get_resolution_at_pixel(s0, px);
          } catch (dxtbx::error) {
            // Known failure: resolution at beam center is undefined
            resolution_(j, i) = 0.0;
          }
        }
      }
    }

    /**
     * Apply the mask
     * @param mask The mask
     * @param d_min The high resolution of the range
     * @param d_max The low resolution of the range
     */
    void apply(af::ref<bool, af::c_grid<2> > mask, double d_min, double d_max) const {
      DIALS_ASSERT(d_min < d_max);
      DIALS_ASSERT(resolution_.accessor()[0] == mask.accessor()[0]);
      DIALS_ASSERT(resolution_.accessor()[1] == mask.accessor()[1]);
      for (std::size_t j = 0; j < resolution_.accessor()[0]; ++j) {
        for (std::size_t i = 0; i < resolution_.accessor()[1]; ++i) {
          double d = resolution_(j, i);
          if (d_min <= d && d <= d_max) {
            mask(j, i) = false;
          }
        }
      }
    }

  private:
    af::versa<double, af::c_grid<2> > resolution_;
  };

}}  // namespace dials::util

#endif /* DIALS_UTIL_MASKING_H */
