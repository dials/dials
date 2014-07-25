/*
 * mask_foreground_2d.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SHOEBOX_MASK_FOREGROUND_2D_H
#define DIALS_ALGORITHMS_SHOEBOX_MASK_FOREGROUND_2D_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>
#include <dials/algorithms/reflection_basis/beam_vector_map.h>
#include <dials/algorithms/shoebox/mask_code.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace shoebox {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using scitbx::af::int4;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dials::model::Shoebox;
  using dials::algorithms::reflection_basis::CoordinateSystem2d;
  using dials::algorithms::reflection_basis::transform::beam_vector_map;

  /**
   * A class to mask foreground/background pixels
   */
  class MaskForeground2d {
  public:

    /**
     * Initialise the stuff needed to create the mask.
     * @param beam The beam model
     * @param detector The detector model
     * @param delta_b nsigma * sigma_divergence
     */
    MaskForeground2d(const Beam &beam, const Detector &detector,
                     double delta_b)
      : s1_map_(beam_vector_map(detector, beam, true)),
        s0_(beam.get_s0()) {
      DIALS_ASSERT(delta_b > 0.0);
      delta_b_r_ = 1.0 / delta_b;
    }

    /**
     * Set all the foreground/background pixels in the shoebox mask.
     * @param mask The mask
     * @param bbox The bounding box
     * @param s1 The beam vector
     */
    void operator()(const af::ref<int, af::c_grid<2> > &mask,
                    int4 bbox, vec3<double> s1) const {

      // Get some bits from the shoebox
      int x0 = bbox[0], x1 = bbox[1];
      int y0 = bbox[2], y1 = bbox[3];
      int xsize = x1 - x0;
      int ysize = y1 - y0;

      double delta_b_r2 = delta_b_r_ * delta_b_r_;

      // Check the size of the mask
      DIALS_ASSERT(mask.accessor()[0] == ysize);
      DIALS_ASSERT(mask.accessor()[1] == xsize);

      // Create the coordinate system and generators
      CoordinateSystem2d cs(s0_, s1);

      // Get the size of the image
      std::size_t width = s1_map_.accessor()[1];
      std::size_t height = s1_map_.accessor()[0];

      // Loop through all the pixels in the shoebox, transform the point
      // to the reciprocal space coordinate system and check that it is
      // within the ellipse defined by:
      // (c1 / delta_b)^2 + (c2 / delta_b)^2 <= 1
      // Mark those points within as Foreground and those without as
      // Background.
      for (int j = 0; j < ysize; ++j) {
        for (int i = 0; i < xsize; ++i) {
          if (x0 + i >= 0 && y0 + j >= 0 &&
              x0 + i < width && y0 + j < height) {
            vec2<double> gxy1 = cs.from_beam_vector(s1_map_(y0 + j, x0 + i));
            vec2<double> gxy2 = cs.from_beam_vector(s1_map_(y0 + j + 1, x0 + i));
            vec2<double> gxy3 = cs.from_beam_vector(s1_map_(y0 + j, x0 + i + 1));
            vec2<double> gxy4 = cs.from_beam_vector(s1_map_(y0 + j + 1, x0 + i + 1));
            double dxy1 = (gxy1[0]*gxy1[0] + gxy1[1]*gxy1[1]) * delta_b_r2;
            double dxy2 = (gxy2[0]*gxy2[0] + gxy2[1]*gxy2[1]) * delta_b_r2;
            double dxy3 = (gxy3[0]*gxy3[0] + gxy3[1]*gxy3[1]) * delta_b_r2;
            double dxy4 = (gxy4[0]*gxy4[0] + gxy4[1]*gxy4[1]) * delta_b_r2;
            double dxy = std::min(std::min(dxy1, dxy2), std::min(dxy3, dxy4));
            mask(j, i) |= (dxy <= 1.0) ? Foreground : Background;
          }
        }
      }
    }

  private:
    af::versa< vec3<double>, af::c_grid<2> > s1_map_;
    vec3<double> s0_;
    double delta_b_r_;
  };

}}} // namespace dials::algorithms::shoebox

#endif /* DIALS_ALGORITHMS_SHOEBOX_MASK_FOREGROUND_2D_H */
