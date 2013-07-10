/*
 * rebin_pixels.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_REFLEXION_BASIS_REBIN_PIXELS_H
#define DIALS_ALGORITHMS_REFLEXION_BASIS_REBIN_PIXELS_H

#include <cmath>
#include <scitbx/vec2.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/algorithms/polygon/clip/clip.h>
#include <dials/algorithms/polygon/area.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace reflexion_basis {

  using std::floor;
  using std::ceil;
  using scitbx::vec2;
  using scitbx::af::flex_double;
  using scitbx::af::double4;
  using scitbx::af::max;
  using scitbx::af::min;
  using dials::algorithms::polygon::simple_area;
  using dials::algorithms::polygon::clip::vert4;
  using dials::algorithms::polygon::clip::vert8;
  using dials::algorithms::polygon::clip::quad_with_convex_quad;

  typedef scitbx::af::flex<vec2<double> >::type flex_vec2_double;

  /**
   * Rebin pixels onto a regular grid
   * @param output The output grid
   * @param input The input grid
   * @param inputxy The input x/y coordinates
   */
  inline
  void rebin_pixels(flex_double &output, const flex_double &input,
                    const flex_vec2_double &inputxy) {

    // Check the sizes
    DIALS_ASSERT(output.accessor().all().size() == 2);
    DIALS_ASSERT(input.accessor().all().size() == 2);
    DIALS_ASSERT(inputxy.accessor().all().size() == 2);
    DIALS_ASSERT(inputxy.accessor().all()[0] == input.accessor().all()[0] + 1);
    DIALS_ASSERT(inputxy.accessor().all()[1] == input.accessor().all()[1] + 1);

    // Get the input and output sizes
    std::size_t output_height = output.accessor().all()[0];
    std::size_t output_width = output.accessor().all()[1];
    std::size_t input_height = input.accessor().all()[0];
    std::size_t input_width = input.accessor().all()[1];

    // Loop through all the input pixels
    for (std::size_t j = 0; j < input_height; ++j) {
      for (std::size_t i = 0; i < input_width; ++i) {

        // Get the x, y coords of the target point
        vec2<double> ixy00 = inputxy(j, i);
        vec2<double> ixy01 = inputxy(j, i+1);
        vec2<double> ixy10 = inputxy(j+1, i);
        vec2<double> ixy11 = inputxy(j+1, i+1);

        // Create the target polygon and calculate its area
        vert4 target(ixy00, ixy01, ixy11, ixy10);
        double target_area = simple_area(target);
        double value = input(j, i);
        DIALS_ASSERT(target_area > 0.0);

        // Get the range of new grid points
        double4 ix(ixy00[0], ixy01[0], ixy10[0], ixy11[0]);
        double4 iy(ixy00[1], ixy01[1], ixy10[1], ixy11[1]);
        int ox0 = (int)floor(min(ix.const_ref()));
        int oy0 = (int)floor(min(iy.const_ref()));
        int ox1 = (int)ceil(max(ix.const_ref()));
        int oy1 = (int)ceil(max(iy.const_ref()));

        // Cap the coordinates within the the output grid
        if (ox0 < 0) ox0 = 0;
        if (oy0 < 0) oy0 = 0;
        if (ox1 > output_width) ox1 = output_width;
        if (oy1 > output_height) oy1 = output_height;

        // Loop over all the pixels within the pixel range
        for (std::size_t jj = oy0; jj < oy1; ++jj) {
          for (std::size_t ii = ox0; ii < ox1; ++ii) {

            // Create the subject polygon
            vert4 subject(vec2<double>(ii, jj),
                          vec2<double>(ii+1, jj),
                          vec2<double>(ii+1,jj+1),
                          vec2<double>(ii,jj+1));

            // clip the polygon with the target polygon and calculate the
            // fraction of the area of the clipped polygon against the target.
            // Then redistribute the values from the target grid to the subject.
            vert8 result = quad_with_convex_quad(subject, target);
            double result_area = simple_area(result);
            double fraction = result_area / target_area;
            output(jj, ii) += fraction * value;
          }
        }
      }
    }
  }

}}} // dials::algorithms::reflexion_basis

#endif /* DIALS_ALGORITHMS_REFLEXION_BASIS_REBIN_PIXELS_H */
