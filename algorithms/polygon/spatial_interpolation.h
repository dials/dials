/*
 * spatial_interpolation.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_POLYGON_SPATIAL_INTERPOLATION_H
#define DIALS_ALGORITHMS_POLYGON_SPATIAL_INTERPOLATION_H

#include <cmath>
#include <scitbx/vec2.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/simple_io.h>
#include <dials/algorithms/polygon/clip/clip.h>
#include <dials/algorithms/polygon/area.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace polygon {
  namespace spatial_interpolation {

    using dials::algorithms::polygon::simple_area;
    using dials::algorithms::polygon::clip::quad_with_convex_quad;
    using dials::algorithms::polygon::clip::vert4;
    using dials::algorithms::polygon::clip::vert8;
    using scitbx::vec2;
    using scitbx::af::double4;
    using scitbx::af::int2;
    using scitbx::af::int4;
    using scitbx::af::max;
    using scitbx::af::min;
    using std::ceil;
    using std::floor;

    /**
     * Struct to hold indices of matches
     */
    struct Match {
      int in;
      int out;
      double fraction;
      Match() : in(0), out(0), fraction(0.0) {}
      Match(int in_, int out_, double f_) : in(in_), out(out_), fraction(f_) {}
    };

    /**
     * Helper function to calculate range of quad on grid
     * @param input The input Quad
     * @param output_size The output size
     * @returns The grid indices (x0, x1, y0, y1)
     */
    inline int4 quad_grid_range(const vert4 &input, af::c_grid<2> output_size) {
      // Get the range of new grid points
      double4 ix(input[0][0], input[1][0], input[2][0], input[3][0]);
      double4 iy(input[0][1], input[1][1], input[2][1], input[3][1]);
      int ox0 = (int)floor(min(ix.const_ref()));
      int oy0 = (int)floor(min(iy.const_ref()));
      int ox1 = (int)ceil(max(ix.const_ref()));
      int oy1 = (int)ceil(max(iy.const_ref()));

      // Cap the coordinates within the the output grid
      if (ox0 < 0) ox0 = 0;
      if (oy0 < 0) oy0 = 0;
      if (ox1 > output_size[1]) ox1 = output_size[1];
      if (oy1 > output_size[0]) oy1 = output_size[0];
      return int4(ox0, ox1, oy0, oy1);
    }

    /**
     * Get the area of the polygon and reverse if it's backwards
     * @param input The input polygon
     * @returns The area
     */
    inline double reverse_quad_inplace_if_backward(vert4 &input) {
      double area = simple_area(input);
      DIALS_ASSERT(area != 0.0);
      if (area < 0.0) {
        std::swap(input[0], input[3]);
        std::swap(input[1], input[2]);
        area = -area;
      }
      return area;
    }

    /**
     * Get the intersection of a quad with a regular grid point
     * @param a The quad
     * @param i The fast grid index
     * @param j The slow grid index
     * @returns The area
     */
    inline double quad_grid_intersection_area(const vert4 &a, int i, int j) {
      vert4 b(vec2<double>(i, j),
              vec2<double>(i + 1, j),
              vec2<double>(i + 1, j + 1),
              vec2<double>(i, j + 1));
      return simple_area(quad_with_convex_quad(a, b));
    }

    /**
     * Get the overlaps between an input quad and the grid points
     * @param input The quad
     * @param output_size The size of the output grid
     * @param index The index of the quad
     * @returns The matches between the quad and the grid
     */
    inline af::shared<Match> quad_to_grid(vert4 input,
                                          af::c_grid<2> output_size,
                                          int index) {
      af::shared<Match> matches;
      int4 range = quad_grid_range(input, output_size);
      if (range[0] >= range[1] || range[2] >= range[3]) return matches;
      double target_area = reverse_quad_inplace_if_backward(input);
      for (std::size_t jj = range[2]; jj < range[3]; ++jj) {
        for (std::size_t ii = range[0]; ii < range[1]; ++ii) {
          double result_area = quad_grid_intersection_area(input, ii, jj);
          if (result_area > 0) {
            double fraction = result_area / target_area;
            matches.push_back(Match(index, ii + jj * output_size[1], fraction));
          }
        }
      }
      return matches;
    }

    /**
     * Get the overlaps between and input grid and output quad
     * @param output The quad
     * @param input_size The size of the output grid
     * @param index The index of the quad
     * @returns The matches between the quad and the grid
     */
    inline af::shared<Match> grid_to_quad(vert4 output,
                                          af::c_grid<2> input_size,
                                          int index) {
      af::shared<Match> matches;
      int4 range = quad_grid_range(output, input_size);
      if (range[0] >= range[1] || range[2] >= range[3]) return matches;
      reverse_quad_inplace_if_backward(output);
      for (std::size_t jj = range[2]; jj < range[3]; ++jj) {
        for (std::size_t ii = range[0]; ii < range[1]; ++ii) {
          double result_area = quad_grid_intersection_area(output, ii, jj);
          if (result_area > 0) {
            double fraction = result_area;
            matches.push_back(Match(ii + jj * input_size[1], index, fraction));
          }
        }
      }
      return matches;
    }

    /**
     * Get the matches between an input irregular grid and an output regular
     * grid.
     * @param inputxy The input x/y coordinates
     * @param output_size The size of the output grid
     * @returns The matches between grid points
     */
    inline af::shared<Match> irregular_grid_to_grid(
      const af::const_ref<vec2<double>, af::c_grid<2> > &inputxy,
      af::tiny<std::size_t, 2> output_size) {
      af::shared<Match> matches;
      DIALS_ASSERT(inputxy.accessor().all_gt(0) && output_size.all_gt(0));
      for (std::size_t j = 0, k = 0; j < inputxy.accessor()[0] - 1; ++j) {
        for (std::size_t i = 0; i < inputxy.accessor()[1] - 1; ++i, ++k) {
          vert4 input(
            inputxy(j, i), inputxy(j, i + 1), inputxy(j + 1, i + 1), inputxy(j + 1, i));
          af::shared<Match> temp = quad_to_grid(input, output_size, k);
          std::copy(temp.begin(), temp.end(), std::back_inserter(matches));
        }
      }
      return matches;
    }

    /**
     * Get the matches between an output irregular grid and an input regular
     * grid.
     * @param outputxy The output x/y coordinates
     * @param input_size The size of the input grid
     * @returns The matches between grid points
     */
    inline af::shared<Match> grid_to_irregular_grid(
      const af::const_ref<vec2<double>, af::c_grid<2> > &outputxy,
      af::tiny<std::size_t, 2> input_size) {
      af::shared<Match> matches;
      DIALS_ASSERT(outputxy.accessor().all_gt(0) && input_size.all_gt(0));
      for (std::size_t j = 0, k = 0; j < outputxy.accessor()[0] - 1; ++j) {
        for (std::size_t i = 0; i < outputxy.accessor()[1] - 1; ++i, ++k) {
          vert4 output(outputxy(j, i),
                       outputxy(j, i + 1),
                       outputxy(j + 1, i + 1),
                       outputxy(j + 1, i));
          af::shared<Match> temp = grid_to_quad(output, input_size, k);
          std::copy(temp.begin(), temp.end(), std::back_inserter(matches));
        }
      }
      return matches;
    }

    /**
     * Regrid the input irregular grid onto a regular grid.
     * @param input The input irregular grid
     * @param inputxy The input grid coordinates
     * @param output_size The output grid size
     * @returns The output grid
     */
    inline af::versa<double, af::c_grid<2> > regrid_irregular_grid_to_grid(
      const af::const_ref<double, af::c_grid<2> > &input,
      const af::const_ref<vec2<double>, af::c_grid<2> > &inputxy,
      af::tiny<std::size_t, 2> output_size) {
      DIALS_ASSERT(inputxy.accessor()[0] == input.accessor()[0] + 1);
      DIALS_ASSERT(inputxy.accessor()[1] == input.accessor()[1] + 1);
      af::c_grid<2> accessor(output_size);
      af::versa<double, af::c_grid<2> > result(accessor, 0.0);
      af::shared<Match> matches = irregular_grid_to_grid(inputxy, output_size);
      for (std::size_t i = 0; i < matches.size(); ++i) {
        Match &m = matches[i];
        result[m.out] += input[m.in] * m.fraction;
      }
      return result;
    }

    /**
     * Regrid the input grid onto an irregular grid.
     * @param input The input regular grid
     * @param inputxy The output grid coordinates
     * @returns The output grid
     */
    inline af::versa<double, af::c_grid<2> > regrid_grid_to_irregular_grid(
      const af::const_ref<double, af::c_grid<2> > &input,
      const af::const_ref<vec2<double>, af::c_grid<2> > &outputxy) {
      af::c_grid<2> accessor(outputxy.accessor()[0] - 1, outputxy.accessor()[1] - 1);
      af::versa<double, af::c_grid<2> > result(accessor, 0.0);
      af::shared<Match> matches = grid_to_irregular_grid(outputxy, input.accessor());
      for (std::size_t i = 0; i < matches.size(); ++i) {
        Match &m = matches[i];
        result[m.out] += input[m.in] * m.fraction;
      }
      return result;
    }

}}}}  // namespace dials::algorithms::polygon::spatial_interpolation

#endif /* DIALS_ALGORITHMS_POLYGON_SPATIAL_INTERPOLATION_H */
