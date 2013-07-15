/*
 * map_pixels.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_REFLEXION_BASIS_MAP_PIXELS_H
#define DIALS_ALGORITHMS_REFLEXION_BASIS_MAP_PIXELS_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/flex_types.h>
#include <dials/algorithms/reflexion_basis/rebin_pixels.h>
#include <dials/algorithms/reflexion_basis/coordinate_system.h>

namespace dials { namespace algorithms { namespace reflexion_basis {
    namespace transform {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;

  typedef scitbx::af::flex< vec3<double> >::type flex_vec3_double;

  /**
   * A class to generate local reciprocal space coordinates
   */
  class GridIndexGenerator {
  public:

    /**
     * Initialise the class
     * @param cs The coordinate system
     * @param bbox The bounding box
     * @param step_size The step size
     * @param grid_half_size The grid size
     * @param s1_map The map of detector s1 vectors.
     */
    GridIndexGenerator(const CoordinateSystem &cs, int6 bbox,
        vec2<double> step_size, std::size_t grid_half_size,
        const flex_vec3_double &s1_map)
      : s1_(cs.s1()),
        e1_(cs.e1_axis() / s1_.length()),
        e2_(cs.e2_axis() / s1_.length()),
        x0_(bbox[0]),
        y0_(bbox[2]),
        step_size_(step_size),
        grid_half_size_(grid_half_size),
        s1_map_(s1_map) {}

    /**
     * Get the e1/e2 coordinate at i, j
     * @param j The slow image index
     * @param i The fast image index
     * @returns The e1/e2 coordinate
     */
    vec2<double> operator()(int j, int i) const {
      double xx = (double)(x0_ + i);
      double yy = (double)(y0_ + j);
      vec3<double> ds = s1_map_(yy, xx) - s1_;
      double c1 = e1_ * ds;
      double c2 = e2_ * ds;
      double gi = grid_half_size_ + c1 / step_size_[0];
      double gj = grid_half_size_ + c2 / step_size_[1];
      return vec2<double>(gj, gi);
    }

  private:
    vec3<double> s1_;
    vec3<double> e1_;
    vec3<double> e2_;
    int x0_, y0_;
    vec2<double> step_size_;
    std::size_t grid_half_size_;
    flex_vec3_double s1_map_;
  };

  class MapPixels {
  public:

    MapPixels(const flex_vec3_double &s1_map, std::size_t grid_half_size,
        vec2<double> step_size)
      : s1_map_(s1_map),
        grid_half_size_(grid_half_size),
        step_size_(step_size) {
      DIALS_ASSERT(step_size_[0] > 0.0 && step_size_[1] > 0.0);
    }

    void operator()(const CoordinateSystem &cs, int6 bbox,
        const flex_double &image, flex_double &grid) {

      // Check array sizes
      DIALS_ASSERT(grid.accessor().all().size() == 3);
      DIALS_ASSERT(grid.accessor().all().all_eq(2 * grid_half_size_ + 1));
      DIALS_ASSERT(image.accessor().all().size() == 2);
      DIALS_ASSERT(image.accessor().all()[0] + 1 == s1_map_.accessor().all()[0]);
      DIALS_ASSERT(image.accessor().all()[1] + 1 == s1_map_.accessor().all()[1]);

      // Create the index generator for each coordinate of the bounding box
      GridIndexGenerator index(cs, bbox, step_size_, grid_half_size_, s1_map_);

      // Rebin the pixels to the grid
      rebin_pixels(grid, image, index);
    }

  private:

    flex_vec3_double s1_map_;
    std::size_t grid_half_size_;
    vec2<double> step_size_;
  };

}}}} // namespace dials::algorithms::reflexion_basis::transform

#endif /* DIALS_ALGORITHMS_REFLEXION_BASIS_MAP_PIXELS_H */
