/*
 * index_generator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_INDEX_GENERATOR_H
#define DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_INDEX_GENERATOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <dials/algorithms/profile_model/gaussian_rs/coordinate_system.h>
#include <dials/array_family/scitbx_shared_and_versa.h>

namespace dials {
  namespace algorithms {
    namespace profile_model {
      namespace gaussian_rs {
  namespace transform {

    using scitbx::vec2;
    using scitbx::vec3;
    using scitbx::af::int6;

    /**
     * A class to generate local reciprocal space coordinates
     */
    class CoordinateGenerator {
    public:
      /**
       * Initialise the class
       * @param cs The coordinate system
       * @param x0 The x offset
       * @param y0 The y offset
       * @param s1_map The map of detector s1 vectors.
       */
      CoordinateGenerator(const CoordinateSystem &cs,
                          int x0,
                          int y0,
                          const af::versa<vec3<double>, af::c_grid<2> > &s1_map)
          : s1_(cs.s1()),
            e1_(cs.e1_axis() / s1_.length()),
            e2_(cs.e2_axis() / s1_.length()),
            x0_(x0),
            y0_(y0),
            s1_map_(s1_map) {}

      /**
       * Get the e1/e2 coordinate at i, j
       * @param j The slow image index
       * @param i The fast image index
       * @returns The e1/e2 coordinate
       */
      vec2<double> operator()(int j, int i) const {
        int xx = x0_ + i;
        int yy = y0_ + j;
        DIALS_ASSERT(yy >= 0 && xx >= 0);
        DIALS_ASSERT(yy < s1_map_.accessor()[0] && xx < s1_map_.accessor()[1]);
        vec3<double> ds = s1_map_(yy, xx) - s1_;
        double c1 = e1_ * ds;
        double c2 = e2_ * ds;
        return vec2<double>(c1, c2);
      }

    private:
      vec3<double> s1_;
      vec3<double> e1_;
      vec3<double> e2_;
      int x0_, y0_;
      af::versa<vec3<double>, af::c_grid<2> > s1_map_;
    };

    /**
     * A class to generate local reciprocal space coordinates
     */
    class GridIndexGenerator {
    public:
      /**
       * Initialise the class
       * @param cs The coordinate system
       * @param x0 The x offset
       * @param y0 The y offset
       * @param step_size The step size
       * @param grid_half_size The grid size
       * @param s1_map The map of detector s1 vectors.
       */
      GridIndexGenerator(const CoordinateSystem &cs,
                         int x0,
                         int y0,
                         vec2<double> step_size,
                         std::size_t grid_half_size,
                         const af::versa<vec3<double>, af::c_grid<2> > &s1_map)
          : coordinate_(cs, x0, y0, s1_map),
            step_size_(step_size),
            grid_half_size_(grid_half_size) {
        DIALS_ASSERT(step_size_.const_ref().all_gt(0.0));
        DIALS_ASSERT(grid_half_size_ > 0.0);
      }

      /**
       * Get the e1/e2 coordinate at i, j
       * @param j The slow image index
       * @param i The fast image index
       * @returns The e1/e2 coordinate
       */
      vec2<double> operator()(int j, int i) const {
        vec2<double> c12 = coordinate_(j, i);
        double gi = grid_half_size_ + c12[0] / step_size_[0] + 0.5;
        double gj = grid_half_size_ + c12[1] / step_size_[1] + 0.5;
        return vec2<double>(gi, gj);
      }

    private:
      CoordinateGenerator coordinate_;
      vec2<double> step_size_;
      std::size_t grid_half_size_;
    };

}}}}}  // namespace dials::algorithms::profile_model::gaussian_rs::transform

#endif /* DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_INDEX_GENERATOR_H */
