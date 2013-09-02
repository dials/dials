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
#include <dials/algorithms/reflection_basis/rebin_pixels.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>

namespace dials { namespace algorithms { namespace reflection_basis {
    namespace transform {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using scitbx::af::flex_double;
  using scitbx::af::flex_bool;
  using scitbx::af::flex_grid;

  typedef scitbx::af::flex< vec3<double> >::type flex_vec3_double;

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
    CoordinateGenerator(const CoordinateSystem &cs, int x0, int y0,
        const flex_vec3_double &s1_map)
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
    flex_vec3_double s1_map_;
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
    GridIndexGenerator(const CoordinateSystem &cs, int x0, int y0,
        vec2<double> step_size, std::size_t grid_half_size,
        const flex_vec3_double &s1_map)
      : coordinate_(cs, x0, y0, s1_map),
        step_size_(step_size),
        grid_half_size_(grid_half_size) {}

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


  /**
   * Class to map the pixels to the reflection basis grid
   */
  class MapPixelsForward {
  public:

    /**
     * Initialise the class
     * @param s1_map The beam vector map
     * @param grid_half_size The grid half size
     * @param step_size The grid step size
     */
    MapPixelsForward(const flex_vec3_double &s1_map, std::size_t grid_half_size,
        vec2<double> step_size)
      : s1_map_(s1_map),
        grid_half_size_(grid_half_size),
        step_size_(step_size) {
      DIALS_ASSERT(grid_half_size > 0);
      DIALS_ASSERT(step_size_[0] > 0.0 && step_size_[1] > 0.0);
      DIALS_ASSERT(s1_map.accessor().all().size() == 2);
      image_size_ = vec2<std::size_t>(s1_map.accessor().all()[0],
                                      s1_map.accessor().all()[1]);
    }

    /**
     * Map the pixels for a reflection
     * @param cs The coordinate system
     * @param bbox The bounding box
     * @param image The image array
     * @param mask The mask array
     * @param z_fraction The z fraction array
     * @returns grid The grid array
     */
    flex_double operator()(const CoordinateSystem &cs, int6 bbox,
        const flex_double &image, const flex_bool &mask,
        const flex_double &z_fraction) const {
      flex_double grid(flex_grid<>(2 * grid_half_size_ + 1,
                                   2 * grid_half_size_ + 1,
                                   2 * grid_half_size_ + 1), 0);
      this->operator()(cs, bbox, image, mask, z_fraction, grid);
      return grid;
    }

    /**
     * Map the pixels for a reflection
     * @param cs The coordinate system
     * @param bbox The bounding box
     * @param image The image array
     * @param background The background array
     * @param mask The mask array
     * @param z_fraction The z fraction array
     * @returns grid The grid array
     */
    std::pair<flex_double, flex_double> operator()(
        const CoordinateSystem &cs, int6 bbox,
        const flex_double &image, const flex_double &background,
        const flex_bool &mask, const flex_double &z_fraction) const {
      flex_double image_grid(flex_grid<>(
        2 * grid_half_size_ + 1,
        2 * grid_half_size_ + 1,
        2 * grid_half_size_ + 1), 0);
      flex_double background_grid(flex_grid<>(
        2 * grid_half_size_ + 1,
        2 * grid_half_size_ + 1,
        2 * grid_half_size_ + 1), 0);
      this->operator()(cs, bbox, image, mask, z_fraction, image_grid);
      return std::make_pair(image_grid, background_grid);
    }

    /**
     * Map the pixels for a reflection
     * @param cs The coordinate system
     * @param bbox The bounding box
     * @param image The image array
     * @param mask The mask array
     * @param z_fraction The z fraction array
     * @param grid The grid array
     */
    void operator()(const CoordinateSystem &cs, int6 bbox,
        const flex_double &image, const flex_bool &mask,
        const flex_double z_fraction, flex_double &grid) const {

      // Check array sizes
      DIALS_ASSERT(grid.accessor().all().size() == 3);
      DIALS_ASSERT(grid.accessor().all().all_eq(2 * grid_half_size_ + 1));
      DIALS_ASSERT(z_fraction.accessor().all()[0] == bbox[5] - bbox[4]);
      DIALS_ASSERT(z_fraction.accessor().all()[1] == 2 * grid_half_size_ + 1);
      DIALS_ASSERT(image.accessor().all().size() == 3);
      DIALS_ASSERT(image.accessor().all()[0] == bbox[5] - bbox[4]);
      DIALS_ASSERT(image.accessor().all()[1] == bbox[3] - bbox[2]);
      DIALS_ASSERT(image.accessor().all()[2] == bbox[1] - bbox[0]);
      DIALS_ASSERT(bbox[0] >= 0 && bbox[1] < image_size_[1]);
      DIALS_ASSERT(bbox[2] >= 0 && bbox[3] < image_size_[0]);

      // Create the index generator for each coordinate of the bounding box
      GridIndexGenerator index(cs, bbox[0], bbox[2], step_size_,
        grid_half_size_, s1_map_);

      // Get the input and output sizes
      std::size_t grid_depth = grid.accessor().all()[0];
      std::size_t grid_height = grid.accessor().all()[1];
      std::size_t grid_width = grid.accessor().all()[2];
      std::size_t image_depth = image.accessor().all()[0];
      std::size_t image_height = image.accessor().all()[1];
      std::size_t image_width = image.accessor().all()[2];

      // Loop through all the input pixels
      for (std::size_t j = 0; j < image_height; ++j) {
        for (std::size_t i = 0; i < image_width; ++i) {

          // Get the x, y coords of the target point
          vec2<double> ixy00 = index(j, i);
          vec2<double> ixy01 = index(j, i+1);
          vec2<double> ixy10 = index(j+1, i);
          vec2<double> ixy11 = index(j+1, i+1);

          // Get the range of new grid points
          double4 ix(ixy00[0], ixy01[0], ixy10[0], ixy11[0]);
          double4 iy(ixy00[1], ixy01[1], ixy10[1], ixy11[1]);
          int gx0 = (int)floor(min(ix.const_ref()));
          int gy0 = (int)floor(min(iy.const_ref()));
          int gx1 = (int)ceil(max(ix.const_ref()));
          int gy1 = (int)ceil(max(iy.const_ref()));

          // Cap the coordinates within the the output grid
          if (gx0 < 0) gx0 = 0;
          if (gy0 < 0) gy0 = 0;
          if (gx1 > grid_width) gx1 = grid_width;
          if (gy1 > grid_height) gy1 = grid_height;
          if (gx0 >= gx1 || gy0 >= gy1) {
            continue;
          }

          // Create the target polygon and calculate its area
          vert4 target(ixy00, ixy01, ixy11, ixy10);
          double target_area = simple_area(target);
          if (target_area < 0.0) {
            std::swap(target[0], target[3]);
            std::swap(target[1], target[2]);
            target_area = -target_area;
          } else if (target_area == 0.0) {
            continue;
          }

          // Loop over all the pixels within the pixel range
          for (std::size_t jj = gy0; jj < gy1; ++jj) {
            for (std::size_t ii = gx0; ii < gx1; ++ii) {

              // Create the subject polygon
              vert4 subject(vec2<double>(ii, jj),
                            vec2<double>(ii+1, jj),
                            vec2<double>(ii+1,jj+1),
                            vec2<double>(ii,jj+1));

              // clip the polygon with the target polygon and calculate the
              // fraction of the area of the clipped polygon against the target.
              // Then redistribute the values from the target grid to the
              // subject.
              vert8 result = quad_with_convex_quad(subject, target);
              double result_area = simple_area(result);
              double xy_fraction = result_area / target_area;

              // Copy the values to the grid
              for (int k = 0; k < image_depth; ++k) {
                if (mask(k, j, i) != 0) {
                  double value = image(k, j, i) * xy_fraction;
                  for (int kk = 0; kk < grid_depth; ++kk) {
                    grid(kk, jj, ii) += value * z_fraction(k, kk);
                  }
                }
              }
            }
          }
        }
      }
    }

  private:

    flex_vec3_double s1_map_;
    std::size_t grid_half_size_;
    vec2<double> step_size_;
    vec2<std::size_t> image_size_;
  };


  /**
   * Class to map the pixels from the reflection basis grid
   */
  class MapPixelsReverse {
  public:

    /**
     * Initialise the class
     * @param s1_map The beam vector map
     * @param grid_half_size The grid half size
     * @param step_size The grid step size
     */
    MapPixelsReverse(const flex_vec3_double &s1_map, std::size_t grid_half_size,
        vec2<double> step_size)
      : s1_map_(s1_map),
        grid_half_size_(grid_half_size),
        step_size_(step_size) {
      DIALS_ASSERT(grid_half_size > 0);
      DIALS_ASSERT(step_size_[0] > 0.0 && step_size_[1] > 0.0);
      image_size_ = vec2<std::size_t>(s1_map.accessor().all()[0],
                                      s1_map.accessor().all()[1]);
    }

    /**
     * Map the pixels for a reflection
     * @param cs The coordinate system
     * @param bbox The bounding box
     * @param grid The grid
     * @param z_fraction The z fraction array
     * @returns The mapped image
     */
    flex_double operator()(const CoordinateSystem &cs, int6 bbox,
        const flex_double &grid,
        const flex_double &z_fraction) const {
      flex_double image(flex_grid<>(bbox[5] - bbox[4],
                                    bbox[3] - bbox[2],
                                    bbox[1] - bbox[0]));
      this->operator()(cs, bbox, grid, z_fraction, image);
      return image;
    }

    /**
     * Map the pixels for a reflection
     * @param cs The coordinate system
     * @param bbox The bounding box
     * @param grid The grid
     * @param z_fraction The z fraction array
     * @param image The mapped image
     */
    void operator()(const CoordinateSystem &cs, int6 bbox,
        const flex_double &grid, const flex_double z_fraction,
        flex_double &image) const {

      // Check array sizes
      DIALS_ASSERT(grid.accessor().all().size() == 3);
      DIALS_ASSERT(grid.accessor().all().all_eq(2 * grid_half_size_ + 1));
      DIALS_ASSERT(z_fraction.accessor().all()[1] == bbox[5] - bbox[4]);
      DIALS_ASSERT(z_fraction.accessor().all()[0] == 2 * grid_half_size_ + 1);
      DIALS_ASSERT(image.accessor().all().size() == 3);
      DIALS_ASSERT(image.accessor().all()[0] == bbox[5] - bbox[4]);
      DIALS_ASSERT(image.accessor().all()[1] == bbox[3] - bbox[2]);
      DIALS_ASSERT(image.accessor().all()[2] == bbox[1] - bbox[0]);
      DIALS_ASSERT(bbox[0] >= 0 && bbox[1] <= image_size_[1]);
      DIALS_ASSERT(bbox[2] >= 0 && bbox[3] <= image_size_[0]);

      // Create the index generator for each coordinate of the bounding box
      GridIndexGenerator index(cs, bbox[0], bbox[2], step_size_,
        grid_half_size_, s1_map_);

      // Get the input and output sizes
      std::size_t grid_depth = grid.accessor().all()[0];
      std::size_t grid_height = grid.accessor().all()[1];
      std::size_t grid_width = grid.accessor().all()[2];
      std::size_t image_depth = image.accessor().all()[0];
      std::size_t image_height = image.accessor().all()[1];
      std::size_t image_width = image.accessor().all()[2];

      // Initialise all output counts to zero
      for (std::size_t j = 0; j < image.size(); ++j) {
        image[j] = 0.0;
      }

      // Loop through all the input pixels
      for (std::size_t j = 0; j < image_height; ++j) {
        for (std::size_t i = 0; i < image_width; ++i) {

          // Get the x, y coords of the target point
          vec2<double> ixy00 = index(j, i);
          vec2<double> ixy01 = index(j, i+1);
          vec2<double> ixy10 = index(j+1, i);
          vec2<double> ixy11 = index(j+1, i+1);

          // Get the range of new grid points
          double4 ix(ixy00[0], ixy01[0], ixy10[0], ixy11[0]);
          double4 iy(ixy00[1], ixy01[1], ixy10[1], ixy11[1]);
          int gx0 = (int)floor(min(ix.const_ref()));
          int gy0 = (int)floor(min(iy.const_ref()));
          int gx1 = (int)ceil(max(ix.const_ref()));
          int gy1 = (int)ceil(max(iy.const_ref()));

          // Cap the coordinates within the the output grid
          if (gx0 < 0) gx0 = 0;
          if (gy0 < 0) gy0 = 0;
          if (gx1 > grid_width) gx1 = grid_width;
          if (gy1 > grid_height) gy1 = grid_height;
          if (gx0 >= gx1 || gy0 >= gy1) {
            continue;
          }

          // Create the target polygon and calculate its area
          vert4 target(ixy00, ixy01, ixy11, ixy10);
          double target_area = simple_area(target);
          if (target_area < 0.0) {
            std::swap(target[0], target[3]);
            std::swap(target[1], target[2]);
            target_area = -target_area;
          } else if (target_area == 0.0) {
            continue;
          }

          // Loop over all the pixels within the pixel range
          for (std::size_t jj = gy0; jj < gy1; ++jj) {
            for (std::size_t ii = gx0; ii < gx1; ++ii) {

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
              double xy_fraction = result_area / 1.0;

              // Copy the values to the grid
              for (int kk = 0; kk < grid_depth; ++kk) {
                double value = grid(kk, jj, ii) * xy_fraction;
                for (int k = 0; k < image_depth; ++k) {
                  image(k, j, i) += value * z_fraction(kk, k);
                }
              }
            }
          }
        }
      }
    }

  private:

    flex_vec3_double s1_map_;
    std::size_t grid_half_size_;
    vec2<double> step_size_;
    vec2<std::size_t> image_size_;
  };

}}}} // namespace dials::algorithms::reflection_basis::transform

#endif /* DIALS_ALGORITHMS_REFLEXION_BASIS_MAP_PIXELS_H */
