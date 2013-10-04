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
#include <dials/algorithms/reflection_basis/rebin_pixels.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>

namespace dials { namespace algorithms { namespace reflection_basis {
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
    CoordinateGenerator(const CoordinateSystem &cs, int x0, int y0,
        const af::versa< vec3<double>, af::c_grid<2> > &s1_map)
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
    af::versa< vec3<double>, af::c_grid<2> > s1_map_;
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
        const af::versa< vec3<double>, af::c_grid<2> > &s1_map)
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
    MapPixelsForward(const af::versa< vec3<double>, af::c_grid<2> > &s1_map,
        std::size_t grid_half_size, vec2<double> step_size)
      : s1_map_(s1_map),
        grid_half_size_(grid_half_size),
        step_size_(step_size) {
      DIALS_ASSERT(grid_half_size > 0);
      DIALS_ASSERT(step_size_[0] > 0.0 && step_size_[1] > 0.0);
      image_size_ = s1_map.accessor();
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
    af::versa< double, af::c_grid<3> > operator()(
        const CoordinateSystem &cs, int6 bbox,
        const af::const_ref< double, af::c_grid<3> > &image,
        const af::const_ref< bool, af::c_grid<3> > &mask,
        const af::const_ref< double, af::c_grid<2> > &z_fraction) const {
      std::size_t size = 2 * grid_half_size_ + 1;
      af::c_grid<3> accessor(size, size, size);
      af::versa< double, af::c_grid<3> > grid(accessor, 0);
      af::ref< double, af::c_grid<3> > grid_ref = grid.ref();
      this->operator()(cs, bbox, image, mask, z_fraction, grid_ref);
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
    std::pair<af::versa< double, af::c_grid<3> >,
              af::versa< double, af::c_grid<3> > > operator()(
        const CoordinateSystem &cs, int6 bbox,
        const af::const_ref<double, af::c_grid<3> > &image,
        const af::const_ref<double, af::c_grid<3> > &background,
        const af::const_ref<bool, af::c_grid<3> > &mask,
        const af::const_ref<double, af::c_grid<2> > &z_fraction) const {
      std::size_t size = 2 * grid_half_size_ + 1;
      af::c_grid<3> accessor(size, size, size);
      af::versa<double, af::c_grid<3> > image_grid(accessor, 0);
      af::versa<double, af::c_grid<3> > background_grid(accessor, 0);
      af::ref<double, af::c_grid<3> > image_grid_ref = image_grid.ref();
      af::ref<double, af::c_grid<3> > background_grid_ref = background_grid.ref();
      this->operator()(cs, bbox, image, mask, background,
        z_fraction, image_grid_ref, background_grid_ref);
      return std::make_pair(image_grid, background_grid);
    }

    /**
     * Map the pixels for a reflection
     * @param cs The coordinate system
     * @param bbox The bounding box
     * @param image The image array
     * @param mask The mask array
     * @param z_fraction The z fraction array
     * @param grid The grid array (assumed to be initialised to zero)
     */
    void operator()(const CoordinateSystem &cs, int6 bbox,
        const af::const_ref<double, af::c_grid<3> > &image,
        const af::const_ref<bool, af::c_grid<3> > &mask,
        const af::const_ref<double, af::c_grid<2> > z_fraction,
        af::ref<double, af::c_grid<3> > grid) const {

      // Check array sizes
      DIALS_ASSERT(grid.accessor().all_eq(2 * grid_half_size_ + 1));
      DIALS_ASSERT(z_fraction.accessor()[0] == bbox[5] - bbox[4]);
      DIALS_ASSERT(z_fraction.accessor()[1] == 2 * grid_half_size_ + 1);
      DIALS_ASSERT(image.accessor()[0] == bbox[5] - bbox[4]);
      DIALS_ASSERT(image.accessor()[1] == bbox[3] - bbox[2]);
      DIALS_ASSERT(image.accessor()[2] == bbox[1] - bbox[0]);
      DIALS_ASSERT(mask.accessor().all_eq(image.accessor()));
      DIALS_ASSERT(bbox[0] >= 0 && bbox[1] < image_size_[1]);
      DIALS_ASSERT(bbox[2] >= 0 && bbox[3] < image_size_[0]);

      // Create the index generator for each coordinate of the bounding box
      GridIndexGenerator index(cs, bbox[0], bbox[2], step_size_,
        grid_half_size_, s1_map_);

      // Call the rebinning routine and map just the signal
      vec2<int> isize(image.accessor()[1], image.accessor()[2]);
      vec2<int> osize(grid.accessor()[1], grid.accessor()[2]);
      rebin_pixels_internal(isize, osize, index,
        MapSignal(mask, image, grid, z_fraction));
    }

    /**
     * Map the pixels for a reflection
     * @param cs The coordinate system
     * @param bbox The bounding box
     * @param image The image array
     * @param mask The mask array
     * @param background The background array
     * @param z_fraction The z fraction array
     * @param image_grid The grid array (assumed to be initialised to zero)
     * @param background_grid The background grid (assumed zero)
     */
    void operator()(const CoordinateSystem &cs, int6 bbox,
        const af::const_ref<double, af::c_grid<3> > &image,
        const af::const_ref<bool, af::c_grid<3> > &mask,
        const af::const_ref<double, af::c_grid<3> >&background,
        const af::const_ref<double, af::c_grid<2> > z_fraction,
        af::ref<double, af::c_grid<3> > image_grid,
        af::ref<double, af::c_grid<3> > background_grid) const {

      // Check array sizes
      DIALS_ASSERT(image_grid.accessor().all_eq(2 * grid_half_size_ + 1));
      DIALS_ASSERT(background_grid.accessor().all_eq(image_grid.accessor()));
      DIALS_ASSERT(z_fraction.accessor()[0] == bbox[5] - bbox[4]);
      DIALS_ASSERT(z_fraction.accessor()[1] == 2 * grid_half_size_ + 1);
      DIALS_ASSERT(image.accessor()[0] == bbox[5] - bbox[4]);
      DIALS_ASSERT(image.accessor()[1] == bbox[3] - bbox[2]);
      DIALS_ASSERT(image.accessor()[2] == bbox[1] - bbox[0]);
      DIALS_ASSERT(mask.accessor().all_eq(image.accessor()));
      DIALS_ASSERT(background.accessor().all_eq(image.accessor()));
      DIALS_ASSERT(bbox[0] >= 0 && bbox[1] < image_size_[1]);
      DIALS_ASSERT(bbox[2] >= 0 && bbox[3] < image_size_[0]);

      // Create the index generator for each coordinate of the bounding box
      GridIndexGenerator index(cs, bbox[0], bbox[2], step_size_,
        grid_half_size_, s1_map_);

      // Call the rebinning routine and map just the signal
      vec2<int> isize(image.accessor()[1], image.accessor()[2]);
      vec2<int> osize(image_grid.accessor()[1], image_grid.accessor()[2]);
      rebin_pixels_internal(isize, osize, index,
        MapSignalAndBackground(mask, image, background, image_grid,
          background_grid, z_fraction));
    }

  private:

    /**
     * A struct used as a callback in the function rebin pixels
     */
    struct MapSignal {
      MapSignal(af::const_ref<bool, af::c_grid<3> > mask_,
                af::const_ref<double, af::c_grid<3> > input_,
                af::ref<double, af::c_grid<3> > output_,
                af::const_ref<double, af::c_grid<2> > z_fraction_)
        : mask(mask_),
          input(input_),
          output(output_),
          z_fraction(z_fraction_),
          input_depth(input.accessor()[0]),
          output_depth(output.accessor()[0]){}

      void operator()(std::size_t j, std::size_t i,
                      std::size_t jj, std::size_t ii,
                      double xy_fraction) {
        // Copy the values to the grid
        for (int k = 0; k < input_depth; ++k) {
          if (mask(k, j, i) != 0) {
            double value = input(k, j, i) * xy_fraction;
            for (int kk = 0; kk < output_depth; ++kk) {
              output(kk, jj, ii) += value * z_fraction(k, kk);
            }
          }
        }
      }

      af::const_ref<bool, af::c_grid<3> > mask;
      af::const_ref<double, af::c_grid<3> > input;
      af::ref<double, af::c_grid<3> > output;
      af::const_ref<double, af::c_grid<2> > z_fraction;
      std::size_t input_depth;
      std::size_t output_depth;
    };

    /**
     * A struct used as a callback in the function rebin pixels
     */
    struct MapSignalAndBackground {
      MapSignalAndBackground(
            const af::const_ref<bool, af::c_grid<3> > &mask_,
            const af::const_ref<double, af::c_grid<3> > &image_,
            const af::const_ref<double, af::c_grid<3> > &background_,
            af::ref<double, af::c_grid<3> > image_grid_,
            af::ref<double, af::c_grid<3> > background_grid_,
            const af::const_ref<double, af::c_grid<2> > &z_fraction_)
        : mask(mask_),
          image(image_),
          background(background_),
          image_grid(image_grid_),
          background_grid(background_grid_),
          z_fraction(z_fraction_),
          image_depth(image.accessor()[0]),
          grid_depth(image_grid.accessor()[0]){}

      void operator()(std::size_t j, std::size_t i,
                      std::size_t jj, std::size_t ii,
                      double xy_fraction) {
        // Copy the values to the grid
        for (int k = 0; k < image_depth; ++k) {
          if (mask(k, j, i) != 0) {
            double ivalue = image(k, j, i) * xy_fraction;
            double bvalue = background(k, j, i) * xy_fraction;
            for (int kk = 0; kk < grid_depth; ++kk) {
              image_grid(kk, jj, ii) += ivalue * z_fraction(k, kk);
              background_grid(kk, jj, ii) += bvalue * z_fraction(k, kk);
            }
          }
        }
      }

      af::const_ref<bool, af::c_grid<3> > mask;
      af::const_ref<double, af::c_grid<3> > image;
      af::const_ref<double, af::c_grid<3> > background;
      af::ref<double, af::c_grid<3> > image_grid;
      af::ref<double, af::c_grid<3> > background_grid;
      af::const_ref<double, af::c_grid<2> > z_fraction;
      std::size_t image_depth;
      std::size_t grid_depth;
    };

    af::versa< vec3<double>, af::c_grid<2> > s1_map_;
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
    MapPixelsReverse(const af::versa< vec3<double>, af::c_grid<2> > &s1_map,
        std::size_t grid_half_size, vec2<double> step_size)
      : s1_map_(s1_map),
        grid_half_size_(grid_half_size),
        step_size_(step_size) {
      DIALS_ASSERT(grid_half_size > 0);
      DIALS_ASSERT(step_size_[0] > 0.0 && step_size_[1] > 0.0);
      image_size_ = s1_map.accessor();
    }

    /**
     * Map the pixels for a reflection
     * @param cs The coordinate system
     * @param bbox The bounding box
     * @param grid The grid
     * @param z_fraction The z fraction array
     * @returns The mapped image
     */
    af::versa<double, af::c_grid<3> > operator()(
        const CoordinateSystem &cs, int6 bbox,
        const af::const_ref<double, af::c_grid<3> > &grid,
        const af::const_ref<double, af::c_grid<2> > &z_fraction) const {
      af::versa<double, af::c_grid<3> > image(
        af::c_grid<3>(bbox[5] - bbox[4],
                      bbox[3] - bbox[2],
                      bbox[1] - bbox[0]), 0);
      af::ref<double, af::c_grid<3> > image_ref = image.ref();
      this->operator()(cs, bbox, grid, z_fraction, image_ref);
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
        const af::const_ref< double, af::c_grid<3> > &grid,
        const af::const_ref< double, af::c_grid<2> > &z_fraction,
        af::ref< double, af::c_grid<3> > image) const {

      // Check array sizes
      DIALS_ASSERT(grid.accessor().all_eq(2 * grid_half_size_ + 1));
      DIALS_ASSERT(z_fraction.accessor()[1] == bbox[5] - bbox[4]);
      DIALS_ASSERT(z_fraction.accessor()[0] == 2 * grid_half_size_ + 1);
      DIALS_ASSERT(image.accessor()[0] == bbox[5] - bbox[4]);
      DIALS_ASSERT(image.accessor()[1] == bbox[3] - bbox[2]);
      DIALS_ASSERT(image.accessor()[2] == bbox[1] - bbox[0]);
      DIALS_ASSERT(bbox[0] >= 0 && bbox[1] <= image_size_[1]);
      DIALS_ASSERT(bbox[2] >= 0 && bbox[3] <= image_size_[0]);

      // Create the index generator for each coordinate of the bounding box
      GridIndexGenerator index(cs, bbox[0], bbox[2], step_size_,
        grid_half_size_, s1_map_);

//      FIXME: To use templated stuff need to divie by subject area and not
//      by target area. Needs a bit of rejigging.
//      // Call the rebinning routine and map the grid
//      vec2<int> isize(image.accessor()[1], image.accessor()[2]);
//      vec2<int> osize(grid.accessor()[1], grid.accessor()[2]);
//      rebin_pixels_internal(isize, osize, index,
//        MapGrid(image, grid, z_fraction));

      // Get the input and output sizes
      std::size_t grid_depth = grid.accessor()[0];
      std::size_t grid_height = grid.accessor()[1];
      std::size_t grid_width = grid.accessor()[2];
      std::size_t image_depth = image.accessor()[0];
      std::size_t image_height = image.accessor()[1];
      std::size_t image_width = image.accessor()[2];

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
              // Then redistribute the values from the target grid to subject.
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

    /**
     * A struct used as a callback in the function rebin pixels
     */
    struct MapGrid {
      MapGrid(af::ref< double, af::c_grid<3> > image_,
              af::const_ref< double, af::c_grid<3> > &grid_,
              af::const_ref< double, af::c_grid<2> > &z_fraction_)
        : image(image_),
          grid(grid_),
          z_fraction(z_fraction_),
          image_depth(image.accessor()[0]),
          grid_depth(grid.accessor()[0]){}

      void operator()(std::size_t j, std::size_t i,
                      std::size_t jj, std::size_t ii,
                      double xy_fraction) {
        // Copy the values to the grid
        for (int kk = 0; kk < grid_depth; ++kk) {
          double value = grid(kk, jj, ii) * xy_fraction;
          for (int k = 0; k < image_depth; ++k) {
            image(k, j, i) += value * z_fraction(kk, k);
          }
        }
      }

      af::ref< double, af::c_grid<3> > image;
      af::const_ref< double, af::c_grid<3> > grid;
      af::const_ref< double, af::c_grid<2> > z_fraction;
      std::size_t image_depth;
      std::size_t grid_depth;
    };

    af::versa< vec3<double>, af::c_grid<2> > s1_map_;
    std::size_t grid_half_size_;
    vec2<double> step_size_;
    vec2<std::size_t> image_size_;
  };

}}}} // namespace dials::algorithms::reflection_basis::transform

#endif /* DIALS_ALGORITHMS_REFLEXION_BASIS_MAP_PIXELS_H */
