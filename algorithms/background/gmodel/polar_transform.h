/*
 * polar_transform.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_BACKGROUND_GMODEL_POLAR_TRANSFORM_H
#define DIALS_ALGORITHMS_BACKGROUND_GMODEL_POLAR_TRANSFORM_H

#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dials/algorithms/polygon/spatial_interpolation.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using dials::algorithms::polygon::clip::vert4;
  using dials::algorithms::polygon::spatial_interpolation::grid_to_quad;
  using dials::algorithms::polygon::spatial_interpolation::Match;
  using dials::algorithms::polygon::spatial_interpolation::quad_to_grid;
  using dxtbx::model::BeamBase;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Panel;

  namespace detail {
    template <typename T>
    int sign(T val) {
      return (T(0) < val) - (val < T(0));
    }
  }  // namespace detail

  class PolarTransformResult {
  public:
    PolarTransformResult(af::versa<double, af::c_grid<2> > data,
                         af::versa<bool, af::c_grid<2> > mask)
        : data_(data), mask_(mask) {}

    af::versa<double, af::c_grid<2> > data() const {
      return data_;
    }

    af::versa<bool, af::c_grid<2> > mask() const {
      return mask_;
    }

  protected:
    af::versa<double, af::c_grid<2> > data_;
    af::versa<bool, af::c_grid<2> > mask_;
  };

  /**
   * A class to do a polar transform along resolution
   */
  class PolarTransform {
  public:
    /**
     * Initialise the class
     * @param beam The beam model
     * @param panel The panel model
     * @param goniometer The goniometer model
     */
    PolarTransform(const BeamBase &beam,
                   const Panel &panel,
                   const Goniometer &goniometer) {
      // Set some image sizes
      vec2<std::size_t> image_size = panel.get_image_size();
      DIALS_ASSERT(image_size[0] > 0);
      DIALS_ASSERT(image_size[1] > 0);
      image_grid_ = af::c_grid<2>(image_size[1], image_size[0]);

      // Allocate map arrays
      image_xmap_ = af::versa<double, af::c_grid<2> >(
        af::c_grid<2>(image_size[1] + 1, image_size[0] + 1));
      image_ymap_ = af::versa<double, af::c_grid<2> >(
        af::c_grid<2>(image_size[1] + 1, image_size[0] + 1));
      discontinuity_ = af::versa<bool, af::c_grid<2> >(
        af::c_grid<2>(image_size[1] + 1, image_size[0] + 1));

      // Setup x, y, z axis for transform
      vec3<double> s0 = beam.get_s0().normalize();
      vec3<double> m2 = goniometer.get_rotation_axis().normalize();
      vec3<double> zaxis = s0;
      vec3<double> yaxis = zaxis.cross(m2);
      vec3<double> xaxis = zaxis.cross(yaxis);

      // Generate temporary array of x, y z coords
      af::versa<double, af::c_grid<2> > temp_x(
        af::c_grid<2>(image_size[1] + 1, image_size[0] + 1));
      af::versa<double, af::c_grid<2> > temp_y(
        af::c_grid<2>(image_size[1] + 1, image_size[0] + 1));
      af::versa<double, af::c_grid<2> > temp_z(
        af::c_grid<2>(image_size[1] + 1, image_size[0] + 1));
      for (std::size_t j = 0; j < image_size[1] + 1; ++j) {
        for (std::size_t i = 0; i < image_size[0] + 1; ++i) {
          vec3<double> s1 = panel.get_pixel_lab_coord(vec2<double>(i, j)).normalize();
          double z = s1 * zaxis;
          double y = s1 * yaxis;
          double x = s1 * xaxis;
          temp_x(j, i) = x;
          temp_y(j, i) = y;
          temp_z(j, i) = z;
        }
      }

      // Generate polar coords from image pixels
      for (std::size_t j = 0; j < image_size[1] + 1; ++j) {
        for (std::size_t i = 0; i < image_size[0] + 1; ++i) {
          double z = temp_z(j, i);
          double y = temp_y(j, i);
          double x = temp_x(j, i);
          image_xmap_(j, i) = std::acos(z);
          image_ymap_(j, i) = std::atan2(y, x);

          // Create an image with the discontinuity shown
          discontinuity_(j, i) = false;
          if (x <= 0) {
            if (j < image_size[1]) {
              int s1 = detail::sign(temp_y(j, i));
              int s2 = detail::sign(temp_y(j + 1, i));
              int s3 = detail::sign(temp_y(j, i + 1));
              int s4 = detail::sign(temp_y(j + 1, i + 1));
              if (s1 != s2 || s1 != s3 || s1 != s4 || s2 != s3 || s2 != s4
                  || s3 != s4) {
                discontinuity_(j, i) = true;
              }
            }
          }
        }
      }

      // Get the min/max x and y
      double map_xmin = image_xmap_[0];
      double map_xmax = image_xmap_[0];
      double map_ymin = image_ymap_[1];
      double map_ymax = image_ymap_[1];
      std::size_t map_xmin_ind = 0;
      for (std::size_t i = 1; i < image_xmap_.size(); ++i) {
        if (image_xmap_[i] < map_xmin) map_xmin_ind = i;
        if (image_xmap_[i] < map_xmin) map_xmin = image_xmap_[i];
        if (image_xmap_[i] > map_xmax) map_xmax = image_xmap_[i];
        if (image_ymap_[i] < map_ymin) map_ymin = image_ymap_[i];
        if (image_ymap_[i] > map_ymax) map_ymax = image_ymap_[i];
      }

      // Compute the polar grid x size
      int xind1 = (map_xmin_ind % (image_size[0] + 1));
      int yind1 = (map_xmin_ind / (image_size[0] + 1));
      int xsize = (int)image_size[0];
      int ysize = (int)image_size[1];
      double d1 =
        std::sqrt((double)((xind1 - 0) * (xind1 - 0) + (yind1 - 0) * (yind1 - 0)));
      double d2 = std::sqrt(
        (double)((xind1 - 0) * (xind1 - 0) + (yind1 - ysize) * (yind1 - ysize)));
      double d3 = std::sqrt(
        (double)((xind1 - xsize) * (xind1 - xsize) + (yind1 - 0) * (yind1 - 0)));
      double d4 = std::sqrt((double)((xind1 - xsize) * (xind1 - xsize)
                                     + (yind1 - ysize) * (yind1 - ysize)));
      std::size_t polar_xsize =
        (std::size_t)std::ceil(std::max(std::max(d1, d2), std::max(d3, d4)));
      std::size_t polar_ysize =
        (std::size_t)std::ceil((double)(image_size[0] * image_size[1]) / polar_xsize);
      polar_grid_ = af::c_grid<2>(polar_ysize, polar_xsize);

      // The polar grid limits
      double polar_xmin = map_xmin;
      double polar_xmax = map_xmax;
      double polar_xstep = (polar_xmax - polar_xmin) / polar_xsize;
      double polar_ymin = map_ymin;
      double polar_ymax = map_ymax;
      double polar_ystep = (polar_ymax - polar_ymin) / polar_ysize;

      // Convert polar coords to indices
      for (std::size_t j = 0; j < image_size[1] + 1; ++j) {
        for (std::size_t i = 0; i < image_size[0] + 1; ++i) {
          double px = image_xmap_(j, i);
          double py = image_ymap_(j, i);
          double pi = ((px - polar_xmin) / polar_xstep);
          double pj = ((py - polar_ymin) / polar_ystep);
          image_xmap_(j, i) = pi;
          image_ymap_(j, i) = pj;
        }
      }
    }

    /**
     * @returns The x map coordinates
     */
    af::versa<double, af::c_grid<2> > image_xmap() const {
      return image_xmap_;
    }

    /**
     * @returns The y map coordinates
     */
    af::versa<double, af::c_grid<2> > image_ymap() const {
      return image_ymap_;
    }

    /**
     * @returns The discontinuity coordinates
     */
    af::versa<bool, af::c_grid<2> > discontinuity() const {
      return discontinuity_;
    }

    /**
     * @param j The y index
     * @param i The x index
     * @returns The grid coordinate
     */
    vec2<double> gc(std::size_t j, std::size_t i) const {
      DIALS_ASSERT(image_xmap_.accessor().all_eq(image_ymap_.accessor()));
      DIALS_ASSERT(image_xmap_.accessor()[0] == image_grid_[0] + 1);
      DIALS_ASSERT(image_xmap_.accessor()[1] == image_grid_[1] + 1);
      DIALS_ASSERT(j <= image_grid_[0]);
      DIALS_ASSERT(i <= image_grid_[1]);
      return vec2<double>(image_xmap_(j, i), image_ymap_(j, i));
    }

    /**
     * transform to polar
     * @param data The image data
     * @param mask The image mask
     * @returns The transformed data
     */
    PolarTransformResult to_polar(
      const af::const_ref<double, af::c_grid<2> > &data,
      const af::const_ref<bool, af::c_grid<2> > &mask) const {
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(data.accessor().all_eq(image_grid_));
      DIALS_ASSERT(data.accessor()[0] + 1 == discontinuity_.accessor()[0]);
      DIALS_ASSERT(data.accessor()[1] + 1 == discontinuity_.accessor()[1]);
      af::versa<double, af::c_grid<2> > data_out(polar_grid_, 0);
      af::versa<bool, af::c_grid<2> > mask_out(polar_grid_, true);
      af::versa<bool, af::c_grid<2> > mask_tmp(polar_grid_, false);

      for (std::size_t j = 0; j < image_grid_[0]; ++j) {
        for (std::size_t i = 0; i < image_grid_[1]; ++i) {
          // FIXME - Discarding pixels where the angle wraps round. Need to
          // handle this better

          if (discontinuity_(j, i)) {
            continue;
          }
          vert4 input(gc(j, i), gc(j, i + 1), gc(j + 1, i + 1), gc(j + 1, i));
          af::shared<Match> matches = quad_to_grid(input, polar_grid_, 0);
          for (int m = 0; m < matches.size(); ++m) {
            double fraction = matches[m].fraction;
            int index = matches[m].out;
            int ii = index % polar_grid_[1];
            int jj = index / polar_grid_[1];
            DIALS_ASSERT(jj >= 0 && jj < polar_grid_[0]);
            DIALS_ASSERT(ii >= 0 && ii < polar_grid_[1]);
            if (mask_out(jj, ii) && mask(j, i)) {
              data_out(jj, ii) += data(j, i) * fraction;
              mask_tmp(jj, ii) = true;
            } else {
              mask_out(jj, ii) = false;
            }
          }
        }
      }

      // Apply both masks
      for (std::size_t i = 0; i < mask_out.size(); ++i) {
        if (mask_out[i] && mask_tmp[i]) {
          mask_out[i] = true;
        } else {
          mask_out[i] = false;
        }
      }

      return PolarTransformResult(data_out, mask_out);
    }

    /**
     * transform from polar
     * @param data The polar data
     * @param mask The polar mask
     * @returns The transformed data
     */
    PolarTransformResult from_polar(
      const af::const_ref<double, af::c_grid<2> > &data,
      const af::const_ref<bool, af::c_grid<2> > &mask) const {
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(data.accessor().all_eq(polar_grid_));
      af::versa<double, af::c_grid<2> > data_out(image_grid_, 0);
      af::versa<bool, af::c_grid<2> > mask_out(image_grid_, true);
      for (std::size_t j = 0; j < image_grid_[0]; ++j) {
        for (std::size_t i = 0; i < image_grid_[1]; ++i) {
          // FIXME - Discarding pixels where the angle wraps round. Need to
          // handle this better

          if (discontinuity_(j, i)) {
            mask_out(j, i) = false;
            continue;
          }
          vert4 input(gc(j, i), gc(j, i + 1), gc(j + 1, i + 1), gc(j + 1, i));
          af::shared<Match> matches = grid_to_quad(input, polar_grid_, 0);
          for (int m = 0; m < matches.size(); ++m) {
            double fraction = matches[m].fraction;
            int index = matches[m].in;
            int ii = index % polar_grid_[1];
            int jj = index / polar_grid_[1];
            DIALS_ASSERT(ii < polar_grid_[1]);
            DIALS_ASSERT(jj < polar_grid_[0]);
            DIALS_ASSERT(ii >= 0);
            DIALS_ASSERT(jj >= 0);
            if (mask_out(j, i) && mask(jj, ii)) {
              data_out(j, i) += data(jj, ii) * fraction;
            } else {
              mask_out(j, i) = false;
            }
          }
        }
      }

      return PolarTransformResult(data_out, mask_out);
    }

  protected:
    af::c_grid<2> image_grid_;
    af::c_grid<2> polar_grid_;
    af::versa<double, af::c_grid<2> > image_xmap_;
    af::versa<double, af::c_grid<2> > image_ymap_;
    af::versa<bool, af::c_grid<2> > discontinuity_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_BACKGROUND_GMODEL_POLAR_TRANSFORM_H
