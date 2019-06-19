/*
 * fft3d.h
 *
 *  Copyright (C) 2014 Diamond Light Source
 *
 *  Author: Richard Gildea
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_FFT3D_H
#define DIALS_ALGORITHMS_INTEGRATION_FFT3D_H
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <scitbx/vec2.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/math/utils.h>

#include <cstdlib>
#include <scitbx/array_family/versa_matrix.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/spot_prediction/rotation_angles.h>
#include <dxtbx/model/scan_helpers.h>

namespace dials { namespace algorithms {

  using dxtbx::model::is_angle_in_range;

  // helper function for sampling_volume_map
  bool are_angles_in_range(af::ref<vec2<double> > const& angle_ranges,
                           vec2<double> const& angles) {
    for (std::size_t i = 0; i < 2; i++) {
      double angle = angles[i];
      for (std::size_t j = 0; j < angle_ranges.size(); j++) {
        if (is_angle_in_range(angle_ranges[j], angle)) {
          return true;
        }
      }
    }
    return false;
  }

  // compute a map of the sampling volume of a scan
  void sampling_volume_map(af::ref<double, af::c_grid<3> > const& data,
                           af::ref<vec2<double> > const& angle_ranges,
                           vec3<double> s0,
                           vec3<double> m2,
                           double const& rl_grid_spacing,
                           double d_min,
                           double b_iso) {
    typedef af::c_grid<3>::index_type index_t;
    index_t const gridding_n_real = index_t(data.accessor());

    RotationAngles calculate_rotation_angles_(s0, m2);

    double one_over_d_sq_min = 1 / (d_min * d_min);

    for (std::size_t i = 0; i < gridding_n_real[0]; i++) {
      double i_rl = (double(i) - double(gridding_n_real[0] / 2.0)) * rl_grid_spacing;
      double i_rl_sq = i_rl * i_rl;
      for (std::size_t j = 0; j < gridding_n_real[1]; j++) {
        double j_rl = (double(j) - double(gridding_n_real[1] / 2.0)) * rl_grid_spacing;
        double j_rl_sq = j_rl * j_rl;
        for (std::size_t k = 0; k < gridding_n_real[2]; k++) {
          double k_rl =
            (double(k) - double(gridding_n_real[2] / 2.0)) * rl_grid_spacing;
          double k_rl_sq = k_rl * k_rl;
          double reciprocal_length_sq = (i_rl_sq + j_rl_sq + k_rl_sq);
          if (reciprocal_length_sq > one_over_d_sq_min) {
            continue;
          }
          vec3<double> pstar0(i_rl, j_rl, k_rl);

          // Try to calculate the diffracting rotation angles
          vec2<double> phi;
          try {
            phi = calculate_rotation_angles_(pstar0);
          } catch (error) {
            continue;
          }

          // Check that the angles are within the rotation range
          if (are_angles_in_range(angle_ranges, phi)) {
            double T;
            if (b_iso != 0) {
              T = std::exp(-b_iso * reciprocal_length_sq / 4);
            } else {
              T = 1;
            }
            data(i, j, k) = T;
          }
        }
      }
    }
  }

  /*
  Peak-finding algorithm inspired by the CLEAN algorithm of
  Högbom, J. A. 1974, A&AS, 15, 417.

  See also:
    https://doi.org/10.1051/0004-6361/200912148
  */
  af::shared<vec3<int> > clean_3d(
    af::const_ref<double, af::c_grid<3> > const& dirty_beam,
    af::ref<double, af::c_grid<3> > const& dirty_map,
    std::size_t n_peaks,
    double gamma = 1) {
    af::shared<vec3<int> > peaks;
    typedef af::c_grid<3>::index_type index_t;
    index_t const gridding_n_real = index_t(dirty_map.accessor());
    DIALS_ASSERT(dirty_map.size() == dirty_beam.size());
    double max_db = af::max(dirty_beam);
    af::c_grid<3> accessor(dirty_map.accessor());  // index_type conversion

    const int height = int(gridding_n_real[0]);
    const int depth = int(gridding_n_real[1]);
    const int width = int(gridding_n_real[2]);
    const long height_depth = height * depth;

    int max_idx = af::max_index(dirty_map);

    for (std::size_t i_peak = 0; i_peak < n_peaks; i_peak++) {
      // Find the maximum value in the map - this is the next "peak"
      const index_t shift = accessor.index_nd(max_idx);
      peaks.push_back(vec3<int>(shift));

      // reposition the dirty beam on the current peak and subtract from
      // the dirty map
      const double max_value = dirty_map[max_idx];
      const double scale = max_value / max_db * gamma;
      max_idx = 0;  // reset for next cycle
#pragma omp parallel for
      for (int i = 0; i < width; i++) {
        int i_db = i - shift[0];
        if (i_db < 0) {
          i_db += width;
        } else if (i_db >= width) {
          i_db -= width;
        }
        // DIALS_ASSERT(i_db >= 0 && i_db < width);
        const long ipart_dm = i * height_depth;
        const long ipart_db = i_db * height_depth;
        for (int j = 0; j < height; j++) {
          int j_db = j - shift[1];
          if (j_db < 0) {
            j_db += height;
          } else if (j_db >= height) {
            j_db -= height;
          }
          // DIALS_ASSERT(j_db >= 0 && j_db < height);
          const long ijpart_dm = ipart_dm + j * depth;
          const long ijpart_db = ipart_db + j_db * depth;
          for (int k = 0; k < depth; k++) {
            int k_db = k - shift[2];
            if (k_db < 0) {
              k_db += depth;
            } else if (k_db >= depth) {
              k_db -= depth;
            }
            // DIALS_ASSERT(k_db >= 0 && k_db < depth);
            const long idx_dm = ijpart_dm + k;
            const long idx_db = ijpart_db + k_db;
            dirty_map[idx_dm] -= dirty_beam[idx_db] * scale;
            if (dirty_map[max_idx] < dirty_map[idx_dm])
#pragma omp critical(max_idx)
            {
              max_idx = idx_dm;
            }
          }
        }
      }
    }
    return peaks;
  }

  void map_centroids_to_reciprocal_space_grid(
    af::ref<double, af::c_grid<3> > const& grid,
    af::const_ref<vec3<double> > const& reciprocal_space_vectors,
    af::ref<bool> const& selection,
    double d_min,
    double b_iso = 0) {
    typedef af::c_grid<3>::index_type index_t;
    index_t const gridding_n_real = index_t(grid.accessor());
    DIALS_ASSERT(d_min >= 0);
    DIALS_ASSERT(gridding_n_real[0] == gridding_n_real[1]);
    DIALS_ASSERT(gridding_n_real[0] == gridding_n_real[2]);

    const int n_points = gridding_n_real[0];
    const double rlgrid = 2 / (d_min * n_points);
    const double one_over_rlgrid = 1 / rlgrid;
    const int half_n_points = n_points / 2;

    for (int i = 0; i < reciprocal_space_vectors.size(); i++) {
      if (!selection[i]) {
        continue;
      }
      const vec3<double> v = reciprocal_space_vectors[i];
      const double v_length = v.length();
      const double d_spacing = 1 / v_length;
      if (d_spacing < d_min) {
        selection[i] = false;
        continue;
      }
      vec3<int> coord;
      for (int j = 0; j < 3; j++) {
        coord[j] = scitbx::math::iround(v[j] * one_over_rlgrid) + half_n_points;
      }
      if ((coord.max() >= n_points) || coord.min() < 0) {
        selection[i] = false;
        continue;
      }
      double T;
      if (b_iso != 0) {
        T = std::exp(-b_iso * v_length * v_length / 4.0);
      } else {
        T = 1;
      }
      grid(coord) = T;
    }
  }

}}  // namespace dials::algorithms

#endif
