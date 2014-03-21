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
#include <scitbx/vec2.h>
#include <cmath>
#include <scitbx/array_family/flex_types.h>

#include <cstdlib>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/array_family/accessors/c_grid_periodic.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/spot_prediction/rotation_angles.h>
#include <dxtbx/model/scan_helpers.h>


namespace dials { namespace algorithms {

  using dxtbx::model::is_angle_in_range;

  // helper function for sampling_volume_map
  bool are_angles_in_range(
    af::ref<vec2<double> > const & angle_ranges,
    vec2<double> const & angles) {
      for (std::size_t i=0; i<2; i++) {
        double angle = angles[i];
        for (std::size_t j=0; j<angle_ranges.size(); j++) {
          if (is_angle_in_range(angle_ranges[j], angle)) {
            return true;
          }
        }
      }
    return false;
  }

  // compute a map of the sampling volume of a scan
  void sampling_volume_map(
    af::ref<double, af::c_grid<3> > const & data,
    af::ref<vec2<double> > const & angle_ranges,
    vec3 <double> s0, vec3 <double> m2,
    double const & rl_grid_spacing,
    double d_min,
    double b_iso)
  {
    typedef af::c_grid<3>::index_type index_t;
    index_t const gridding_n_real = index_t(data.accessor());

    RotationAngles calculate_rotation_angles_(s0, m2);

    double one_over_d_sq_min = 1/(d_min*d_min);

    for (std::size_t i=0; i<gridding_n_real[0]; i++) {
      double i_rl = (double(i) - double(gridding_n_real[0]/2.0)) * rl_grid_spacing;
      double i_rl_sq = i_rl * i_rl;
      for (std::size_t j=0; j<gridding_n_real[1]; j++) {
        double j_rl = (double(j) - double(gridding_n_real[1]/2.0)) * rl_grid_spacing;
        double j_rl_sq = j_rl * j_rl;
        for (std::size_t k=0; k<gridding_n_real[2]; k++) {
          double k_rl = (double(k) - double(gridding_n_real[2]/2.0)) * rl_grid_spacing;
          double k_rl_sq = k_rl * k_rl;
          double reciprocal_length_sq = (i_rl_sq + j_rl_sq + k_rl_sq);
          if (reciprocal_length_sq > one_over_d_sq_min) {
            continue;
          }
          vec3 <double> pstar0(i_rl, j_rl, k_rl);

          // Try to calculate the diffracting rotation angles
          vec2 <double> phi;
          try {
            phi = calculate_rotation_angles_(pstar0);
          } catch(error) {
            continue;
          }

          // Check that the angles are within the rotation range
          if (are_angles_in_range(angle_ranges, phi)) {
            double T;
            if (b_iso != 0) {
              T = std::exp(-200 * reciprocal_length_sq / 4);
            }
            else {
              T = 1;
            }
            data(i,j,k) = T;
          }
        }
      }
    }
  }


  /*
  Peak-finding algorithm inspired by the CLEAN algorithm of
  Högbom, J. A. 1974, A&AS, 15, 417.

  See also:
    http://dx.doi.org/10.1051/0004-6361/200912148
  */
  af::shared<vec3<int> > clean_3d(
    af::const_ref<double, af::c_grid_periodic<3> > const & dirty_beam,
    af::ref<double, af::c_grid<3> > const & dirty_map,
    std::size_t n_peaks,
    double gamma=1)
  {
    af::shared<vec3<int> > peaks;
    typedef af::c_grid<3>::index_type index_t;
    index_t const gridding_n_real = index_t(dirty_map.accessor());
    DIALS_ASSERT(dirty_map.size() == dirty_beam.size());
    double max_db = af::max(dirty_beam);
    for (std::size_t i_peak=0; i_peak<n_peaks; i_peak++) {
      vec3<int> max_index = vec3<int>(0,0,0);
      double max_value = dirty_map(max_index);
      for (int i=0; i<gridding_n_real[0]; i++){
        for (int j=0; j<gridding_n_real[1]; j++){
          for (int k=0; k<gridding_n_real[2]; k++){
            double v = dirty_map(i,j,k);
            if (v > max_value) {
              max_index = vec3<int>(i,j,k);
              max_value = v;
            }
          }
        }
      }
      vec3<int> shift = max_index;
      peaks.push_back(shift);
      double scale = max_value/max_db * gamma;
      for (int i=0; i<gridding_n_real[0]; i++){
        for (int j=0; j<gridding_n_real[1]; j++){
          for (int k=0; k<gridding_n_real[2]; k++){
            dirty_map(i,j,k) -= dirty_beam(i-shift[0], j-shift[1], k-shift[2]) * scale;
          }
        }
      }
    }
    return peaks;
  }


}}

#endif
