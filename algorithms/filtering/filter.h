/*
 * filter.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_FILTER_H
#define DIALS_ALGORITHMS_FILTER_H

#include <cmath>
#include <scitbx/vec3.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/beam.h>
#include <dials/algorithms/image/threshold/unimodal.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace algorithms { namespace filter {

  using scitbx::vec3;
  using scitbx::af::small;
  using scitbx::af::int6;
  using scitbx::af::flex_int;
  using scitbx::af::flex_bool;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Beam;
  using dials::algorithms::reflection_basis::CoordinateSystem;
  using dials::algorithms::reflection_basis::zeta_factor;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /**
   * Calculate the zeta factor and check its absolute value is above the
   * minimum specified value.
   * @param m2 The rotation axis (normalized)
   * @param s0 The incident beam vector
   * @param s1 The diffracted beam vector
   * @param zeta_min The minimum allowed zeta value
   * @returns True/False, zeta is valid
   */
  inline
  bool is_zeta_valid(vec3<double> m2, vec3<double> s0, vec3<double> s1,
      double zeta_min) {
    return std::abs(zeta_factor(m2, s0, s1)) >= zeta_min;
  }

  /**
   * Calculate the zeta factor and check its absolute value is above the
   * minimum specified value.
   * @param cs The local reflection coordinate system
   * @param zeta_min The minimum allowed zeta value
   * @returns True/False, zeta is valid
   */
  inline
  bool is_zeta_valid(const CoordinateSystem &cs, double zeta_min) {
    return std::abs(cs.zeta()) >= zeta_min;
  }

  /**
   * Calculate the zeta factor and check its absolute value is above the
   * minimum specified value.
   * @param m2 The rotation axis (normalized)
   * @param s0 The incident beam vector
   * @param r The reflection
   * @param zeta_min The minimum allowed zeta value
   * @returns True/False, zeta is valid
   */
  inline
  bool is_zeta_valid(vec3<double> m2, vec3<double> s0, const Reflection &r,
      double zeta_min) {
    return is_zeta_valid(m2, s0, r.get_beam_vector(), zeta_min);
  }

  /**
   * Calculate the zeta factor and check its absolute value is above the
   * minimum specified value.
   * @param g The goniometer
   * @param b The beam
   * @param r The reflection
   * @param zeta_min The minimum allowed zeta value
   * @returns True/False, zeta is valid
   */
  inline
  bool is_zeta_valid(const Goniometer &g, const Beam &b, const Reflection &r,
      double zeta_min) {
    return is_zeta_valid(g.get_rotation_axis(), b.get_s0(), r, zeta_min);
  }

  /**
   * Check if the XDS small angle approximation holds for the local
   * reflection transform.
   *
   * Check that the following condition holds:
   *  (m2.e1)^2 + 2*c3*(m2.e3)*(m2.p*) - c3^2 >= 0
   *
   * @param m2 The rotation axis
   * @param s0 The incident beam vector
   * @param s1 The diffracted beam vector
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the small angle approximation is valid
   */
  inline
  bool is_xds_small_angle_valid(vec3<double> m2, vec3<double> s0,
      vec3<double> s1, double delta_m) {
    vec3<double> ps = (s1 - s0).normalize();
    vec3<double> e1 = s1.cross(s0).normalize();
    vec3<double> e3 = (s1 + s0).normalize();
    double m2e1 = m2 * e1;
    double m2e3 = m2 * e3;
    double m2ps = m2 * ps;
    double c3 = -std::abs(delta_m);
    return (m2e1*m2e1 + 2.0*c3*m2e3*m2ps - c3*c3) >= 0.0;
  }

  /**
   * Check if the XDS small angle approximation holds for the local
   * reflection transform.
   *
   * Check that the following condition holds:
   *  (m2.e1)^2 + 2*c3*(m2.e3)*(m2.p*) - c3^2 >= 0
   *
   * @param cs The local reflection coordinate system
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the small angle approximation is valid
   */
  inline
  bool is_xds_small_angle_valid(const CoordinateSystem &cs, double delta_m) {
    vec3<double> m2 = cs.m2();
    vec3<double> ps = cs.p_star().normalize();
    vec3<double> e1 = cs.e1_axis();
    vec3<double> e3 = cs.e3_axis();
    double m2e1 = m2 * e1;
    double m2e3 = m2 * e3;
    double m2ps = m2 * ps;
    double c3 = -std::abs(delta_m);
    return (m2e1*m2e1 + 2.0*c3*m2e3*m2ps - c3*c3) >= 0.0;
  }

  /**
   * Check if the XDS small angle approximation holds for the local
   * reflection transform.
   *
   * Check that the following condition holds:
   *  (m2.e1)^2 + 2*c3*(m2.e3)*(m2.p*) - c3^2 >= 0
   *
   * @param m2 The rotation axis
   * @param s0 The incident beam vector
   * @param r The reflection
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the small angle approximation is valid
   */
  inline
  bool is_xds_small_angle_valid(vec3<double> m2, vec3<double> s0,
      const Reflection &r, double delta_m) {
    return is_xds_small_angle_valid(m2, s0, r.get_beam_vector(), delta_m);
  }

  /**
   * Check if the XDS small angle approximation holds for the local
   * reflection transform.
   *
   * Check that the following condition holds:
   *  (m2.e1)^2 + 2*c3*(m2.e3)*(m2.p*) - c3^2 >= 0
   *
   * @param g The goniometer
   * @param b The beam
   * @param r The reflection
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the small angle approximation is valid
   */
  inline
  bool is_xds_small_angle_valid(const Goniometer &g, const Beam &b,
      const Reflection &r, double delta_m) {
    return is_xds_small_angle_valid(g.get_rotation_axis(), b.get_s0(),
      r, delta_m);
  }

  /**
   * Check that the angle can be mapped to the local reflection coordinate
   * system
   * @param m2 The rotation axis
   * @param s0 The incident beam vector
   * @param s1 The diffracted beam vector
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the angle is valid
   */
  inline
  bool is_xds_angle_valid(vec3<double> m2, vec3<double> s0, vec3<double> s1,
      double delta_m) {
    vec3<double> ps = (s1 - s0).normalize();
    vec3<double> e1 = s1.cross(s0).normalize();
    vec3<double> e3 = (s1 + s0).normalize();
    double m2e1 = m2 * e1;
    double m2e3 = m2 * e3;
    double m2ps = m2 * ps;
    double m2e3_m2ps = m2e3 * m2ps;
    if (m2e1 == 0) {
      return false;
    }
    double rt = std::sqrt(m2e1*m2e1 + m2e3_m2ps*m2e3_m2ps);
    double tandphi0 = (m2e3_m2ps + rt) / m2e1;
    double tandphi1 = (m2e3_m2ps - rt) / m2e1;
    double dphi0 = 2.0 * std::atan(tandphi0);
    double dphi1 = 2.0 * std::atan(tandphi1);
    if (dphi0 > dphi1) {
      std::swap(dphi0, dphi1);
    }
    delta_m = std::abs(delta_m);
    return dphi0 <= -delta_m && dphi1 >= delta_m;
  }

  /**
   * Check that the angle can be mapped to the local reflection coordinate
   * system
   * @param cs The local reflection coordinate system
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the angle is valid
   */
  inline
  bool is_xds_angle_valid(const CoordinateSystem &cs, double delta_m) {
    vec3<double> m2 = cs.m2();
    vec3<double> ps = cs.p_star().normalize();
    vec3<double> e1 = cs.e1_axis();
    vec3<double> e3 = cs.e3_axis();
    double m2e1 = m2 * e1;
    double m2e3 = m2 * e3;
    double m2ps = m2 * ps;
    double m2e3_m2ps = m2e3 * m2ps;
    if (m2e1 == 0) {
      return false;
    }
    double rt = std::sqrt(m2e1*m2e1 + m2e3_m2ps*m2e3_m2ps);
    double tandphi0 = (m2e3_m2ps + rt) / m2e1;
    double tandphi1 = (m2e3_m2ps - rt) / m2e1;
    double dphi0 = 2.0 * std::atan(tandphi0);
    double dphi1 = 2.0 * std::atan(tandphi1);
    if (dphi0 > dphi1) {
      std::swap(dphi0, dphi1);
    }
    delta_m = std::abs(delta_m);
    return dphi0 <= -delta_m && dphi1 >= delta_m;
  }

  /**
   * Check that the angle can be mapped to the local reflection coordinate
   * system
   * @param m2 The rotation axis
   * @param s0 The incident beam vector
   * @param r The reflection
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the angle is valid
   */
  inline
  bool is_xds_angle_valid(vec3<double> m2, vec3<double> s0,
      const Reflection &r, double delta_m) {
    return is_xds_angle_valid(m2, s0, r.get_beam_vector(), delta_m);
  }

  /**
   * Check that the angle can be mapped to the local reflection coordinate
   * system
   * @param g The goniometer
   * @param s0 The beam
   * @param r The reflection
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the angle is valid
   */
  inline
  bool is_xds_angle_valid(const Goniometer &g, const Beam &b,
      const Reflection &r, double delta_m) {
    return is_xds_angle_valid(g.get_rotation_axis(), b.get_s0(), r, delta_m);
  }

  /**
   * Filter the reflection list by the value of zeta. Set any reflections
   * below the value to invalid.
   * @param g The goniometer
   * @param b The beam
   * @param r The reflection list
   * @param min_zeta The minimum zeta value
   */
  void by_zeta(const Goniometer &g, const Beam &b, ReflectionList &r,
      double min_zeta) {
    for (std::size_t i = 0; i < r.size(); ++i) {
      if (!is_zeta_valid(g, b, r[i], min_zeta)) {
        r[i].set_valid(false);
      }
    }
  }

  /**
   * Filter the reflection list by the validity of the xds small angle approx.
   * Set any reflections for which its invalid to invalid.
   * @param g The goniometer
   * @param b The beam
   * @param r The reflection list
   * @param delta_m The mosaicity * n_sigma
   */
  void by_xds_small_angle(const Goniometer &g, const Beam &b,
      ReflectionList &r, double delta_m) {
    for (std::size_t i = 0; i < r.size(); ++i) {
      if (!is_xds_small_angle_valid(g, b, r[i], delta_m)) {
        r[i].set_valid(false);
      }
    }
  }

  /**
   * Filter the reflection list by the validity of the xds angle.
   * Set any reflections for which its invalid to invalid.
   * @param g The goniometer
   * @param b The beam
   * @param r The reflection list
   * @param delta_m The mosaicity * n_sigma
   */
  void by_xds_angle(const Goniometer &g, const Beam &b,
      ReflectionList &r, double delta_m) {
    for (std::size_t i = 0; i < r.size(); ++i) {
      if (!is_xds_angle_valid(g, b, r[i], delta_m)) {
        r[i].set_valid(false);
      }
    }
  }

  /**
   * Filter the reflection list based on the bounding box volume
   * @param reflections The list of reflections
   */
  void by_bbox_volume(ReflectionList &reflections, std::size_t num_bins) {

    // Check the bins are correct
    DIALS_ASSERT(num_bins > 0);

    // Calculate the bounding box volume for all reflections and then
    // find the minimum and maximum volumes
    flex_int volume(reflections.size());
    int min_volume = std::numeric_limits<int>::max(), max_volume = 0;
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      int6 bbox = reflections[i].get_bounding_box();
      volume[i] = (bbox[1]-bbox[0]) * (bbox[3]-bbox[2]) * (bbox[5]-bbox[4]);
      if (volume[i] < min_volume) min_volume = volume[i];
      if (volume[i] > max_volume) max_volume = volume[i];
    }

    // Check that the volumes are valid
    DIALS_ASSERT(max_volume > min_volume && min_volume > 0 && max_volume > 0);

    // Create the volume histogram
    flex_double histo(num_bins);
    double bin_size = (float)(max_volume - min_volume) / (float)(num_bins - 1);
    for (std::size_t i = 0; i < volume.size(); ++i) {
      int index = (int)((volume[i] - min_volume) / bin_size);
      if (index < 0) index = 0;
      if (index >= num_bins) index = num_bins - 1;
      histo[(int)((volume[i] - min_volume) / bin_size)]++;
    }

    // Calculate the threshold and set any reflections with bounding
    // box size greater than the threshold to be invalid.
    double threshold = maximum_deviation(histo) * bin_size;
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      if (volume[i] > threshold) {
        reflections[i].set_valid(false);
      }
    }
  }

  /**
   * Filter the reflections by the bounding box volume. Use a histogram with
   * nbins = cube_root(nref)
   * @param reflections The reflections list
   */
  void by_bbox_volume(ReflectionList &reflections) {
    std::size_t num = std::exp((1.0/3.0) * std::log(reflections.size()));
    return by_bbox_volume(reflections, num);
  }

  /**
   * Check if the bounding box has points outside the image range.
   * @param bbox The bounding box
   * @param image_size The image size
   * @param scan_range The scan range
   * @returns True/False
   */
  bool is_bbox_outside_image_range(int6 bbox, small<long,10> image_size,
      int2 scan_range) {
    DIALS_ASSERT(image_size.size() == 2);
    return bbox[0] < 0 || bbox[1] >= image_size[1] ||
           bbox[2] < 0 || bbox[3] >= image_size[0] ||
           bbox[4] < scan_range[0] || bbox[5] >= scan_range[1];
  }

  /**
   * Check if the bounding box has points that cover bad pixels
   * @param bbox The bounding box
   * @param mask The mask array
   * @returns True/False
   */
  bool does_bbox_contain_bad_pixels(int6 bbox, const flex_bool &mask) {
    DIALS_ASSERT(mask.accessor().all().size() == 2);
    for (int j = bbox[2]; j < bbox[3]; ++j) {
      for (int i = bbox[0]; i < bbox[1]; ++i) {
        if (mask(j, i) == false) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Check if the bounding box is valid in terms of the detector mask
   * @param bbox The bounding box
   * @param mask The mask array
   * @param scan_range The scan range
   * @returns True/False
   */
  bool is_bbox_valid(int6 bbox, const flex_bool &mask, int2 scan_range) {
    return !(is_bbox_outside_image_range(
        bbox, mask.accessor().all(), scan_range) ||
        does_bbox_contain_bad_pixels(bbox, mask));
  }

  /**
   * Filter the reflection based on the detector mask
   * @param reflection The reflection
   * @param mask The detector mask
   * @param scan_range The scan range
   */
  void by_detector_mask(Reflection &reflection, const flex_bool &mask,
      int2 scan_range) {
    // get the bounding box
    int6 bbox = reflection.get_bounding_box();

    // Set whether the reflection is valid or not
    reflection.set_valid(is_bbox_valid(bbox, mask, scan_range));
  }

  /**
   * Filter the reflection list based on the detector mask
   * @param reflection The reflection
   * @param mask The detector mask
   * @param scan_range The scan range
   */
  void by_detector_mask(ReflectionList &reflections,
      const flex_bool &mask, int2 scan_range) {
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      by_detector_mask(reflections[i], mask, scan_range);
    }
  }

}}} // namespace dials::algorithms::filter

#endif /* DIALS_ALGORITHMS_FILTER_H */
