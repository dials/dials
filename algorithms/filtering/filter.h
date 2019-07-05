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
#include <limits>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/image/threshold/unimodal.h>
#include <dials/algorithms/profile_model/gaussian_rs/coordinate_system.h>

namespace dials { namespace algorithms { namespace filter {

  using dials::algorithms::profile_model::gaussian_rs::CoordinateSystem;
  using dials::algorithms::profile_model::gaussian_rs::zeta_factor;
  using dials::model::Foreground;
  using dials::model::Shoebox;
  using dials::model::Valid;
  using dxtbx::model::BeamBase;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using scitbx::af::small;
  using scitbx::af::tiny;

  /**
   * Calculate the zeta factor and check its absolute value is above the
   * minimum specified value.
   * @param m2 The rotation axis (normalized)
   * @param s0 The incident beam vector
   * @param s1 The diffracted beam vector
   * @param zeta_min The minimum allowed zeta value
   * @returns True/False, zeta is valid
   */
  inline bool is_zeta_valid(vec3<double> m2,
                            vec3<double> s0,
                            vec3<double> s1,
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
  inline bool is_zeta_valid(const CoordinateSystem &cs, double zeta_min) {
    return std::abs(cs.zeta()) >= zeta_min;
  }

  /**
   * Calculate the zeta factor and check its absolute value is above the
   * minimum specified value.
   * @param g The goniometer
   * @param b The beam
   * @param s1 The beam vector
   * @param zeta_min The minimum allowed zeta value
   * @returns True/False, zeta is valid
   */
  inline bool is_zeta_valid(const Goniometer &g,
                            const BeamBase &b,
                            vec3<double> s1,
                            double zeta_min) {
    return is_zeta_valid(g.get_rotation_axis(), b.get_s0(), s1, zeta_min);
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
  inline bool is_xds_small_angle_valid(vec3<double> m2,
                                       vec3<double> s0,
                                       vec3<double> s1,
                                       double delta_m) {
    vec3<double> ps = (s1 - s0).normalize();
    vec3<double> e1 = s1.cross(s0).normalize();
    vec3<double> e3 = (s1 + s0).normalize();
    double m2e1 = m2 * e1;
    double m2e3 = m2 * e3;
    double m2ps = m2 * ps;
    double c3 = -std::abs(delta_m);
    return (m2e1 * m2e1 + 2.0 * c3 * m2e3 * m2ps - c3 * c3) >= 0.0;
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
  inline bool is_xds_small_angle_valid(const CoordinateSystem &cs, double delta_m) {
    vec3<double> m2 = cs.m2();
    vec3<double> ps = cs.p_star().normalize();
    vec3<double> e1 = cs.e1_axis();
    vec3<double> e3 = cs.e3_axis();
    double m2e1 = m2 * e1;
    double m2e3 = m2 * e3;
    double m2ps = m2 * ps;
    double c3 = -std::abs(delta_m);
    return (m2e1 * m2e1 + 2.0 * c3 * m2e3 * m2ps - c3 * c3) >= 0.0;
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
   * @param s1 The beam vector
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the small angle approximation is valid
   */
  inline bool is_xds_small_angle_valid(const Goniometer &g,
                                       const BeamBase &b,
                                       vec3<double> s1,
                                       double delta_m) {
    return is_xds_small_angle_valid(g.get_rotation_axis(), b.get_s0(), s1, delta_m);
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
  inline bool is_xds_angle_valid(vec3<double> m2,
                                 vec3<double> s0,
                                 vec3<double> s1,
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
    double rt = std::sqrt(m2e1 * m2e1 + m2e3_m2ps * m2e3_m2ps);
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
  inline bool is_xds_angle_valid(const CoordinateSystem &cs, double delta_m) {
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
    double rt = std::sqrt(m2e1 * m2e1 + m2e3_m2ps * m2e3_m2ps);
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
   * @param g The goniometer
   * @param s0 The beam
   * @param s1 The diffracted beam vector
   * @param delta_m The mosaicity * n_sigma
   * @returns True/False, the angle is valid
   */
  inline bool is_xds_angle_valid(const Goniometer &g,
                                 const BeamBase &b,
                                 vec3<double> s1,
                                 double delta_m) {
    return is_xds_angle_valid(g.get_rotation_axis(), b.get_s0(), s1, delta_m);
  }

  /**
   * Filter the reflection list by the value of zeta. Set any reflections
   * below the value to invalid.
   * @param g The goniometer
   * @param b The beam
   * @param s1 The list of beam vectors
   * @param min_zeta The minimum zeta value
   */
  inline af::shared<bool> by_zeta(const Goniometer &g,
                                  const BeamBase &b,
                                  const af::const_ref<vec3<double> > &s1,
                                  double min_zeta) {
    af::shared<bool> result(s1.size(), true);
    for (std::size_t i = 0; i < result.size(); ++i) {
      if (!is_zeta_valid(g, b, s1[i], min_zeta)) {
        result[i] = false;
      }
    }
    return result;
  }

  /**
   * Filter the reflection list by the validity of the xds small angle approx.
   * Set any reflections for which its invalid to invalid.
   * @param g The goniometer
   * @param b The beam
   * @param s1 The list of beam vector
   * @param delta_m The mosaicity * n_sigma
   */
  inline af::shared<bool> by_xds_small_angle(const Goniometer &g,
                                             const BeamBase &b,
                                             const af::const_ref<vec3<double> > s1,
                                             double delta_m) {
    af::shared<bool> result(s1.size(), true);
    for (std::size_t i = 0; i < result.size(); ++i) {
      if (!is_xds_small_angle_valid(g, b, s1[i], delta_m)) {
        result[i] = false;
      }
    }
    return result;
  }

  /**
   * Filter the reflection list by the validity of the xds angle.
   * Set any reflections for which its invalid to invalid.
   * @param g The goniometer
   * @param b The beam
   * @param s1 The list of beam vectors
   * @param delta_m The mosaicity * n_sigma
   */
  inline af::shared<bool> by_xds_angle(const Goniometer &g,
                                       const BeamBase &b,
                                       const af::const_ref<vec3<double> > s1,
                                       double delta_m) {
    af::shared<bool> result(s1.size(), true);
    for (std::size_t i = 0; i < result.size(); ++i) {
      if (!is_xds_angle_valid(g, b, s1[i], delta_m)) {
        result[i] = false;
      }
    }
    return result;
  }

  /**
   * Filter the reflection list based on the bounding box volume
   * @param bboxes The list of bounding boxes
   */
  inline af::shared<bool> by_bbox_volume(const af::const_ref<int6> &bboxes,
                                         std::size_t num_bins) {
    // Check the bins are correct
    DIALS_ASSERT(num_bins > 0);

    // Calculate the bounding box volume for all reflections and then
    // find the minimum and maximum volumes
    af::shared<int> volume(bboxes.size(), af::init_functor_null<int>());
    int min_volume = std::numeric_limits<int>::max(), max_volume = 0;
    for (std::size_t i = 0; i < bboxes.size(); ++i) {
      int6 bbox = bboxes[i];
      volume[i] = (bbox[1] - bbox[0]) * (bbox[3] - bbox[2]) * (bbox[5] - bbox[4]);
      if (volume[i] < min_volume) min_volume = volume[i];
      if (volume[i] > max_volume) max_volume = volume[i];
    }

    // Check that the volumes are valid
    DIALS_ASSERT(max_volume > min_volume && min_volume > 0 && max_volume > 0);

    // Create the volume histogram
    af::shared<double> histo(num_bins, af::init_functor_null<double>());
    double bin_size = (float)(max_volume - min_volume) / (float)(num_bins - 1);
    for (std::size_t i = 0; i < volume.size(); ++i) {
      int index = (int)((volume[i] - min_volume) / bin_size);
      if (index < 0) index = 0;
      if (index >= num_bins) index = num_bins - 1;
      histo[(int)((volume[i] - min_volume) / bin_size)]++;
    }

    // Calculate the threshold and set any reflections with bounding
    // box size greater than the threshold to be invalid.
    double threshold = maximum_deviation(histo.const_ref()) * bin_size;
    af::shared<bool> result(bboxes.size(), true);
    for (std::size_t i = 0; i < bboxes.size(); ++i) {
      if (volume[i] > threshold) {
        result[i] = false;
      }
    }
    return result;
  }

  /**
   * Filter the reflections by the bounding box volume. Use a histogram with
   * nbins = cube_root(nref)
   * @param bboxes The list of bounding boxes
   */
  inline af::shared<bool> by_bbox_volume(const af::const_ref<int6> &bboxes) {
    std::size_t num =
      (std::size_t)(std::exp((1.0 / 3.0) * std::log((double)bboxes.size())));
    return by_bbox_volume(bboxes, num);
  }

  /**
   * Check if the bounding box has points outside the image range.
   * @param bbox The bounding box
   * @param image_size The image size
   * @param scan_range The scan range
   * @returns True/False
   */
  inline bool is_bbox_outside_image_range(int6 bbox,
                                          tiny<std::size_t, 2> image_size,
                                          int2 scan_range) {
    DIALS_ASSERT(image_size.size() == 2);
    return bbox[0] < 0 || bbox[1] > (int)image_size[1] || bbox[2] < 0
           || bbox[3] > (int)image_size[0] || bbox[4] < scan_range[0]
           || bbox[5] > scan_range[1];
  }

  /**
   * Check if the bounding box has points that cover bad pixels
   * @param bbox The bounding box
   * @param mask The mask array
   * @returns True/False
   */
  inline bool does_bbox_contain_bad_pixels(
    int6 bbox,
    const af::const_ref<bool, af::c_grid<2> > &mask) {
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
  inline bool is_bbox_valid(int6 bbox,
                            const af::const_ref<bool, af::c_grid<2> > &mask,
                            int2 scan_range) {
    return !(is_bbox_outside_image_range(bbox, mask.accessor(), scan_range)
             || does_bbox_contain_bad_pixels(bbox, mask));
  }

  /**
   * Filter the reflection list based on the detector mask
   * @param bboxes The list of bounding boxes
   * @param mask The detector mask
   * @param scan_range The scan range
   */
  inline af::shared<bool> by_detector_mask(
    const af::const_ref<int6> bboxes,
    const af::const_ref<bool, af::c_grid<2> > &mask,
    int2 scan_range) {
    af::shared<bool> result(bboxes.size());
    for (std::size_t i = 0; i < bboxes.size(); ++i) {
      result[i] = is_bbox_valid(bboxes[i], mask, scan_range);
    }
    return result;
  }

  inline af::shared<bool> by_detector_mask_multipanel(
    const af::const_ref<std::size_t> &panel,
    const af::const_ref<int6> bboxes,
    const af::const_ref<af::const_ref<bool, af::c_grid<2> > > &mask,
    int2 scan_range) {
    DIALS_ASSERT(panel.size() == bboxes.size());
    af::shared<bool> result(bboxes.size());
    for (std::size_t i = 0; i < bboxes.size(); ++i) {
      DIALS_ASSERT(panel[i] < mask.size());
      result[i] = is_bbox_valid(bboxes[i], mask[panel[i]], scan_range);
    }
    return result;
  }

  /**
   * Filter the reflection list by the distance between the centroid
   * position and predicted position.
   * @param reflections The reflection list
   * @param max_separation The maximum allowed separation
   */
  inline af::shared<bool> by_centroid_prediction_separation(
    const af::const_ref<vec3<double> > &xyzobs,
    const af::const_ref<vec3<double> > &xyzcal,
    double max_separation) {
    DIALS_ASSERT(xyzobs.size() == xyzcal.size());
    af::shared<bool> result(xyzobs.size(), true);
    for (std::size_t i = 0; i < result.size(); ++i) {
      vec3<double> obs = xyzobs[i];
      vec3<double> cal = xyzcal[i];
      double sep = std::sqrt((obs[0] - cal[0]) * (obs[0] - cal[0])
                             + (obs[1] - cal[1]) * (obs[1] - cal[1])
                             + (obs[2] - cal[2]) * (obs[2] - cal[2]));
      if (sep > max_separation) {
        result[i] = false;
      }
    }
    return result;
  }

  /**
   * Filter the reflection list by resolution
   * position and predicted position.
   * @param reflections The reflection list
   * @param beam The beam model
   * @param detector The detector model
   * @param d_min The maximum resolution
   * @param d_max The minimum resolution
   */
  inline af::shared<bool> by_resolution_at_centroid(
    const af::const_ref<std::size_t> &panel,
    const af::const_ref<vec3<double> > &xyz,
    const BeamBase &beam,
    const Detector &detector,
    double d_min,
    double d_max) {
    DIALS_ASSERT(panel.size() == xyz.size());
    vec3<double> s0 = beam.get_s0();
    if (d_max < 0) {
      d_max = std::numeric_limits<double>::max();
    }
    af::shared<bool> result(xyz.size(), true);
    for (std::size_t i = 0; i < result.size(); ++i) {
      DIALS_ASSERT(panel[i] < detector.size());
      double resolution = detector[panel[i]].get_resolution_at_pixel(
        s0, vec2<double>(xyz[i][0], xyz[i][1]));
      if (resolution < d_min || resolution > d_max) {
        result[i] = false;
      }
    }
    return result;
  }

  inline af::shared<bool> by_shoebox_mask(const af::const_ref<Shoebox<> > shoeboxes,
                                          int2 image_size,
                                          int2 scan_range) {
    af::shared<bool> result(shoeboxes.size(), true);
    for (std::size_t i = 0; i < shoeboxes.size(); ++i) {
      Shoebox<> sbox = shoeboxes[i];
      DIALS_ASSERT(sbox.is_consistent());
      if (is_bbox_outside_image_range(sbox.bbox, image_size, scan_range)) {
        result[i] = false;
        continue;
      }
      for (std::size_t j = 0; j < sbox.mask.size(); ++j) {
        if (sbox.mask[j] & Foreground) {
          if (!(sbox.mask[j] & Valid)) {
            result[i] = false;
            break;
          }
        }
      }
    }
    return result;
  }

}}}  // namespace dials::algorithms::filter

#endif /* DIALS_ALGORITHMS_FILTER_H */
