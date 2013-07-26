/*
 * forward.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_REFLEXION_BASIS_FORWARD_H
#define DIALS_ALGORITHMS_REFLEXION_BASIS_FORWARD_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>
#include <dials/algorithms/reflection_basis/beam_vector_map.h>
#include <dials/algorithms/reflection_basis/map_frames.h>
#include <dials/algorithms/reflection_basis/map_pixels.h>

namespace dials { namespace algorithms { namespace reflection_basis {
  namespace transform {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using scitbx::af::flex_double;
  using scitbx::af::flex_bool;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;

  /**
   * Class to perform the forward reflection basis transform
   */
  class Forward {
  public:

    /**
     * Initialise the class
     * @param beam The beam model
     * @param detector The detector model
     * @param gonio The goniometer model
     * @param scan The scan model
     * @param mosaicity The crystal mosaicity
     * @param n_sigma The number of standard deviations
     * @param grid_size The size of the reflection basis grid
     */
    Forward(const Beam &beam, const Detector &detector, const Goniometer &gonio,
            const Scan &scan, double mosaicity, double n_sigma,
            std::size_t grid_size)
      : s0_(beam.get_s0()),
        m2_(gonio.get_rotation_axis()),
        map_frames_(
          scan.get_oscillation()[0],
          scan.get_oscillation()[1],
          mosaicity,
          n_sigma,
          grid_size),
        map_pixels_(
          beam_vector_map(detector, beam, true),
          grid_size,
          vec2<double>(
            beam.get_sigma_divergence() * n_sigma / grid_size,
            beam.get_sigma_divergence() * n_sigma / grid_size)) {
      DIALS_ASSERT(n_sigma > 0);
    }

    /**
     * Transform the reflection into its own recirprocal space basis grid
     * @param cs The coordinate system
     * @param bbox The bounding box
     * @param image The image
     * @param mask The mask
     * @returns The grid
     */
    flex_double operator()(const CoordinateSystem &cs, int6 bbox,
        const flex_double &image, const flex_bool &mask) const {

      // Calculate the fraction of intensity contributed from each data
      // frame to each grid coordinate
      flex_double zfraction = map_frames_(vec2<int>(bbox[4], bbox[5]),
        cs.phi(), cs.zeta());

      // Map the image pixels to the grid points
      return map_pixels_(cs, bbox, image, mask, zfraction);
    }

    /**
     * Transform the reflection into its own recirprocal space basis grid
     * @param s1 The beam vector
     * @param phi The rotation angle
     * @param bbox The bounding box
     * @param image The image
     * @param mask The mask
     * @returns The grid
     */
    flex_double operator()(vec3<double> s1, double phi, int6 bbox,
        const flex_double &image, const flex_bool &mask) const {
      CoordinateSystem cs(m2_, s0_, s1, phi);
      return this->operator()(cs, bbox, image, mask);
    }

  private:
    vec3<double> s0_;
    vec3<double> m2_;
    MapFramesForward map_frames_;
    MapPixelsForward map_pixels_;
  };

  /**
   * Class to perform the reverse reflection basis transform
   */
  class Reverse {
  public:

    /**
     * Initialise the class
     * @param beam The beam model
     * @param detector The detector model
     * @param gonio The goniometer model
     * @param scan The scan model
     * @param mosaicity The crystal mosaicity
     * @param n_sigma The number of standard deviations
     * @param grid_size The size of the reflection basis grid
     */
    Reverse(const Beam &beam, const Detector &detector, const Goniometer &gonio,
            const Scan &scan, double mosaicity, double n_sigma,
            std::size_t grid_size)
      : s0_(beam.get_s0()),
        m2_(gonio.get_rotation_axis()),
        map_frames_(
          scan.get_oscillation()[0],
          scan.get_oscillation()[1],
          mosaicity,
          n_sigma,
          grid_size),
        map_pixels_(
          beam_vector_map(detector, beam, true),
          grid_size,
          vec2<double>(
            beam.get_sigma_divergence() * n_sigma / grid_size,
            beam.get_sigma_divergence() * n_sigma / grid_size)) {
      DIALS_ASSERT(n_sigma > 0);
    }

    /**
     * Transform the reflection into its own recirprocal space basis grid
     * @param cs The coordinate system
     * @param bbox The bounding box
     * @param grid The grid
     * @returns The image
     */
    flex_double operator()(const CoordinateSystem &cs, int6 bbox,
        const flex_double &grid) const {

      // Calculate the fraction of intensity contributed from each data
      // frame to each grid coordinate
      flex_double zfraction = map_frames_(vec2<int>(bbox[4], bbox[5]),
        cs.phi(), cs.zeta());

      // Map the image pixels from the grid points
      return map_pixels_(cs, bbox, grid, zfraction);
    }

    /**
     * Transform the reflection into its own recirprocal space basis grid
     * @param s1 The beam vector
     * @param phi The rotation angle
     * @param bbox The bounding box
     * @param grid The grid
     * @returns The image
     */
    flex_double operator()(vec3<double> s1, double phi, int6 bbox,
        const flex_double &grid) const {
      CoordinateSystem cs(m2_, s0_, s1, phi);
      return this->operator()(cs, bbox, grid);
    }

  private:
    vec3<double> s0_;
    vec3<double> m2_;
    MapFramesReverse map_frames_;
    MapPixelsReverse map_pixels_;
  };

}}}} // namespace dials::algorithms::reflection_basis::transform

#endif /* DIALS_ALGORITHMS_REFLEXION_BASIS_FORWARD_H */
