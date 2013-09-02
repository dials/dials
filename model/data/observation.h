/*
 * observation.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_DATA_OBSERVATION_H
#define DIALS_MODEL_DATA_OBSERVATION_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <dials/error.h>

namespace dials { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using dxtbx::model::Detector;
  using dxtbx::model::Scan;

  /**
   * A structure to hold the intensity data we want for both raw and
   * corrected intensities.
   */
  struct Intensity {

    /**
     * A struct with the intensity value and variance
     */
    struct IntensityData {
      double value;
      double variance;

      /** Default construct */
      IntensityData()
        : value(0.0),
          variance(0.0) {}

      /** Construct with values */
      IntensityData(double value_, double variance_)
        : value(value_),
          variance(variance_) {}
    };

    IntensityData observed;
    IntensityData corrected;

    /** Default construct */
    Intensity() {}

    /** Construct with observed */
    Intensity(double observed_value, double observed_variance)
      : observed(observed_value, observed_variance) {}

    /** Construct with observed and corrected */
    Intensity(double observed_value, double observed_variance,
              double corrected_value, double corrected_variance)
      : observed(observed_value, observed_variance),
        corrected(corrected_value, corrected_variance) {}

    /** Construct with observed */
    Intensity(const IntensityData &observed_)
      : observed(observed_) {}

    /** Construct with observed and corrected */
    Intensity(const IntensityData &observed_, const IntensityData &corrected_)
      : observed(observed_),
        corrected(corrected_) {}
  };

  /**
   * A struct containing the centroid position and variance in both pixel
   * and millimeter coordinates.
   */
  struct Position {

    /**
     * The centroid data
     */
    struct PositionData {
      vec3<double> position;
      vec3<double> variance;
      vec3<double> std_err_sq;

      /** Default construct */
      PositionData()
        : position(0, 0, 0),
          variance(0, 0, 0),
          std_err_sq(0, 0, 0) {}

      /** Construct from values */
      PositionData(vec3<double> position_,
                   vec3<double> variance_,
                   vec3<double> std_err_sq_)
        : position(position_),
          variance(variance_),
          std_err_sq(std_err_sq_) {}
    };

    PositionData px;
    PositionData mm;

    /** Default construct */
    Position() {}

    /** Construct with pixel coordinates */
    Position(vec3<double> px_position,
             vec3<double> px_variance,
             vec3<double> px_std_err_sq)
      : px(px_position, px_variance, px_std_err_sq) {}

    /** Construct with pixel and millimeter coordinates */
    Position(vec3<double> px_position,
             vec3<double> px_variance,
             vec3<double> px_std_err_sq,
             vec3<double> mm_position,
             vec3<double> mm_variance,
             vec3<double> mm_std_err_sq)
      : px(px_position, px_variance, px_std_err_sq),
        mm(mm_position, mm_variance, mm_std_err_sq) {}

    /** Construct with the pixel coordinate */
    Position(const PositionData &px_)
      : px(px_) {}

    /** Construct the with pixel and millimetre position */
    Position(const PositionData &px_, const PositionData &mm_)
      : px(px_),
        mm(mm_) {}

    /**
     * Update the millimeter centroid position from the pixel centroid
     * position using the given geometry.
     * @param panel The panel
     * @param d The detector model
     * @param s The scan model
     */
    void update_mm(std::size_t panel, const Detector &d, const Scan &s) {

      // Check the panel number
      DIALS_ASSERT(panel < d.num_panels());

      // Get the milliemeter x, y, z coordinates
      vec2<double> px_xy = vec2<double>(px.position[0], px.position[1]);
      vec2<double> mm_xy = d[panel].pixel_to_millimeter(px_xy);
      double px_z = px.position[2];
      double mm_z = s.get_angle_from_array_index(px_z);

      // Set the millimeter position
      mm.position[0] = mm_xy[0];
      mm.position[1] = mm_xy[1];
      mm.position[2] = mm_z;

      // Scale the variance and standard error squared
      vec2<double> pixel_size = d[panel].get_pixel_size();
      vec2<double> oscillation = s.get_oscillation();
      vec3<double> scale(pixel_size[0], pixel_size[1], oscillation[1]);
      for (std::size_t i = 0; i < 3; ++i) {
        mm.variance[i] = px.variance[i] * scale[i];
        mm.std_err_sq[i] = px.std_err_sq[i] * scale[i];
      }
    }

    /**
     * Update the millimeter centroid position from the pixel centroid
     * position using the given geometry.
     * @param d The detector model
     * @param s The scan model
     */
    void update_mm(const Detector &d, const Scan &s) {
      update_mm(0, d, s);
    }
  };

  /**
   * A struct holding details about an observation
   */
  struct Observation {

    std::size_t panel;
    Position centroid;
    Intensity intensity;

    /** Default construct */
    Observation()
      : panel(0) {}

    /** Construct with position */
    Observation(const Position &centroid_)
      : panel(0),
        centroid(centroid_) {}

    /** Construct with intensity */
    Observation(const Intensity &intensity_)
      : panel(0),
        intensity(intensity_) {}

    /** Construct with position and intensity */
    Observation(const Position &centroid_, const Intensity &intensity_)
      : panel(0),
        centroid(centroid_),
        intensity(intensity_) {}

    /** Construct with panel */
    Observation(std::size_t panel_)
      : panel(panel_) {}

    /** Construct with position */
    Observation(std::size_t panel_, const Position &centroid_)
      : panel(panel_),
        centroid(centroid_) {}

    /** Construct with intensity */
    Observation(std::size_t panel_, const Intensity &intensity_)
      : panel(panel_),
        intensity(intensity_) {}

    /** Construct with position and intensity */
    Observation(std::size_t panel_, const Position &centroid_,
      const Intensity &intensity_)
      : panel(panel_),
        centroid(centroid_),
        intensity(intensity_) {}

    /**
     * Update the millimeter centroid position from the pixel centroid
     * position using the given geometry.
     * @param d The detector model
     * @param s The scan model
     */
    void update_centroid_mm(const Detector &d, const Scan &s) {
      centroid.update_mm(d, s);
    }
  };

}}; // namespace dials::model

#endif /* DIALS_MODEL_DATA_OBSERVATION_H */
