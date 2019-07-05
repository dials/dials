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
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/scan.h>
#include <dials/error.h>

namespace dials { namespace model {

  using dxtbx::model::BeamBase;
  using dxtbx::model::Detector;
  using dxtbx::model::Scan;
  using scitbx::vec2;
  using scitbx::vec3;

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
      bool success;

      /** Default construct */
      IntensityData() : value(0.0), variance(0.0), success(false) {}

      /** Construct with values */
      IntensityData(double value_, double variance_, bool success_)
          : value(value_), variance(variance_), success(success_) {}

      /**
       * Test to see if intensity contain the same data
       * @param rhs The other intensity
       * @returns True/False. They are the same
       */
      bool operator==(const IntensityData &rhs) const {
        const double eps = 1e-7;
        return (std::abs(value - rhs.value) < eps
                && std::abs(variance - rhs.variance) < eps && success == rhs.success);
      }

      /**
       * Test to see if intensity contain the same data
       * @param rhs The other intensity
       * @returns True/False. They are the same
       */
      bool operator!=(const IntensityData &rhs) const {
        return !(*this == rhs);
      }
    };

    IntensityData observed;
    IntensityData corrected;
    IntensityData background;

    /** Default construct */
    Intensity() {}

    /** Construct with observed */
    Intensity(double observed_value, double observed_variance, bool observed_success)
        : observed(observed_value, observed_variance, observed_success) {}

    /** Construct with observed and corrected */
    Intensity(double observed_value,
              double observed_variance,
              bool observed_success,
              double corrected_value,
              double corrected_variance,
              bool corrected_success)
        : observed(observed_value, observed_variance, observed_success),
          corrected(corrected_value, corrected_variance, corrected_success) {}

    /** Construct with observed */
    Intensity(const IntensityData &observed_) : observed(observed_) {}

    /** Construct with observed and corrected */
    Intensity(const IntensityData &observed_, const IntensityData &corrected_)
        : observed(observed_), corrected(corrected_) {}

    /**
     * Test to see if intensity contain the same data
     * @param rhs The other intensity
     * @returns True/False. They are the same
     */
    bool operator==(const Intensity &rhs) const {
      return observed == rhs.observed && corrected == rhs.corrected;
    }

    /**
     * Test to see if intensity contain the same data
     * @param rhs The other intensity
     * @returns True/False. They are the same
     */
    bool operator!=(const Intensity &rhs) const {
      return !(*this == rhs);
    }
  };

  /**
   * A struct containing the centroid position and variance in both pixel
   * and millimeter coordinates.
   */
  struct Centroid {
    /**
     * The centroid data
     */
    struct CentroidData {
      vec3<double> position;
      vec3<double> variance;
      vec3<double> std_err_sq;

      /** Default construct */
      CentroidData() : position(0, 0, 0), variance(0, 0, 0), std_err_sq(0, 0, 0) {}

      /** Construct from values */
      CentroidData(vec3<double> position_,
                   vec3<double> variance_,
                   vec3<double> std_err_sq_)
          : position(position_), variance(variance_), std_err_sq(std_err_sq_) {}

      /**
       * Test to see if centroids contain the same data
       * @param rhs The other centroid
       * @returns True/False. They are the same
       */
      bool operator==(const CentroidData &rhs) const {
        const double eps = 1e-7;
        return ((std::abs(position[0] - rhs.position[0]) < eps)
                && (std::abs(position[1] - rhs.position[1]) < eps)
                && (std::abs(position[2] - rhs.position[2]) < eps)
                && (std::abs(variance[0] - rhs.variance[0]) < eps)
                && (std::abs(variance[1] - rhs.variance[1]) < eps)
                && (std::abs(variance[2] - rhs.variance[2]) < eps)
                && (std::abs(std_err_sq[0] - rhs.std_err_sq[0]) < eps)
                && (std::abs(std_err_sq[1] - rhs.std_err_sq[1]) < eps)
                && (std::abs(std_err_sq[2] - rhs.std_err_sq[2]) < eps));
      }

      /**
       * Test to see if centroids contain the same data
       * @param rhs The other centroid
       * @returns True/False. They are the same
       */
      bool operator!=(const CentroidData &rhs) const {
        return !(*this == rhs);
      }
    };

    CentroidData px;
    CentroidData mm;

    /** Default construct */
    Centroid() {}

    /** Construct with pixel coordinates */
    Centroid(vec3<double> px_position,
             vec3<double> px_variance,
             vec3<double> px_std_err_sq)
        : px(px_position, px_variance, px_std_err_sq) {}

    /** Construct with pixel and millimeter coordinates */
    Centroid(vec3<double> px_position,
             vec3<double> px_variance,
             vec3<double> px_std_err_sq,
             vec3<double> mm_position,
             vec3<double> mm_variance,
             vec3<double> mm_std_err_sq)
        : px(px_position, px_variance, px_std_err_sq),
          mm(mm_position, mm_variance, mm_std_err_sq) {}

    /** Construct with the pixel coordinate */
    Centroid(const CentroidData &px_) : px(px_) {}

    /** Construct the with pixel and millimetre position */
    Centroid(const CentroidData &px_, const CentroidData &mm_) : px(px_), mm(mm_) {}

    /**
     * Test to see if centroids contain the same data
     * @param rhs The other centroid
     * @returns True/False. They are the same
     */
    bool operator==(const Centroid &rhs) const {
      return px == rhs.px && mm == rhs.mm;
    }

    /**
     * Test to see if centroids contain the same data
     * @param rhs The other centroid
     * @returns True/False. They are the same
     */
    bool operator!=(const Centroid &rhs) const {
      return !(*this == rhs);
    }

    /**
     * Update the millimeter centroid position from the pixel centroid
     * position using the given geometry.
     * @param panel The panel
     * @param d The detector model
     * @param s The scan model
     */
    void update_mm(std::size_t panel, const Detector &d, const Scan &s) {
      // Check the panel number
      DIALS_ASSERT(panel < d.size());

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

    /**
     * Calculate the resolution of the observed reflection
     * @param panel The panel number
     * @param d The detector model
     * @returns The resolution
     */
    double resolution(std::size_t panel, const BeamBase &b, const Detector &d) const {
      return d[panel].get_resolution_at_pixel(
        b.get_s0(), vec2<double>(px.position[0], px.position[1]));
    }

    /**
     * Calculate the resolution of the observed reflection
     * @param d The detector model
     * @returns The resolution
     */
    double resolution(const BeamBase &b, const Detector &d) const {
      return resolution(0, b, d);
    }
  };

  /**
   * A struct holding details about an observation
   */
  struct Observation {
    std::size_t panel;
    Centroid centroid;
    Intensity intensity;

    /** Default construct */
    Observation() : panel(0) {}

    /** Construct with position */
    Observation(const Centroid &centroid_) : panel(0), centroid(centroid_) {}

    /** Construct with intensity */
    Observation(const Intensity &intensity_) : panel(0), intensity(intensity_) {}

    /** Construct with position and intensity */
    Observation(const Centroid &centroid_, const Intensity &intensity_)
        : panel(0), centroid(centroid_), intensity(intensity_) {}

    /** Construct with panel */
    Observation(std::size_t panel_) : panel(panel_) {}

    /** Construct with position */
    Observation(std::size_t panel_, const Centroid &centroid_)
        : panel(panel_), centroid(centroid_) {}

    /** Construct with intensity */
    Observation(std::size_t panel_, const Intensity &intensity_)
        : panel(panel_), intensity(intensity_) {}

    /** Construct with position and intensity */
    Observation(std::size_t panel_,
                const Centroid &centroid_,
                const Intensity &intensity_)
        : panel(panel_), centroid(centroid_), intensity(intensity_) {}

    /**
     * Update the millimeter centroid position from the pixel centroid
     * position using the given geometry.
     * @param d The detector model
     * @param s The scan model
     */
    void update_centroid_mm(const Detector &d, const Scan &s) {
      centroid.update_mm(panel, d, s);
    }

    /**
     * Calculate the resolution of the observed reflection
     * @param d The detector model
     * @returns The resolution
     */
    double resolution(const BeamBase &b, const Detector &d) const {
      return centroid.resolution(panel, b, d);
    }

    /**
     * Test to see if observations contain the same data
     * @param rhs The other observation
     * @returns True/False. They are the same
     */
    bool operator==(const Observation &rhs) const {
      return (panel == rhs.panel && centroid == rhs.centroid
              && intensity == rhs.intensity);
    }

    /**
     * Test to see if observations contain the same data
     * @param rhs The other observation
     * @returns True/False. They are not the same
     */
    bool operator!=(const Observation &rhs) const {
      return !(*this == rhs);
    }
  };

}};  // namespace dials::model

#endif /* DIALS_MODEL_DATA_OBSERVATION_H */
