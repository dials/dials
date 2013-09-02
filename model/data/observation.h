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

#include <scitbx/vec3.h>

namespace dials { namespace model {

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
  };

  /**
   * A struct holding details about an observation
   */
  struct Observation {
    Position position;
    Intensity intensity;

    /** Default construct */
    Observation() {}

    /** Construct with position */
    Observation(const Position &position_)
      : position(position_) {}

    /** Construct with intensity */
    Observation(const Intensity &intensity_)
      : intensity(intensity_) {}

    /** Construct with position and intensity */
    Observation(const Position &position_, const Intensity &intensity_)
      : position(position_),
        intensity(intensity_) {}
  };

}}; // namespace dials::model

#endif /* DIALS_MODEL_DATA_OBSERVATION_H */
