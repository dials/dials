/*
 * prediction.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_DATA_PREDICTION_H
#define DIALS_MODEL_DATA_PREDICTION_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <cctbx/miller.h>

namespace dials { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;

  typedef cctbx::miller::index<> MillerIndex;

  /**
   * A struct to contain the predictions
   */
  struct Prediction {

    /** Position data */
    struct PositionData {
      vec3<double> px;
      vec3<double> mm;

      /** Default construct */
      PositionData()
        : px(0.0, 0.0, 0.0),
          mm(0.0, 0.0, 0.0) {}

      /** Construct with values */
      PositionData(vec3<double> px_, vec3<double> mm_)
        : px(px_),
          mm(mm_) {}

      /** Construct with seperate values */
      PositionData(vec2<double> px_, double frame_,
                   vec2<double> mm_, double angle_)
        : px(px_[0], px_[1], frame_),
          mm(mm_[0], mm_[1], angle_) {}
    };

    MillerIndex miller_index;
    vec3<double> beam_vector;
    PositionData position;
    std::size_t panel;
    bool entering;

    /** Default construct */
    Prediction()
      : miller_index(0, 0, 0),
        beam_vector(0.0, 0.0, 0.0),
        panel(0),
        entering(false) {}

    /** Construct with miller index */
    Prediction(MillerIndex miller_index_)
      : miller_index(miller_index_),
        beam_vector(0.0, 0.0, 0.0),
        panel(0),
        entering(false) {}

    /**
     * Construct with values.
     */
    Prediction(MillerIndex miller_index_,
               vec3<double> beam_vector_,
               const PositionData &position_,
               std::size_t panel_,
               bool entering_)
      : miller_index(miller_index_),
        beam_vector(beam_vector_),
        position(position_),
        panel(panel_),
        entering(entering_) {}

    /**
     * Construct with values.
     */
    Prediction(MillerIndex miller_index_,
               vec3<double> beam_vector_,
               vec3<double> px_position_,
               vec3<double> mm_position_,
               std::size_t panel_,
               bool entering_)
      : miller_index(miller_index_),
        beam_vector(beam_vector_),
        position(px_position_, mm_position_),
        panel(panel_),
        entering(entering_) {}

    /**
     * Construct with values.
     */
    Prediction(MillerIndex miller_index_,
               vec3<double> beam_vector_,
               vec2<double> px_position_,
               double frame_,
               vec2<double> mm_position_,
               double angle_,
               std::size_t panel_,
               bool entering_)
      : miller_index(miller_index_),
        beam_vector(beam_vector_),
        position(px_position_, frame_, mm_position_, angle_),
        panel(panel_),
        entering(entering_) {}
  };

}}; // namespace dials::model

#endif /* DIALS_MODEL_DATA_PREDICTION_H */
