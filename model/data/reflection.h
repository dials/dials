/*
 * reflection.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_DATA_REFLECTION_H
#define DIALS_MODEL_DATA_REFLECTION_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/shared.h>
#include <cctbx/miller.h>

namespace dials { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;

  typedef cctbx::miller::index <> miller_index_type;

  class ReflectionBase {
  public:

    ReflectionBase() {}

    ReflectionBase(miller_index_type miller_index)
      : miller_index_(miller_index) {}

    virtual ~ReflectionBase() {}

    miller_index_type get_miller_index() const {
      return miller_index_;
    }

    void set_miller_index(miller_index_type miller_index) {
      miller_index_ = miller_index;
    }

    bool is_zero() {
      return miller_index_.is_zero();
    }

  protected:
    miller_index_type miller_index_;
  };

  class Reflection : public ReflectionBase {
  public:

    Reflection() {}

    Reflection(miller_index_type miller_index)
      : ReflectionBase(miller_index) {}

    Reflection(miller_index_type miller_index,
               double rotation_angle,
               vec3 <double> beam_vector)
      : ReflectionBase(miller_index),
        rotation_angle_(rotation_angle),
        beam_vector_(beam_vector){}

    virtual ~Reflection() {}

    double get_rotation_angle() const {
      return rotation_angle_;
    }

    vec3 <double> get_beam_vector() const {
      return beam_vector_;
    }

    vec2 <double> get_image_coord_mm() const {
      return image_coord_mm_;
    }

    vec2 <double> get_image_coord_px() const {
      return image_coord_px_;
    }

    double get_frame_number() const {
      return frame_number_;
    }

    int get_panel_number() const {
      return panel_number_;
    }

    void set_rotation_angle(double rotation_angle) {
      rotation_angle_ = rotation_angle;
    }

    void set_beam_vector(vec3 <double> beam_vector) {
      beam_vector_ = beam_vector;
    }

    void set_image_coord_mm(vec2 <double> image_coord_mm) {
      image_coord_mm_ = image_coord_mm;
    }

    void set_image_coord_px(vec2 <double> image_coord_px) {
      image_coord_px_ = image_coord_px;
    }

    void set_frame_number(double frame_number) {
      frame_number_ = frame_number;
    }

    void set_panel_number(int panel_number) {
      panel_number_ = panel_number;
    }

  protected:

    double rotation_angle_;
    vec3 <double> beam_vector_;
    vec2 <double> image_coord_px_;
    vec2 <double> image_coord_mm_;
    double frame_number_;
    int panel_number_;
  };

  typedef scitbx::af::flex <Reflection>::type ReflectionList;

}} // namespace dials::model

#endif /* DIALS_MODEL_DATA_REFLECTION_H */
