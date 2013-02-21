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

  class Reflection {
  public:

    Reflection()
      : miller_index_(0, 0, 0),
        rotation_angle_(0.0),
        beam_vector_(0.0, 0.0, 0.0),
        image_coord_(0.0, 0.0) {}

    Reflection(cctbx::miller::index <> miller_index,
               double rotation_angle,
               vec3 <double> beam_vector,
               vec2 <double> image_coord)
      : miller_index_(miller_index),
        rotation_angle_(rotation_angle),
        beam_vector_(beam_vector),
        image_coord_(image_coord) {}

  public:

    cctbx::miller::index <> get_miller_index() const {
      return miller_index_;
    }

    double get_rotation_angle() const {
      return rotation_angle_;
    }

    vec3 <double> get_beam_vector() const {
      return beam_vector_;
    }

    vec2 <double> get_image_coord() const {
      return image_coord_;
    }

    void set_miller_index(cctbx::miller::index <> miller_index) {
      miller_index_ = miller_index;
    }

    void set_rotation_angle(double rotation_angle) {
      rotation_angle_ = rotation_angle;
    }

    void set_beam_vector(vec3 <double> beam_vector) {
      beam_vector_ = beam_vector;
    }

    void set_image_coord(vec2 <double>  image_coord) {
      image_coord_ = image_coord;
    }

    bool is_zero() {
      return miller_index_.is_zero();
    }

  private:

    cctbx::miller::index <> miller_index_;
    double                  rotation_angle_;
    vec3 <double>   beam_vector_;
    vec2 <double>   image_coord_;
  };

  typedef scitbx::af::flex <Reflection>::type ReflectionList;

}} // namespace dials::model

#endif /* DIALS_MODEL_DATA_REFLECTION_H */
