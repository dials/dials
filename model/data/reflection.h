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
#include <scitbx/array_family/flex_types.h>
#include <cctbx/miller.h>

namespace dials { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;

  // Typedef the miller type
  typedef cctbx::miller::index <> miller_index_type;

  /** A base class for reflections. */
  class ReflectionBase {
  public:

    /** Initialise the reflection */
    ReflectionBase()
      : miller_index_(0, 0, 0) {}

    /**
     * Initialise the reflection with the miller index
     * @param miller_index The miller index
     */
    ReflectionBase(miller_index_type miller_index)
      : miller_index_(miller_index) {}

    /** Virtual destructor */
    virtual ~ReflectionBase() {}

    /** Get the miller index */
    miller_index_type get_miller_index() const {
      return miller_index_;
    }

    /** Set the miller index */
    void set_miller_index(miller_index_type miller_index) {
      miller_index_ = miller_index;
    }

    /** True/False the miller index is (0, 0, 0) */
    bool is_zero() {
      return miller_index_.is_zero();
    }

  protected:
    miller_index_type miller_index_;
  };

  /**
   * A class for reflections containing some additional data
   */
  class Reflection : public ReflectionBase {
  public:

    /** Default initialisation */
    Reflection()
      : ReflectionBase(),
        rotation_angle_(0.0),
        beam_vector_(0.0, 0.0, 0.0),
        image_coord_px_(0.0, 0.0),
        image_coord_mm_(0.0, 0.0),
        frame_number_(0),
        panel_number_(0),
        shoebox_(0, 0, 0, 0, 0, 0) {}

    /**
     * Initialise the reflection with the miller index
     * @param miller_index The miller index
     */
    Reflection(miller_index_type miller_index)
      : ReflectionBase(miller_index),
        rotation_angle_(0.0),
        beam_vector_(0.0, 0.0, 0.0),
        image_coord_px_(0.0, 0.0),
        image_coord_mm_(0.0, 0.0),
        frame_number_(0),
        panel_number_(0),
        shoebox_(0, 0, 0, 0, 0, 0) {}

    /**
     * Initialise the reflection with the miller index, rotation angle and
     * beam parameters
     * @param miller_index The miller index
     * @param rotation_angle The rotation angle
     * @param beam_vector The beam vector
     */
    Reflection(miller_index_type miller_index,
               double rotation_angle,
               vec3 <double> beam_vector)
      : ReflectionBase(miller_index),
        rotation_angle_(rotation_angle),
        beam_vector_(beam_vector),
        image_coord_px_(0.0, 0.0),
        image_coord_mm_(0.0, 0.0),
        frame_number_(0),
        panel_number_(0),
        shoebox_(0, 0, 0, 0, 0, 0) {}

    /** Virtual destructor */
    virtual ~Reflection() {}

    /** Get the rotation angle */
    double get_rotation_angle() const {
      return rotation_angle_;
    }

    /** Get the beam vector */
    vec3 <double> get_beam_vector() const {
      return beam_vector_;
    }

    /** Get the image coordinate in mm */
    vec2 <double> get_image_coord_mm() const {
      return image_coord_mm_;
    }

    /** Get the image coordinate in pixels */
    vec2 <double> get_image_coord_px() const {
      return image_coord_px_;
    }

    /** Get the frame number */
    double get_frame_number() const {
      return frame_number_;
    }

    /** Get the panel number */
    int get_panel_number() const {
      return panel_number_;
    }

    /** Get the shoebox */
    int6 get_shoebox() const {
      return shoebox_;
    }

    /** Set the rotation angle */
    void set_rotation_angle(double rotation_angle) {
      rotation_angle_ = rotation_angle;
    }

    /** Set the beam vector */
    void set_beam_vector(vec3 <double> beam_vector) {
      beam_vector_ = beam_vector;
    }

    /** Set the image coordinate in mm */
    void set_image_coord_mm(vec2 <double> image_coord_mm) {
      image_coord_mm_ = image_coord_mm;
    }

    /** Set the image coordinate in pixels */
    void set_image_coord_px(vec2 <double> image_coord_px) {
      image_coord_px_ = image_coord_px;
    }

    /** Set the frame number */
    void set_frame_number(double frame_number) {
      frame_number_ = frame_number;
    }

    /** Set the panel number */
    void set_panel_number(int panel_number) {
      panel_number_ = panel_number;
    }

    /** Set the shoebox */
    void set_shoebox(int6 shoebox) {
      shoebox_ = shoebox;
    }

  protected:

    double rotation_angle_;
    vec3 <double> beam_vector_;
    vec2 <double> image_coord_px_;
    vec2 <double> image_coord_mm_;
    double frame_number_;
    int panel_number_;
    int6 shoebox_;
  };

  typedef scitbx::af::flex <Reflection>::type ReflectionList;

}} // namespace dials::model

#endif /* DIALS_MODEL_DATA_REFLECTION_H */
