/*
 * detector.h
 *
 *   Copyright (C) 2013 Diamond Light Source, James Parkhurst
 *
 *   This code is distributed under the BSD license, a copy of which is
 *   included in the root directory of this package.
 */
#ifndef DIALS_MODEL_EXPERIMENT_DETECTOR_H
#define DIALS_MODEL_EXPERIMENT_DETECTOR_H

#include <string>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/shared.h>

namespace dials { namespace model { namespace experiment {

  using scitbx::vec2;
  using scitbx::vec3;

  /** A base class for detectors */
  class DetectorBase {};

  /**
   * A class representing a detector panel. A detector can have multiple
   * panels which are each represented by this class.
   *
   * The class contains the following members accessible through getter and
   * setter methods:
   *   type A string representing the type of the detector panel (i.e the
   *     manufactors name or some such identifier).
   *
   *   x_axis A unit vector pointing along the x (fast) axis of the panel. The
   *     vector is given with respect to the laboratory coordinate frame.
   *
   *   y_axis A unit vector pointing along the y (slow) axis of the panel. The
   *     vector is given with respect to the laboratory coordinate frame.
   *
   *   normal A unit vector giving the normal to the panel plane. The
   *     vector is given with respect to the laboratory coordinate frame.
   *
   *   origin The origin of the detector plane. (i.e. the laboratory coordinate
   *     of the edge of the zeroth pixel
   *
   *   pixel_size The size of the pixels in mm. The convention used is that
   *     of (y, x) i.e. (slow, fast)
   *
   *   image_size The size of the panel in pixels. The convention used is that
   *     of (y, x) i.e. (slow, fast)
   *
   *   trusted_range The range of counts that are considered reliable, given
   *     in the range [min, max].
   *
   *   distance The signed distance from the source to the detector
   *
   * In the document detailing the conventions used:
   *    type -> *unspecified*
   *    x_axis -> d1
   *    y_axis -> d2
   *    normal -> d3
   *    origin -> *unspecified*
   *    pixel_size -> *unspecified*
   *    image_size -> *unspecified*
   *    trusted_range -> *unspecified*
   *    distance -> *unspecified*
   */
  class FlatPanelDetector : public DetectorBase {
  public:

    /** The default constructor */
    FlatPanelDetector()
      : type_("Unknown"),
        x_axis_(1.0, 0.0, 0.0),
        y_axis_(0.0, 1.0, 0.0),
        normal_(0.0, 0.0, 1.0),
        origin_(0.0, 0.0, 0.0),
        pixel_size_(0.0, 0.0),
        image_size_(0, 0),
        distance_(0.0) {}

    /**
     * Initialise the detector panel.
     * @param type The type of the detector panel
     * @param x_axis The x axis of the detector. The given vector is normalized.
     * @param y_axis The y axis of the detector. The given vector is normalized.
     * @param normal The detector normal. The given vector is normalized.
     * @param origin The detector origin
     * @param pixel_size The size of the individual pixels
     * @param image_size The size of the detector panel (in pixels)
     * @param distance The distance from the detector to the crystal origin
     */
    FlatPanelDetector(std::string type,
                      vec3 <double> x_axis,
                      vec3 <double> y_axis,
                      vec3 <double> normal,
                      vec3 <double> origin,
                      vec2 <double> pixel_size,
                      vec2 <std::size_t> image_size,
                      double distance)
      : type_(type),
        x_axis_(x_axis.normalize()),
        y_axis_(y_axis.normalize()),
        normal_(normal.normalize()),
        origin_(origin),
        pixel_size_(pixel_size),
        image_size_(image_size),
        distance_(distance) {}

    /** Get the sensor type */
    std::string get_type() const {
      return type_;
    }

    /** Get the x axis */
    vec3 <double> get_x_axis() const {
      return x_axis_;
    }

    /** Get the y axis */
    vec3 <double> get_y_axis() const {
      return y_axis_;
    }

    /** Get the normal */
    vec3 <double> get_normal() const {
      return normal_;
    }

    /** Get the pixel origin */
    vec3 <double> get_origin() const {
      return origin_;
    }

    /** Get the pixel size */
    vec2 <double> get_pixel_size() const {
      return pixel_size_;
    }

    /** Get the image size */
    vec2 <std::size_t> get_image_size() const {
      return image_size_;
    }

    /** Get the distance from the crystal */
    double get_distance() const {
      return distance_;
    }

    /** Set the detector panel type */
    void set_type(std::string type) {
      type_ = type;
    }

    /** Set the x axis */
    void set_x_axis(vec3 <double> x_axis) {
      x_axis_ = x_axis;
    }

    /** Set the y axis */
    void set_y_axis(vec3 <double> y_axis) {
      y_axis_ = y_axis;
    }

    /** Set the normal */
    void set_normal(vec3 <double> normal) {
      normal_ = normal;
    }

    /** Set the origin */
    void set_origin(vec3 <double> origin) {
      origin_ = origin;
    }

    /** Set the pixel size */
    void set_pixel_size(vec2 <double> pixel_size) {
      pixel_size_ = pixel_size;
    }

    /** Set the image size */
    void set_image_size(vec2 <std::size_t> image_size) {
      image_size_ = image_size;
    }

    /* Set the distance from the crystal */
    void set_distance(double distance) {
      distance_ = distance;
    }

  protected:

    std::string type_;
    vec3 <double> x_axis_;
    vec3 <double> y_axis_;
    vec3 <double> normal_;
    vec3 <double> origin_;
    vec2 <double> pixel_size_;
    vec2 <std::size_t> image_size_;
    double distance_;
  };

  class MultiFlatPanelDetector : public DetectorBase {
  public:
    
    typedef FlatPanelDetector panel_type;
    typedef scitbx::af::shared <panel_type> panel_list_type;
    typedef panel_list_type::iterator iterator;

    MultiFlatPanelDetector()
      : type_("Unknown") {}

    MultiFlatPanelDetector(std::string type)
      :type_(type) {}

    iterator begin() {
      return panel_list_.begin();
    }

    iterator end() {
      return panel_list_.end();
    }

    void add_panel(const panel_type &panel) {
      panel_list_.push_back(panel);
    }

    void remove_panels() {
      panel_list_.erase(panel_list_.begin(), panel_list_.end());
    }

    void remove_panel(std::size_t i) {
      panel_list_.erase(panel_list_.begin() + i);
    }

    std::size_t num_panels() const {
      return panel_list_.size();
    }

    panel_type& operator[](std::size_t index) {
      return panel_list_[index];
    }

    const panel_type& operator[](std::size_t index) const {
      return panel_list_[index];
    }

  protected:

    std::string type_;
    panel_list_type panel_list_;
  };

}}} // namespace dials::model::experiment

#endif // DIALS_MODEL_EXPERIMENT_DETECTOR_H
