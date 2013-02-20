/*
 * scan.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_EXPERIMENT_SCAN_H
#define DIALS_MODEL_EXPERIMENT_SCAN_H

#include <scitbx/vec2.h>
#include <dials/error.h>

namespace dials { namespace model {

  using scitbx::vec2;

  /** A scan base class */
  class ScanBase {};

  /** A class to represent a scan */
  class Scan : public ScanBase {
  public:

    /** The default constructor */
    Scan() {}

    /**
     * Initialise the class
     * @param image_range The range of images covered by the scan
     * @param starting_angle The starting rotation angle
     * @param oscillation_range The oscillation range of each frame
     */
    Scan(vec2 <int> image_range,
         double starting_angle,
         double oscillation_range)
      : image_range_(image_range),
        starting_angle_(starting_angle),
        oscillation_range_(oscillation_range),
        num_images_(image_range_[1] - image_range_[0]) {
      DIALS_ASSERT(num_images_ >= 0);
    }

    /** Virtual destructor */
    virtual ~Scan() {}

    /** Get the image range */
    vec2 <int> get_image_range() const {
      return image_range_;
    }

    /** Get the starting angle */
    double get_starting_angle() const {
      return starting_angle_;
    }

    /** Get the oscillation range */
    double get_oscillation_range() const {
      return oscillation_range_;
    }
    
    /** Get the number of images */
    int get_num_images() const {
      return num_images_;
    }

    /** Get the total oscillation range of the scan */
    double get_total_oscillation_range() const {
      return num_images_ * oscillation_range_;
    }
    
    /** Set the image range */
    void set_image_range(vec2 <int> image_range) {
      image_range_ = image_range;
      num_images_ = image_range_[1] - image_range_[0];
      DIALS_ASSERT(num_images_ >= 0);
    }

    /** Set the starting angle */
    void set_starting_angle(double starting_angle) {
      starting_angle_ = starting_angle;
    }

    /** Set the oscillation range */
    void set_oscillation_range(double oscillation_range) {
      oscillation_range_ = oscillation_range;
    }
    
    /** Check the scans are the same */
    bool operator==(const Scan &scan) {
      double eps = 1.0e-6;
      double d_angle = std::abs(starting_angle_ - scan.starting_angle_);
      double d_range = std::abs(oscillation_range_ - scan.oscillation_range_);
      return image_range_ == scan.image_range_ &&
             d_angle <= eps && d_range <= eps;
    }
    
    /** Check the scans are not the same */
    bool operator!=(const Scan &scan) {
      return !(*this == scan);
    }
    
  private:

    vec2 <int> image_range_;
    double starting_angle_;
    double oscillation_range_;
    int num_images_;
  };

}} // namespace dials::model

#endif // DIALS_MODEL_EXPERIMENT_SCAN_H
