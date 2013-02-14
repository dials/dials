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

namespace dials { namespace model { namespace experiment {

/** A scan base class */
class ScanBase {};

/** A class to represent a scan */
class Scan : public ScanBase {
public:

  /** The default constructor */
  Scan() {}

  /**
   * Initialise the class
   * @param starting_frame The starting frame number
   * @param starting_angle The starting rotation angle
   * @param oscillation_range The oscillation range of each frame
   * @param num_frames The number of frames in the scan.
   */
  Scan(int starting_frame,
       double starting_angle,
       double oscillation_range,
       int num_frames)
    : starting_frame_(starting_frame),
      starting_angle_(starting_angle),
      oscillation_range_(oscillation_range),
      num_frames_(num_frames) {}

  /** Get the starting frame */
  int get_starting_frame() const {
    return starting_frame_;
  }

  /** Get the starting angle */
  double get_starting_angle() const {
    return starting_angle_;
  }

  /** Get the oscillation range */
  double get_oscillation_range() const {
    return oscillation_range_;
  }

  /** Get the number of frames */
  int get_num_frames() const {
    return num_frames_;
  }

  /** Get the total oscillation range of the scan */
  double get_total_oscillation_range() const {
    return num_frames_ * oscillation_range_;
  }

  /** Set the starting frame */
  void set_starting_frame(int starting_frame) {
    starting_frame_ = starting_frame;
  }

  /** Set the starting angle */
  void set_starting_angle(double starting_angle) {
    starting_angle_ = starting_angle;
  }

  /** Set the oscillation range */
  void set_oscillation_range(double oscillation_range) {
    oscillation_range_ = oscillation_range;
  }

  /** Set the number of frames */
  void set_num_frames(int num_frames) {
    num_frames_ = num_frames;
  }

private:

  int starting_frame_;
  double starting_angle_;
  double oscillation_range_;
  int num_frames_;
};

}}} // namespace dials::model::experiment

#endif // DIALS_MODEL_EXPERIMENT_SCAN_H
