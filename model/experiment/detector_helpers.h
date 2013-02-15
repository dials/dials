/*
 * detector_helpers.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_EXPERIMENT_DETECTOR_HELPERS_H
#define DIALS_MODEL_EXPERIMENT_DETECTOR_HELPERS_H

#include "detector.h"

namespace dials { namespace model {

  bool is_coordinate_valid(const FlatPanelDetector &detector,
      vec2 <double> coord) {
    return (coord[0] >= 0 && coord[0] < detector.get_image_size()[0])
        && (coord[1] >= 0 && coord[1] < detector.get_image_size()[1]);
  }

//  vec3 <double>
//  pixel_to_mm(const FlatPanelDetector &detector, vec2 <double> xy) {
//    return detector.get_d_matrix() * vec3 <double> (xy[0], xy[1], 1);
//  }

//  vec2 <double>
//  mm_to_pixel(const FlatPanelDetector &detector, vec3 <double> beam_vector) {
//    vec3 <double> v = detector.get_inverse_d_matrix() * beam_vector;
//    return vec2 <double> (v[0] / v[2], v[1] / v[2]);
//  }


  bool is_coordinate_valid(const MultiFlatPanelDetector &detector,
      vec3 <double> coord) {
    int panel = (int)coord[0];
    return (coord[0] >= 0 && coord[0] < detector.num_panels())
        && (coord[1] >= 0 && coord[1] < detector[panel].get_image_size()[0])
        && (coord[2] >= 0 && coord[2] < detector[panel].get_image_size()[1]);
  }

}} // namespace dials::model

#endif // DIALS_MODEL_EXPERIMENT_DETECTOR_HELPERS_H
