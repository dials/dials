
#ifndef DIALS_MODEL_EXPERIMENT_DETECTOR_HELPERS_H
#define DIALS_MODEL_EXPERIMENT_DETECTOR_HELPERS_H

#include "detector.h"

namespace dials { namespace model { namespace experiment {

  bool is_coordinate_valid(const FlatPanelDetector &detector, 
      vec2 <double> coord) {
    return (coord[0] >= 0 && coord[0] < detector.get_image_size()[0])
        && (coord[1] >= 0 && coord[1] < detector.get_image_size()[1]);
  }

  bool is_coordinate_valid(const MultiFlatPanelDetector &detector, 
      vec3 <double> coord) {
    int panel = (int)coord[0];
    return (coord[0] >= 0 && coord[0] < detector.get_num_panels())
        && (coord[1] >= 0 && coord[1] < detector[panel].get_image_size()[0])
        && (coord[2] >= 0 && coord[2] < detector[panel].get_image_size()[1]);
  }

}}} // namespace dials::model::experiment

#endif // DIALS_MODEL_EXPERIMENT_DETECTOR_HELPERS_H