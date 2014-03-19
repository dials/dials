/*
 * calculate_reflection_frames.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_CALCULATE_REFLECTION_FRAMES_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_CALCULATE_REFLECTION_FRAMES_H

#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <dxtbx/model/scan.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace algorithms {

  using scitbx::vec2;
  using scitbx::constants::two_pi;
  using dxtbx::model::Scan;
  using dials::model::Reflection;

  /**
   * For the given reflection, calculate the list of frames on which the
   * reflection will be observed and return an array of reflections with all
   * values duplicated and the frame number updated.
   * @param scan The scan model object
   * @param reflection The reflection object
   * @returns A list of reflections
   */
  inline
  af::shared <Reflection> reflection_frames(const Scan &scan,
      const Reflection &reflection) {

    double phi = reflection.get_rotation_angle();

    // Get the frames that a reflection with this angle will be observed at
    af::shared< vec2<double> > frames = scan.get_array_indices_with_angle(phi);

    // Loop through all the frames and duplicate the reflection for each
    af::shared <Reflection> reflections_new;
    for (std::size_t j = 0; j < frames.size(); ++j) {
      Reflection r(reflection);
      r.set_rotation_angle(frames[j][0]);
      r.set_frame_number(frames[j][1]);
      reflections_new.push_back(r);
    }

    // Return the reflection list
    return reflections_new;
  }

  /**
   * For each reflection, calculate the list of frames on which the
   * reflection will be observed and return an array of reflections with all
   * values duplicated and the frame number updated.
   * @param scan The scan model object
   * @param reflections The reflection list
   * @returns A list of reflections
   */
  inline
  af::shared <Reflection> reflection_frames(const Scan &scan,
      const af::const_ref<Reflection> &reflections) {

    // Loop through all the given reflections
    af::shared <Reflection> reflections_new;
    for (std::size_t i = 0; i < reflections.size(); ++i) {

      // Get the duplicate reflections at different frames and for each, append
      // to the list of reflections.
      af::shared <Reflection> r = reflection_frames(scan, reflections[i]);
      for (std::size_t j = 0; j < r.size(); ++j) {
        reflections_new.push_back(r[j]);
      }
    }

    // Return the new list of reflectiond
    return reflections_new;
  }

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_CALCULATE_REFLECTION_FRAMES_H
