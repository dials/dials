/*
 * reflection_predictor.h
 *
 *  Copyright (C) 2015 Diamond Light Source, CCP4
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef IDY_ALGORITHMS_SPOT_PREDICTION_REFLECTION_PREDICTOR_H
#define IDY_ALGORITHMS_SPOT_PREDICTION_REFLECTION_PREDICTOR_H

#include <algorithm>
#include <scitbx/constants.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dials/algorithms/spot_prediction/reeke_index_generator.h>
#include <dials/algorithms/spot_prediction/ray_predictor.h>
#include <dials/scratch/idy/algorithms/spot_prediction/stills_ray_predictor.h>
#include <dials/algorithms/spot_prediction/reflection_predictor.h>

namespace dials { namespace algorithms {

  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dials::model::Ray;

  class StillsExperimentalReflectionPredictor : public StillsDeltaPsiReflectionPredictor {
  public:

    StillsExperimentalReflectionPredictor(
        const Beam &beam,
        const Detector &detector,
        mat3<double> ub,
        const cctbx::uctbx::unit_cell &unit_cell,
        const cctbx::sgtbx::space_group_type &space_group_type,
        const double &dmin,
        double const& mos)
      : StillsDeltaPsiReflectionPredictor(beam,detector,ub,unit_cell,space_group_type,dmin),
        half_mosaicity_rad_(mos),
        ci_predict_ray_(beam.get_s0(),mos) {}

    virtual void append_for_index(stills_prediction_data &p,
        const mat3<double> ub,
        const miller_index &h, int panel=-1) {
      Ray ray;
      ray = ci_predict_ray_(h, ub);
      double delpsi = ci_predict_ray_.get_delpsi();
      append_for_ray(p, h, ray, panel, delpsi);
    }

  private:
    const double half_mosaicity_rad_;
    StillsCentralImpactRayPredictor ci_predict_ray_;
  };

}} // namespace dials::algorithms

#endif // IDY_ALGORITHMS_SPOT_PREDICTION_REFLECTION_PREDICTOR_H
