/*
 * reflection_predictor.h
 *
 *  Copyright (C) 2013 Diamond Light Source, CCP4
 *
 *  Author: David Waterman
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_REFLECTION_PREDICTOR_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_REFLECTION_PREDICTOR_H

#include <scitbx/math/r3_rotation.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dials/array_family/reflection_table.h>
#include <dials/algorithms/spot_prediction/index_generator.h>
#include <dials/algorithms/spot_prediction/reeke_index_generator.h>
#include <dials/algorithms/spot_prediction/ray_predictor.h>
#include <dials/algorithms/spot_prediction/scan_varying_ray_predictor.h>
#include <dials/algorithms/spot_prediction/stills_ray_predictor.h>
#include <dials/algorithms/spot_prediction/ray_intersection.h>

namespace dials { namespace algorithms {

  using boost::shared_ptr;
  using scitbx::math::r3_rotation::axis_and_angle_as_matrix;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using dials::model::Ray;

  /**
   * Helper struct for holding prediction data internally/
   */
  struct prediction_data {

    af::shared< miller_index > hkl;
    af::shared< std::size_t > panel;
    af::shared< bool > enter;
    af::shared< vec3<double> > s1;
    af::shared< vec3<double> > xyz_px;
    af::shared< vec3<double> > xyz_mm;

    prediction_data(af::reflection_table &table) {
      hkl    = table.get< miller_index >("miller_index");
      panel  = table.get< std::size_t >("panel");
      enter  = table.get< bool >("entering");
      s1     = table.get< vec3<double> >("s1");
      xyz_px = table.get< vec3<double> >("xyzcal.px");
      xyz_mm = table.get< vec3<double> >("xyzcal.mm");
    }
  };

  /**
   * A reflection predictor for scan static prediction.
   */
  class ScanStaticReflectionPredictor {

    typedef cctbx::miller::index<> miller_index;

  public:

    /**
     * Keep a reference to all the models.
     */
    ScanStaticReflectionPredictor(
          const Beam &beam,
          const Detector &detector,
          const Goniometer &goniometer,
          const Scan &scan,
          const cctbx::uctbx::unit_cell &unit_cell,
          const cctbx::sgtbx::space_group_type &space_group_type,
          mat3<double> ub,
          double dmin)
      : beam_(beam),
        detector_(detector),
        goniometer_(goniometer),
        scan_(scan),
        unit_cell_(unit_cell),
        space_group_type_(space_group_type),
        ub_(ub),
        dmin_(dmin),
        predict_rays_(
            beam.get_s0(),
            goniometer.get_rotation_axis(),
            scan.get_oscillation_range()){}

    /**
     * Predict all the reflections.
     * @returns A reflection table
     */
    af::reflection_table operator()() const {

      // Create the reflection table and the local container
      af::reflection_table table;
      prediction_data predictions(table);

      // Create the index generate and loop through the indices. For each index,
      // predict the rays and append to the reflection table
      IndexGenerator indices(unit_cell_, space_group_type_, dmin_);
      for (;;) {
        miller_index h = indices.next();
        if (h.is_zero()) {
          break;
        }
        append_for_index(predictions, h);
      }

      // Return the reflection table
      return table;
    }

    /**
     * Predict the reflections with the given miller indices.
     * @param h The miller indices.
     * @returns The reflection table.
     */
    af::reflection_table operator()(
        const af::const_ref< miller_index > &h) const {

      // Create the reflection table and local container
      af::reflection_table table;
      prediction_data predictions(table);

      // Create the ray predictor and loop through all the reflection
      // and predict the rays.
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predictions, h[i]);
      }

      // Return the reflection table
      return table;
    }

    /**
     * Predict the reflection on a particular panel.
     * @param h The miller index.
     * @param panel The panel number
     * @returns The reflection table.
     */
    af::reflection_table operator()(
        const af::const_ref< miller_index > &h,
        std::size_t panel) const {
      af::shared<std::size_t> panels(h.size(), panel);
      return (*this)(h, panels.const_ref());
    }

    /**
     * Predict the reflection on a particular panel.
     * @param h The miller index.
     * @param panel The panel number
     * @returns The reflection table.
     */
    af::reflection_table operator()(
        const af::const_ref< miller_index > &h,
        const af::const_ref< std::size_t > &panel) const {
      DIALS_ASSERT(h.size() == panel.size());

      // Create the reflection table and local container
      af::reflection_table table;
      prediction_data predictions(table);

      // Create the ray predictor and loop through all the reflection
      // and predict the rays.
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predictions, h[i], (int)panel[i]);
      }

      // Return the reflection table
      return table;
    }

  private:

    /**
     * For the given index, predict the rays and append the reflections.
     */
    void append_for_index(prediction_data &p,
        const miller_index &h, int panel=-1) const {
      af::small<Ray, 2> rays = predict_rays_(h, ub_);
      for (std::size_t i = 0; i < rays.size(); ++i) {
        DIALS_ASSERT(panel < (int)detector_.size());
        append_for_ray(p, h, rays[i], panel);
      }
    }

    /**
     * For each ray, get the impact and append to the reflection table.
     */
    void append_for_ray(prediction_data &p, const miller_index &h,
        const Ray &ray, int panel) const {
      try {

        // Get the impact on the detector
        Detector::coord_type impact = get_ray_intersection(ray.s1, panel);
        std::size_t pid = impact.first;
        vec2<double> mm = impact.second;
        vec2<double> px = detector_[pid].millimeter_to_pixel(mm);

        // Get the frames that a reflection with this angle will be observed at
        af::shared<double> frames =
          scan_.get_array_indices_with_angle(ray.angle);
        for (std::size_t i = 0; i < frames.size(); ++i) {
          p.hkl.push_back(h);
          p.enter.push_back(ray.entering);
          p.s1.push_back(ray.s1);
          p.xyz_mm.push_back(vec3<double>(mm[0], mm[1], ray.angle));
          p.xyz_px.push_back(vec3<double>(px[0], px[1], frames[i]));
          p.panel.push_back(pid);
        }

      } catch(dxtbx::error) {
        // do nothing
      }
    }

    /**
     * Helper function to do ray intersection with/without panel set.
     */
    Detector::coord_type get_ray_intersection(vec3<double> s1, int panel) const {
      Detector::coord_type coord;
      if (panel < 0) {
        coord = detector_.get_ray_intersection(s1);
      } else {
        coord.first = panel;
        coord.second = detector_[panel].get_ray_intersection(s1);
      }
      return coord;
    }

    Beam beam_;
    Detector detector_;
    Goniometer goniometer_;
    Scan scan_;
    cctbx::uctbx::unit_cell unit_cell_;
    cctbx::sgtbx::space_group_type space_group_type_;
    mat3<double> ub_;
    double dmin_;
    RayPredictor2 predict_rays_;
  };



  class ScanVaryingReflectionPredictor {

    typedef cctbx::miller::index<> miller_index;

  public:

    ScanVaryingReflectionPredictor(
        shared_ptr<Beam> beam,
        shared_ptr<Detector> detector,
        shared_ptr<Goniometer> goniometer,
        shared_ptr<Scan> scan)
      : beam_(beam),
        detector_(detector),
        goniometer_(goniometer),
        scan_(scan) {}

    af::reflection_table all_observable() const {
      af::reflection_table table;
      prediction_data predictions(table);
      ScanVaryingRayPredictor predict_rays = init_ray_predictor();
      vec2<int> array_range = scan_->get_array_range();
      for (int frame = array_range[0]; frame < array_range[1]; ++frame) {
        append_for_image(predict_rays, predictions, frame);
      }
      return table;
    }

  private:

    ScanVaryingRayPredictor init_ray_predictor() const {
      return ScanVaryingRayPredictor(
          beam_->get_s0(),
          goniometer_->get_rotation_axis(),
          scan_->get_oscillation(),
          dmin_);
    }

    void append_for_image(const ScanVaryingRayPredictor &predict_rays,
        prediction_data &p, int frame) const {

      vec3<double> m2 = goniometer_->get_rotation_axis();
      vec3<double> s0 = beam_->get_s0();

      int frame_0 = scan_->get_array_range()[0];
      double phi_beg = scan_->get_angle_from_array_index(frame);
      double phi_end = scan_->get_angle_from_array_index(frame + 1);
      mat3<double> r_beg = axis_and_angle_as_matrix(m2, phi_beg);
      mat3<double> r_end = axis_and_angle_as_matrix(m2, phi_end);
      mat3<double> A1 = r_beg * A_[frame - frame_0];
      mat3<double> A2 = r_end * A_[frame - frame_0 + 1];
      ReekeIndexGenerator indices(A1, A2, m2, s0, dmin_, margin_);
      for (;;) {
        miller_index h = indices.next();
        if (h.is_zero()) {
          break;
        }
        append_for_index(predict_rays, p, A1, A2, frame, h);
      }
    }

    void append_for_index(const ScanVaryingRayPredictor &predict_rays,
        prediction_data &p, mat3<double> A1, mat3<double> A2, std::size_t frame,
        const miller_index &h) const {
      boost::optional<Ray> ray = predict_rays(h, A1, A2, frame, 1);
      if (ray) {
        append_for_ray(p, h, *ray);
      }
    }

    void append_for_ray(prediction_data &p, const miller_index &h, const Ray &ray) const {
      try {

        // Get the impact on the detector
        Detector::coord_type impact = (*detector_).get_ray_intersection(ray.s1);
        std::size_t panel = impact.first;
        vec2<double> mm = impact.second;
        vec2<double> px = (*detector_)[panel].millimeter_to_pixel(mm);

        double frame = scan_->get_array_index_from_angle(ray.angle);

        // Get the frames that a reflection with this angle will be observed at
        p.hkl.push_back(h);
        p.enter.push_back(ray.entering);
        p.s1.push_back(ray.s1);
        p.xyz_mm.push_back(vec3<double>(mm[0], mm[1], ray.angle));
        p.xyz_px.push_back(vec3<double>(px[0], px[1], frame));
        p.panel.push_back(panel);

      } catch(dxtbx::error) {
        // do nothing
      }
    }

    cctbx::uctbx::unit_cell unit_cell_;
    cctbx::sgtbx::space_group_type space_group_type_;
    af::shared< mat3<double> > A_;
    double dmin_;
    shared_ptr<Beam> beam_;
    shared_ptr<Detector> detector_;
    shared_ptr<Goniometer> goniometer_;
    shared_ptr<Scan> scan_;
    std::size_t margin_;
  };


  class StillsReflectionPredictor {

    typedef cctbx::miller::index<> miller_index;

  public:

    StillsReflectionPredictor(
        shared_ptr<Beam> beam,
        shared_ptr<Detector> detector,
        mat3<double> ub)
      : beam_(beam),
        detector_(detector),
        ub_(ub) { }

    af::reflection_table all_observable() const {
      DIALS_ERROR("Not implemented");
      return af::reflection_table();
    }

    af::reflection_table observed(
        const af::const_ref< miller_index > &h,
        const af::const_ref< std::size_t > id) const {
      af::reflection_table table;
      prediction_data predictions(table);
      StillsRayPredictor predict_rays = init_ray_predictor();
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predict_rays, predictions, h[i]);
      }
      return table;
    }

  private:

    StillsRayPredictor init_ray_predictor() const {
      return StillsRayPredictor(beam_->get_s0());
    }

    void append_for_index(const StillsRayPredictor &predict_rays,
        prediction_data &p, const miller_index &h) const {
      af::small<Ray, 2> rays = predict_rays(h, ub_);
      for (std::size_t i = 0; i < rays.size(); ++i) {
        append_for_ray(p, h, rays[i]);
      }
    }

    void append_for_ray(prediction_data &p, const miller_index &h, const Ray &ray) const {
      try {

        // Get the impact on the detector
        Detector::coord_type impact = (*detector_).get_ray_intersection(ray.s1);
        std::size_t panel = impact.first;
        vec2<double> mm = impact.second;
        vec2<double> px = (*detector_)[panel].millimeter_to_pixel(mm);

        // Add the reflections to the table
        p.hkl.push_back(h);
        p.enter.push_back(ray.entering);
        p.s1.push_back(ray.s1);
        p.xyz_mm.push_back(vec3<double>(mm[0], mm[1], 0.0));
        p.xyz_px.push_back(vec3<double>(px[0], px[1], 0.0));
        p.panel.push_back(panel);

      } catch(dxtbx::error) {
        // do nothing
      }
    }

    shared_ptr<Beam> beam_;
    shared_ptr<Detector> detector_;
    mat3<double> ub_;
  };
}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_REFLECTION_PREDICTOR_H
