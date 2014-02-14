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

#include <algorithm>
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
          const mat3<double>& ub,
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
     * Predict reflections for UB
     * @param ub The UB matrix
     * @returns A reflection table.
     */
    af::reflection_table operator()() const {
      return (*this)(ub_);
    }

    /**
     * Predict reflections for UB
     * @param ub The UB matrix
     * @returns A reflection table.
     */
    af::reflection_table operator()(const mat3<double> &ub) const {

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
        append_for_index(predictions, ub, h);
      }

      // Return the reflection table
      return table;
    }

    /**
     * Predict reflections for UB
     * @param ub The array of UB matrices
     * @param h The array of miller indices
     * @returns A reflection table.
     */
    af::reflection_table operator()(
        const af::const_ref< mat3<double> > &ub,
        const af::const_ref<miller_index> &h) const {
      DIALS_ASSERT(ub.size() == h.size());
      af::reflection_table table;
      prediction_data predictions(table);
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predictions, ub[i], h[i]);
      }
      return table;
    }


    /**
     * Predict reflections for UB
     * @param ub The array of UB matrices
     * @param h The array of miller indices
     * @param panel The array of panels
     * @returns A reflection table.
     */
    af::reflection_table operator()(
        const af::const_ref< mat3<double> >&ub,
        const af::const_ref<miller_index> &h,
        const af::const_ref<std::size_t> &panel) const {
      DIALS_ASSERT(ub.size() == h.size());
      DIALS_ASSERT(ub.size() == panel.size());
      af::reflection_table table;
      prediction_data predictions(table);
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predictions, ub[i], h[i], panel[i]);
      }
      return table;
    }

    /**
     * Predict reflections for UB
     * @param ub The array of UB matrices
     * @param h The array of miller indices
     * @param panel The array of panels
     * @param entering The array of entering flags
     * @returns A reflection table.
     */
    af::reflection_table operator()(
        const af::const_ref< mat3<double> >&ub,
        const af::const_ref<miller_index> &h,
        const af::const_ref<std::size_t> &panel,
        const af::const_ref<bool> &entering) const {
      DIALS_ASSERT(ub.size() == h.size());
      DIALS_ASSERT(ub.size() == panel.size());
      DIALS_ASSERT(ub.size() == entering.size());
      af::reflection_table table;
      prediction_data predictions(table);
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predictions, ub[i], h[i], panel[i], entering[i]);
      }
      return table;
    }

  private:

    void append_for_index(
        prediction_data &p,
        const mat3<double> ub,
        const miller_index &h) const {
      af::small<Ray, 2> rays = predict_rays_(h, ub);
      for (std::size_t i = 0; i < rays.size(); ++i) {
        append_for_ray(p, h, rays[i]);
      }
    }

    void append_for_index(
        prediction_data &p,
        const mat3<double> ub,
        const miller_index &h,
        std::size_t panel) const {
      af::small<Ray, 2> rays = predict_rays_(h, ub);
      for (std::size_t i = 0; i < rays.size(); ++i) {
        append_for_ray(p, h, rays[i], panel);
      }
    }

    void append_for_index(
        prediction_data &p,
        const mat3<double> ub,
        const miller_index &h,
        std::size_t panel,
        bool entering) const {
      af::small<Ray, 2> rays = predict_rays_(h, ub);
      for (std::size_t i = 0; i < rays.size(); ++i) {
        if (rays[i].entering == entering) {
          append_for_ray(p, h, rays[i], panel);
        }
      }
    }

    /**
     * For each ray, get the impact and append to the reflection table.
     */
    void append_for_ray(prediction_data &p, const miller_index &h,
        const Ray &ray) const {
      try {
        Detector::coord_type impact = detector_.get_ray_intersection(ray.s1);
        std::size_t panel = impact.first;
        vec2<double> mm = impact.second;
        vec2<double> px = detector_[panel].millimeter_to_pixel(mm);
        append_for_impact(p, h, ray, panel, mm, px);
      } catch(dxtbx::error) {
        // do nothing
      }
    }

    /**
     * For each ray, get the impact and append to the reflection table.
     */
    void append_for_ray(prediction_data &p, const miller_index &h,
        const Ray &ray, std::size_t panel) const {
      try {
        vec2<double> mm = detector_[panel].get_ray_intersection(ray.s1);
        vec2<double> px = detector_[panel].millimeter_to_pixel(mm);
        append_for_impact(p, h, ray, panel, mm, px);
      } catch(dxtbx::error) {
        // do nothing
      }
    }

   void append_for_impact(prediction_data &p,
       const miller_index &h,
       const Ray &ray,
       std::size_t panel,
       vec2<double> mm,
       vec2<double> px) const {
      af::shared<double> frames =
        scan_.get_array_indices_with_angle(ray.angle);
      for (std::size_t i = 0; i < frames.size(); ++i) {
        p.hkl.push_back(h);
        p.enter.push_back(ray.entering);
        p.s1.push_back(ray.s1);
        p.xyz_mm.push_back(vec3<double>(mm[0], mm[1], ray.angle));
        p.xyz_px.push_back(vec3<double>(px[0], px[1], frames[i]));
        p.panel.push_back(panel);
      }
   }

    Beam beam_;
    Detector detector_;
    Goniometer goniometer_;
    Scan scan_;
    cctbx::uctbx::unit_cell unit_cell_;
    cctbx::sgtbx::space_group_type space_group_type_;
    mat3<double> ub_;
    double dmin_;
    ScanStaticRayPredictor predict_rays_;
  };


  /**
   * A reflection predictor for scan varying prediction.
   */
  class ScanVaryingReflectionPredictor {

    typedef cctbx::miller::index<> miller_index;

    /** Struct to help with sorting. */
    struct sort_by_frame {
      af::const_ref<int> frame;
      sort_by_frame(af::const_ref<int> frame_)
        : frame(frame_) {}

      template <typename T>
      bool operator()(T a, T b) {
        return frame[a] < frame[b];
      }
    };

  public:

    /**
     * Initialise the predictor
     */
    ScanVaryingReflectionPredictor(
        const Beam &beam,
        const Detector &detector,
        const Goniometer &goniometer,
        const Scan &scan,
        const af::const_ref< mat3<double> > A,
        double dmin,
        std::size_t margin)
      : beam_(beam),
        detector_(detector),
        goniometer_(goniometer),
        scan_(scan),
        A_(A.begin(), A.end()),
        dmin_(dmin),
        margin_(margin),
        predict_rays_(
            beam.get_s0(),
            goniometer.get_rotation_axis(),
            scan.get_oscillation(),
            dmin) {
      DIALS_ASSERT(A_.size() == scan_.get_num_images() + 1);
    }

    /**
     * Predict all the reflections for this model.
     * @returns The reflection table
     */
    af::reflection_table operator()() const {

      // Create the table and local stuff
      af::reflection_table table;
      prediction_data predictions(table);

      // Get the array range and loop through all the images
      vec2<int> array_range = scan_.get_array_range();
      for (int frame = array_range[0]; frame < array_range[1]; ++frame) {
        append_for_image(predictions, frame);
      }

      // Return the reflection table
      return table;
    }

    /**
     * Do the predictions for the given miller indices on the given frames.
     * @param h The miller indices
     * @param frame The frames
     * @returns The reflection table
     */
    af::reflection_table operator()(
        const af::const_ref<miller_index> &h,
        const af::const_ref<int> &frame) const {
      DIALS_ASSERT(h.size() == frame.size());

      // Create the table and local stuff
      af::reflection_table table;
      prediction_data predictions(table);

      // Sort the miller indices by frame
      af::shared<std::size_t> index(frame.size());
      for (std::size_t i = 0; i < index.size(); ++i) {
        index[i] = i;
      }
      std::sort(index.begin(), index.end(), sort_by_frame(frame));

      // For each hkl, do the prediction on the frame, because all the indices
      // are now sorted, we can calculate the matrix only once
      std::size_t i = 0;
      while (i < index.size()) {
        int current_frame = frame[index[i]];
        mat3<double> A1, A2;
        compute_setting_matrices(A1, A2, current_frame);
        while (i < index.size() && frame[index[i]] == current_frame) {
          append_for_index(predictions, A1, A2, current_frame, h[index[i]]);
          ++i;
        }
      }

      // Return the reflection table
      return table;
    }

    /**
     * Do the predictions for the given miller indices on the given frames.
     * @param h The miller indices
     * @param frame The frames
     * @param panel The panel
     * @returns The reflection table
     */
    af::reflection_table operator()(
        const af::const_ref<miller_index> &h,
        const af::const_ref<int> &frame,
        std::size_t panel) const {
      af::shared<std::size_t> panels(h.size(), panel);
      return (*this)(h, frame, panels.const_ref());
    }

    /**
     * Do the predictions for the given miller indices on the given frames.
     * @param h The miller indices
     * @param frame The frames
     * @param panel The panel
     * @returns The reflection table
     */
    af::reflection_table operator()(
        const af::const_ref<miller_index> &h,
        const af::const_ref<int> &frame,
        const af::const_ref<std::size_t> &panel) const {
      DIALS_ASSERT(h.size() == frame.size());
      DIALS_ASSERT(h.size() == panel.size());

      // Create the table and local stuff
      af::reflection_table table;
      prediction_data predictions(table);

      // Sort the miller indices by frame
      af::shared<std::size_t> index(frame.size());
      for (std::size_t i = 0; i < index.size(); ++i) {
        index[i] = i;
      }
      std::sort(index.begin(), index.end(), sort_by_frame(frame));

      // For each hkl, do the prediction on the frame, because all the indices
      // are now sorted, we can calculate the matrix only once
      std::size_t i = 0;
      while (i < index.size()) {
        int current_frame = frame[index[i]];
        mat3<double> A1, A2;
        compute_setting_matrices(A1, A2, current_frame);
        while (i < index.size() && frame[index[i]] == current_frame) {
          append_for_index(predictions, A1, A2, current_frame,
              h[index[i]], (int)panel[index[i]]);
          ++i;
        }
      }

      // Return the reflection table
      return table;
    }

  private:

    /**
     * Helper function to compute the setting matrix and the beginning and end
     * of a frame.
     */
    void compute_setting_matrices(
        mat3<double> &A1, mat3<double> &A2,
        int frame) const {

      // Get the rotation axis and beam vector
      vec3<double> m2 = goniometer_.get_rotation_axis();
      int frame_0 = scan_.get_array_range()[0];

      // Calculate the setting matrix at the beginning and end
      double phi_beg = scan_.get_angle_from_array_index(frame);
      double phi_end = scan_.get_angle_from_array_index(frame + 1);
      mat3<double> r_beg = axis_and_angle_as_matrix(m2, phi_beg);
      mat3<double> r_end = axis_and_angle_as_matrix(m2, phi_end);
      A1 = r_beg * A_[frame - frame_0];
      A2 = r_end * A_[frame - frame_0 + 1];
    }

    /**
     * For the given image, generate the indices and do the prediction.
     * @param p The reflection data
     * @param frame The image frame to predict on.
     */
    void append_for_image(prediction_data &p, int frame) const {

      // Get the rotation axis and beam vector
      vec3<double> m2 = goniometer_.get_rotation_axis();
      vec3<double> s0 = beam_.get_s0();
      mat3<double> A1, A2;
      compute_setting_matrices(A1, A2, frame);

      // Construct the index generator and do the predictions for each index
      ReekeIndexGenerator indices(A1, A2, m2, s0, dmin_, margin_);
      for (;;) {
        miller_index h = indices.next();
        if (h.is_zero()) {
          break;
        }
        append_for_index(p, A1, A2, frame, h);
      }
    }

    /**
     * Do the prediction for a miller index on a frame.
     * @param p The reflection data
     * @param A1 The beginning setting matrix.
     * @param A2 The end setting matrix
     * @param frame The frame to predict on
     * @param h The miller index
     */
    void append_for_index(prediction_data &p, mat3<double> A1, mat3<double> A2,
        std::size_t frame, const miller_index &h, int panel=-1) const {
      boost::optional<Ray> ray = predict_rays_(h, A1, A2, frame, 1);
      if (ray) {
        append_for_ray(p, h, *ray, panel);
      }
    }

    /**
     * Do the prediction for a given ray.
     * @param p The reflection data
     * @param h The miller index
     * @param ray The ray
     */
    void append_for_ray(prediction_data &p,
        const miller_index &h, const Ray &ray, int panel) const {
      try {

        // Get the impact on the detector
        Detector::coord_type impact = get_ray_intersection(ray.s1, panel);
        std::size_t panel = impact.first;
        vec2<double> mm = impact.second;
        vec2<double> px = detector_[panel].millimeter_to_pixel(mm);

        // Get the frame
        double frame = scan_.get_array_index_from_angle(ray.angle);

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
    af::shared< mat3<double> > A_;
    double dmin_;
    std::size_t margin_;
    ScanVaryingRayPredictor predict_rays_;
  };


  /**
   * A class to do naive stills prediction.
   */
  class StillsReflectionPredictor {

    typedef cctbx::miller::index<> miller_index;

  public:

    /**
     * Initialise the predictor
     */
    StillsReflectionPredictor(
        const Beam &beam,
        const Detector &detector,
        mat3<double> ub)
      : beam_(beam),
        detector_(detector),
        ub_(ub),
        predict_rays_(beam.get_s0()) {}

    /**
     * Predict all reflection.
     * @returns reflection table.
     */
    af::reflection_table operator()() const {
      DIALS_ERROR("Not implemented");
      return af::reflection_table();
    }

    /**
     * Predict the reflections with given HKL.
     * @param h The miller index
     * @returns The reflection list
     */
    af::reflection_table operator()(
        const af::const_ref< miller_index > &h) const {
      af::reflection_table table;
      prediction_data predictions(table);
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predictions, h[i]);
      }
      return table;
    }

    /**
     * Predict for given hkl and panel.
     * @param h The miller index
     * @param panel The panel
     * @returns The reflection table
     */
    af::reflection_table operator()(
        const af::const_ref< miller_index > &h,
        std::size_t panel) const {
      af::shared<std::size_t> panels(h.size(), panel);
      return (*this)(h, panels.const_ref());
    }

    /**
     * Predict for given hkl and panel.
     * @param h The miller index
     * @param panel The panel
     * @returns The reflection table
     */
    af::reflection_table operator()(
        const af::const_ref< miller_index > &h,
        const af::const_ref<std::size_t> &panel) const {
      DIALS_ASSERT(h.size() == panel.size());
      af::reflection_table table;
      prediction_data predictions(table);
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predictions, h[i], (int)panel[i]);
      }
      return table;
    }

  private:

    /**
     * Predict reflections for the given HKL.
     * @param p The reflection data
     * @param h The miller index
     */
    void append_for_index(prediction_data &p,
        const miller_index &h, int panel=-1) const {
      af::small<Ray, 2> rays = predict_rays_(h, ub_);
      for (std::size_t i = 0; i < rays.size(); ++i) {
        append_for_ray(p, h, rays[i], panel);
      }
    }

    /**
     * Predict the reflection for the given ray data.
     * @param p The reflection data
     * @param h The miller index
     * @param ray The ray data
     */
    void append_for_ray(prediction_data &p,
        const miller_index &h, const Ray &ray, int panel) const {
      try {

        // Get the impact on the detector
        Detector::coord_type impact = get_ray_intersection(ray.s1, panel);
        std::size_t panel = impact.first;
        vec2<double> mm = impact.second;
        vec2<double> px = detector_[panel].millimeter_to_pixel(mm);

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
    mat3<double> ub_;
    StillsRayPredictor predict_rays_;
  };
}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_REFLECTION_PREDICTOR_H
