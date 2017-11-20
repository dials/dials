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
#include <scitbx/constants.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/scan_helpers.h>
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
  using scitbx::constants::pi;
  using scitbx::constants::two_pi;
  using scitbx::constants::pi_180;
  using dxtbx::model::BeamBase;
  using dxtbx::model::Detector;
  using dxtbx::model::Panel;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;
  using dxtbx::model::is_angle_in_range;
  using dials::model::Ray;
  using dxtbx::model::plane_ray_intersection;

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
    af::shared< std::size_t > flags;

    prediction_data(af::reflection_table &table) {
      hkl    = table.get< miller_index >("miller_index");
      panel  = table.get< std::size_t >("panel");
      enter  = table.get< bool >("entering");
      s1     = table.get< vec3<double> >("s1");
      xyz_px = table.get< vec3<double> >("xyzcal.px");
      xyz_mm = table.get< vec3<double> >("xyzcal.mm");
      flags  = table.get< std::size_t >("flags");
    }
  };

  struct stills_prediction_data : prediction_data {

    af::shared< double > delpsi;

    stills_prediction_data(af::reflection_table &table)
      : prediction_data(table) {
      delpsi = table.get< double >("delpsical.rad");
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
          const boost::shared_ptr<BeamBase> beam,
          const Detector &detector,
          const Goniometer &goniometer,
          const Scan &scan,
          const cctbx::uctbx::unit_cell &unit_cell,
          const cctbx::sgtbx::space_group_type &space_group_type,
          double dmin,
          double padding)
      : beam_(beam),
        detector_(detector),
        goniometer_(goniometer),
        scan_(scan),
        unit_cell_(unit_cell),
        space_group_type_(space_group_type),
        dmin_(dmin),
        padding_(padding),
        predict_rays_(
            beam->get_s0(),
            goniometer.get_rotation_axis_datum(),
            goniometer.get_fixed_rotation(),
            goniometer.get_setting_rotation(),
            vec2<double>(0.0, two_pi)){
      DIALS_ASSERT(padding >= 0);
    }


    /**
     * Predict reflections for UB
     * @param ub The UB matrix
     * @returns A reflection table.
     */
    af::reflection_table for_ub(const mat3<double> &ub) const {

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
     * @param h The array of miller indices
     * @param entering The array of entering flags
     * @param panel The array of panels
     * @param ub The array of UB matrices
     * @returns A reflection table.
     */
    af::reflection_table for_hkl(
        const af::const_ref<miller_index> &h,
        const af::const_ref<bool> &entering,
        const af::const_ref<std::size_t> &panel,
        const mat3<double> &ub) const {
      af::shared< mat3<double> > uba(h.size(), ub);
      return for_hkl_with_individual_ub(h, entering, panel, uba.const_ref());
    }

    /**
     * Predict reflections for UB
     * @param h The array of miller indices
     * @param entering The array of entering flags
     * @param panel The array of panels
     * @param ub The array of UB matrices
     * @returns A reflection table.
     */
    af::reflection_table for_hkl_with_individual_ub(
        const af::const_ref<miller_index> &h,
        const af::const_ref<bool> &entering,
        const af::const_ref<std::size_t> &panel,
        const af::const_ref< mat3<double> >&ub) const {
      DIALS_ASSERT(ub.size() == h.size());
      DIALS_ASSERT(ub.size() == panel.size());
      DIALS_ASSERT(ub.size() == entering.size());
      DIALS_ASSERT(scan_.get_oscillation()[1] > 0.0);
      af::reflection_table table;
      prediction_data predictions(table);
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predictions, ub[i], h[i], entering[i], panel[i]);
      }
      DIALS_ASSERT(table.nrows() == h.size());
      return table;
    }

    /**
     * Predict reflections to the entries in the table.
     * @param table The reflection table
     * @param ub The ub matrix
     */
    void for_reflection_table(
        af::reflection_table table,
        const mat3<double> &ub) const {
      af::shared< mat3<double> > uba(table.nrows(), ub);
      for_reflection_table_with_individual_ub(table, uba.const_ref());
    }

    /**
     * Predict reflections to the entries in the table.
     * @param table The reflection table
     * @param ub The list of ub matrices
     */
    void for_reflection_table_with_individual_ub(
        af::reflection_table table,
        const af::const_ref< mat3<double> > &ub) const {
      DIALS_ASSERT(ub.size() == table.nrows());
      af::reflection_table new_table = for_hkl_with_individual_ub(
        table["miller_index"],
        table["entering"],
        table["panel"],
        ub);
      DIALS_ASSERT(new_table.nrows() == table.nrows());
      table["miller_index"] = new_table["miller_index"];
      table["entering"] = new_table["entering"];
      table["panel"] = new_table["panel"];
      table["s1"] = new_table["s1"];
      table["xyzcal.px"] = new_table["xyzcal.px"];
      table["xyzcal.mm"] = new_table["xyzcal.mm"];
      af::shared<std::size_t> flags = table["flags"];
      af::shared<std::size_t> new_flags = new_table["flags"];
      for (std::size_t i = 0; i < flags.size(); ++i) {
        flags[i] &= ~af::Predicted;
        flags[i] |= new_flags[i];
      }
      DIALS_ASSERT(table.is_consistent());
    }

  private:

    void append_for_index(
        prediction_data &p,
        const mat3<double> ub,
        const miller_index &h) const {
      af::small<Ray, 2> rays = predict_rays_(h, ub);
      for (std::size_t i = 0; i < rays.size(); ++i) {
        try {
          Detector::coord_type impact = detector_.get_ray_intersection(rays[i].s1);
          std::size_t panel = impact.first;
          vec2<double> mm = impact.second;
          vec2<double> px = detector_[panel].millimeter_to_pixel(mm);
          af::shared< vec2<double> > frames =
            scan_.get_array_indices_with_angle(rays[i].angle, padding_, true);
          for (std::size_t j = 0; j < frames.size(); ++j) {
            p.hkl.push_back(h);
            p.enter.push_back(rays[i].entering);
            p.s1.push_back(rays[i].s1);
            p.panel.push_back(panel);
            p.flags.push_back(af::Predicted);
            p.xyz_mm.push_back(vec3<double>(mm[0], mm[1], frames[j][0]));
            p.xyz_px.push_back(vec3<double>(px[0], px[1], frames[j][1]));
          }
        } catch(dxtbx::error) {
          // do nothing
        }
      }
    }

    void append_for_index(
        prediction_data &p,
        const mat3<double> ub,
        const miller_index &h,
        bool entering,
        std::size_t panel) const {
      p.hkl.push_back(h);
      p.enter.push_back(entering);
      p.panel.push_back(panel);
      af::small<Ray, 2> rays = predict_rays_(h, ub);
      for (std::size_t i = 0; i < rays.size(); ++i) {
        if (rays[i].entering == entering) {
          p.s1.push_back(rays[i].s1);
          double frame = scan_.get_array_index_from_angle(rays[i].angle);
          try {
            vec2<double> mm = detector_[panel].get_ray_intersection(rays[i].s1);
            vec2<double> px = detector_[panel].millimeter_to_pixel(mm);
            p.xyz_mm.push_back(vec3<double>(mm[0], mm[1], rays[i].angle));
            p.xyz_px.push_back(vec3<double>(px[0], px[1], frame));
            p.flags.push_back(af::Predicted);
          } catch(dxtbx::error) {
            p.xyz_mm.push_back(vec3<double>(0, 0, rays[i].angle));
            p.xyz_px.push_back(vec3<double>(0, 0, frame));
            p.flags.push_back(0);
          }
          return;
        }
      }
      p.s1.push_back(vec3<double>(0, 0, 0));
      p.xyz_mm.push_back(vec3<double>(0, 0, 0));
      p.xyz_px.push_back(vec3<double>(0, 0, 0));
      p.flags.push_back(0);
    }

    boost::shared_ptr<BeamBase> beam_;
    Detector detector_;
    Goniometer goniometer_;
    Scan scan_;
    cctbx::uctbx::unit_cell unit_cell_;
    cctbx::sgtbx::space_group_type space_group_type_;
    double dmin_;
    double padding_;
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
        const boost::shared_ptr<BeamBase> beam,
        const Detector &detector,
        const Goniometer &goniometer,
        const Scan &scan,
        const cctbx::sgtbx::space_group_type &space_group_type,
        double dmin,
        std::size_t margin,
        double padding)
      : beam_(beam),
        detector_(detector),
        goniometer_(goniometer),
        scan_(scan),
        space_group_type_(space_group_type),
        dmin_(dmin),
        margin_(margin),
        padding_(padding),
        predict_rays_(
            beam->get_s0(),
            goniometer.get_rotation_axis_datum(),
            scan.get_array_range()[0],
            scan.get_oscillation(),
            dmin) {}

    /**
     * Return the beam model
     */
    boost::shared_ptr<BeamBase> beam() const {
      return beam_;
    }

    /**
     * Return the detector model
     */
    Detector detector() const {
      return detector_;
    }

    /**
     * Return the goniometer model
     */
    Goniometer goniometer() const {
      return goniometer_;
    }

    /**
     * Return the scan model
     */
    Scan scan() const {
      return scan_;
    }

    /**
     * Predict all the reflections for this model.
     * @returns The reflection table
     */
    af::reflection_table for_ub(const af::const_ref< mat3<double> > &A) const {
      DIALS_ASSERT(A.size() == scan_.get_num_images() + 1);

      // Create the table and local stuff
      af::reflection_table table;
      prediction_data predictions(table);

      // Get the array range and loop through all the images
      vec2<int> array_range = scan_.get_array_range();
      double a0 = scan_.get_oscillation_range()[0];
      double a1 = scan_.get_oscillation_range()[1];
      int z0 = std::floor(scan_.get_array_index_from_angle(a0-padding_*pi/180.0) + 0.5);
      int z1 = std::floor(scan_.get_array_index_from_angle(a1+padding_*pi/180.0) + 0.5);
      const int offset = array_range[0];
      for (int frame = z0; frame < z1; ++frame) {
        int i = frame - offset;
        if (i < 0) i = 0;
        if (i >= A.size()-1) i = A.size() - 2;
        append_for_image(predictions, frame, A[i], A[i+1]);
      }

      // Return the reflection table
      return table;
    }

    /**
     * Predict all the reflections for this model.
     * @returns The reflection table
     */
    af::reflection_table for_ub_on_single_image(
        int frame,
        const mat3<double> &A1,
        const mat3<double> &A2) const {
      vec2<int> array_range = scan_.get_array_range();
      DIALS_ASSERT(frame >= array_range[0] && frame < array_range[1]);

      // Create the table and local stuff
      af::reflection_table table;
      prediction_data predictions(table);

      // Get the array range and loop through all the images
      append_for_image(predictions, frame, A1, A2);

      // Return the reflection table
      return table;
    }

    /**
     * Predict reflections for UB
     * @param h The array of miller indices
     * @param entering The array of entering flags
     * @param panel The array of panels
     * @param ub The array of UB matrices
     * @param s0 The array of s0 vectors
     * @param d The array of d matrices
     * @returns A reflection table.
     */
    af::reflection_table for_hkl_with_individual_ub(
        const af::const_ref<miller_index> &h,
        const af::const_ref<bool> &entering,
        const af::const_ref<std::size_t> &panel,
        const af::const_ref< mat3<double> >&ub,
        const af::const_ref< vec3<double> > &s0,
        const af::const_ref< mat3<double> > &d) const {
      DIALS_ASSERT(ub.size() == h.size());
      DIALS_ASSERT(ub.size() == panel.size());
      DIALS_ASSERT(ub.size() == entering.size());
      DIALS_ASSERT(ub.size() == s0.size());
      DIALS_ASSERT(ub.size() == d.size());
      DIALS_ASSERT(scan_.get_oscillation()[1] > 0.0);
      af::reflection_table table;
      prediction_data predictions(table);
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predictions, ub[i], s0[i], d[i], h[i], entering[i], panel[i]);
      }
      DIALS_ASSERT(table.nrows() == h.size());
      return table;
    }

    /**
     * Predict reflections to the entries in the table.
     * @param table The reflection table
     * @param ub The ub matrix array
     * @param s0 The s0 vector array
     * @param d The d matrix array
     */
    void for_reflection_table(
        af::reflection_table table,
        const af::const_ref< mat3<double> > &ub,
        const af::const_ref< vec3<double> > &s0,
        const af::const_ref< mat3<double> > &d
        ) const {
      DIALS_ASSERT(ub.size() == table.nrows());
      DIALS_ASSERT(s0.size() == table.nrows());
      DIALS_ASSERT(d.size() == table.nrows());
      af::reflection_table new_table = for_hkl_with_individual_ub(
        table["miller_index"],
        table["entering"],
        table["panel"],
        ub,
        s0,
        d);
      DIALS_ASSERT(new_table.nrows() == table.nrows());
      table["miller_index"] = new_table["miller_index"];
      table["entering"] = new_table["entering"];
      table["panel"] = new_table["panel"];
      table["s1"] = new_table["s1"];
      table["xyzcal.px"] = new_table["xyzcal.px"];
      table["xyzcal.mm"] = new_table["xyzcal.mm"];
      af::shared<std::size_t> flags = table["flags"];
      af::shared<std::size_t> new_flags = new_table["flags"];
      for (std::size_t i = 0; i < flags.size(); ++i) {
        flags[i] &= ~af::Predicted;
        flags[i] |= new_flags[i];
      }
      DIALS_ASSERT(table.is_consistent());

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
      vec3<double> m2 = goniometer_.get_rotation_axis_datum();

      // Calculate the setting matrix at the beginning and end
      double phi_beg = scan_.get_angle_from_array_index(frame);
      double phi_end = scan_.get_angle_from_array_index(frame + 1);
      mat3<double> r_fixed =  goniometer_.get_fixed_rotation();
      mat3<double> r_setting =  goniometer_.get_setting_rotation();
      mat3<double> r_beg = axis_and_angle_as_matrix(m2, phi_beg);
      mat3<double> r_end = axis_and_angle_as_matrix(m2, phi_end);
      A1 = r_setting * r_beg * r_fixed * A1;
      A2 = r_setting * r_end * r_fixed * A2;
    }

    /**
     * For the given image, generate the indices and do the prediction.
     * @param p The reflection data
     * @param frame The image frame to predict on.
     */
    void append_for_image(
        prediction_data &p, int frame,
        mat3<double> A1,
        mat3<double> A2) const {

      // Get the rotation axis and beam vector
      vec3<double> m2 = goniometer_.get_rotation_axis_datum();
      vec3<double> s0 = beam_->get_s0();
      compute_setting_matrices(A1, A2, frame);

      // Construct the index generator and do the predictions for each index
      ReekeIndexGenerator indices(A1, A2, space_group_type_, m2, s0, dmin_, margin_);
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
        Detector::coord_type impact = detector_.get_ray_intersection(ray.s1);
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
        p.flags.push_back(af::Predicted);

      } catch(dxtbx::error) {
        // do nothing
      }
    }

    /**
     * Do the prediction for a miller index given all model states.
     * @param p The reflection data
     * @param ub The ub matrix array
     * @param s0 The s0 vector array
     * @param d The d matrix array
     * @param h The miller index array
     * @param entering The entering flag array
     * @param panel The array of panel numbers
     */
    void append_for_index(
        prediction_data &p,
        const mat3<double> ub,
        const vec3<double> s0,
        const mat3<double> d,
        const miller_index &h,
        bool entering,
        std::size_t panel) const {
      p.hkl.push_back(h);
      p.enter.push_back(entering);
      p.panel.push_back(panel);
      // Need a local ray predictor for just this reflection's s0
      ScanStaticRayPredictor local_predict_rays_(
        s0,
        goniometer_.get_rotation_axis_datum(),
        goniometer_.get_fixed_rotation(),
        goniometer_.get_setting_rotation(),
        vec2<double>(0.0, two_pi));
      af::small<Ray, 2> rays = local_predict_rays_(h, ub);
      for (std::size_t i = 0; i < rays.size(); ++i) {
        if (rays[i].entering == entering) {
          p.s1.push_back(rays[i].s1);
          double frame = scan_.get_array_index_from_angle(rays[i].angle);
          try {
            // Need a local panel with the right D matrix
            Panel local_panel(detector_[panel]);
            local_panel.set_frame(d.get_column(0), d.get_column(1), d.get_column(2));
            vec2<double> mm = local_panel.get_ray_intersection(rays[i].s1);
            vec2<double> px = local_panel.millimeter_to_pixel(mm);
            p.xyz_mm.push_back(vec3<double>(mm[0], mm[1], rays[i].angle));
            p.xyz_px.push_back(vec3<double>(px[0], px[1], frame));
            p.flags.push_back(af::Predicted);
          } catch(dxtbx::error) {
            p.xyz_mm.push_back(vec3<double>(0, 0, rays[i].angle));
            p.xyz_px.push_back(vec3<double>(0, 0, frame));
            p.flags.push_back(0);
          }
          return;
        }
      }
      p.s1.push_back(vec3<double>(0, 0, 0));
      p.xyz_mm.push_back(vec3<double>(0, 0, 0));
      p.xyz_px.push_back(vec3<double>(0, 0, 0));
      p.flags.push_back(0);
    }

    boost::shared_ptr<BeamBase> beam_;
    Detector detector_;
    Goniometer goniometer_;
    Scan scan_;
    cctbx::sgtbx::space_group_type space_group_type_;
    double dmin_;
    std::size_t margin_;
    double padding_;
    ScanVaryingRayPredictor predict_rays_;
  };


  /**
   * A class to do reflection prediction for either static or scan varying
   */
  class ScansReflectionPredictor {
  public:

    /**
     * Initialise the predictor
     */
    ScansReflectionPredictor(
            const boost::shared_ptr<BeamBase> beam,
            const Detector &detector,
            const Goniometer &goniometer,
            const Scan &scan,
            const cctbx::sgtbx::space_group_type &space_group_type,
            double dmin,
            std::size_t margin,
            double padding)
      : predict_(
          beam,
          detector,
          goniometer,
          scan,
          space_group_type,
          dmin,
          margin,
          padding) {}



    af::reflection_table for_static_ub(
        const mat3<double> &A) const {
      af::shared< mat3<double> > A_list(predict_.scan().get_num_images() + 1, A);
      return predict_.for_ub(A_list.const_ref());
    }

    af::reflection_table for_varying_ub(
        const af::const_ref< mat3<double> > &A) const {
      return predict_.for_ub(A);
    }

    af::reflection_table for_static_ub_on_single_image(
        int frame,
        const mat3<double> &A) const {
      return predict_.for_ub_on_single_image(frame, A, A);
    }

    af::reflection_table for_varying_ub_on_single_image(
        int frame,
        const mat3<double> &A1,
        const mat3<double> &A2) const {
      return predict_.for_ub_on_single_image(frame, A1, A2);
    }

    af::reflection_table for_hkl_with_static_ub(
        const af::const_ref<miller_index> &h,
        const af::const_ref<bool> &entering,
        const af::const_ref<std::size_t> &panel,
        const mat3<double> &A) const {
      af::shared< mat3<double> > A_list(h.size(), A);
      af::shared< vec3<double> > s0_list(h.size(), predict_.beam()->get_s0());
      af::shared< mat3<double> > d_list(h.size());
      for (std::size_t i = 0; i < d_list.size(); ++i) {
        d_list[i] = predict_.detector()[panel[i]].get_d_matrix();
      }
      return predict_.for_hkl_with_individual_ub(
          h,
          entering,
          panel,
          A_list.const_ref(),
          s0_list.const_ref(),
          d_list.const_ref());
    }

    af::reflection_table for_hkl_with_individual_ub(
        const af::const_ref<miller_index> &h,
        const af::const_ref<bool> &entering,
        const af::const_ref<std::size_t> &panel,
        const af::const_ref< mat3<double> > &A) const {
      af::shared< vec3<double> > s0_list(h.size(), predict_.beam()->get_s0());
      af::shared< mat3<double> > d_list(h.size());
      for (std::size_t i = 0; i < d_list.size(); ++i) {
        d_list[i] = predict_.detector()[panel[i]].get_d_matrix();
      }
      return predict_.for_hkl_with_individual_ub(
          h,
          entering,
          panel,
          A,
          s0_list.const_ref(),
          d_list.const_ref());
    }

    af::reflection_table for_hkl_with_individual_ub_beam_and_d_matrix(
        const af::const_ref<miller_index> &h,
        const af::const_ref<bool> &entering,
        const af::const_ref<std::size_t> &panel,
        const af::const_ref< mat3<double> >&A,
        const af::const_ref< vec3<double> > &s0,
        const af::const_ref< mat3<double> > &d) const {
      return predict_.for_hkl_with_individual_ub(h, entering, panel, A, s0, d);
    }

    void for_reflection_table_with_static_ub(
        af::reflection_table table,
        const mat3<double> &A) const {
      af::shared< mat3<double> > A_list(table.size(), A);
      af::shared< vec3<double> > s0_list(table.size(), predict_.beam()->get_s0());
      af::shared< mat3<double> > d_list(table.size());
      af::const_ref<std::size_t> panel = table["panel"];
      for (std::size_t i = 0; i < d_list.size(); ++i) {
        d_list[i] = predict_.detector()[panel[i]].get_d_matrix();
      }
      predict_.for_reflection_table(
          table,
          A_list.const_ref(),
          s0_list.const_ref(),
          d_list.const_ref());
    }

    void for_reflection_table_with_individual_ub(
        af::reflection_table table,
        const af::const_ref< mat3<double> > &A) const {
      af::shared< vec3<double> > s0_list(table.size(), predict_.beam()->get_s0());
      af::shared< mat3<double> > d_list(table.size());
      af::const_ref<std::size_t> panel = table["panel"];
      for (std::size_t i = 0; i < d_list.size(); ++i) {
        d_list[i] = predict_.detector()[panel[i]].get_d_matrix();
      }
      predict_.for_reflection_table(
          table,
          A,
          s0_list.const_ref(),
          d_list.const_ref());
    }

    void for_reflection_table_with_individual_ub_beam_and_d_matrix(
        af::reflection_table table,
        const af::const_ref< mat3<double> > &A,
        const af::const_ref< vec3<double> > &s0,
        const af::const_ref< mat3<double> > &d) const {
      predict_.for_reflection_table(table, A, s0, d);
    }

  protected:

    ScanVaryingReflectionPredictor predict_;

  };


  /**
   * A class to do stills prediction.
   */
  class StillsDeltaPsiReflectionPredictor {
  public:
    typedef cctbx::miller::index<> miller_index;

    /**
     * Initialise the predictor
     */
    StillsDeltaPsiReflectionPredictor(
        const boost::shared_ptr<BeamBase> beam,
        const Detector &detector,
        mat3<double> ub,
        const cctbx::uctbx::unit_cell &unit_cell,
        const cctbx::sgtbx::space_group_type &space_group_type,
        const double &dmin)
      : beam_(beam),
        detector_(detector),
        ub_(ub),
        unit_cell_(unit_cell),
        space_group_type_(space_group_type),
        dmin_(dmin),
        predict_ray_(beam->get_s0()) {}

    /**
     * Predict all reflection.
     * @returns reflection table.
     */
    af::reflection_table operator()() const {
      DIALS_ERROR("Not implemented");
      return af::reflection_table();
    }

    /**
     * Predict reflections for UB. Also filters based on ewald sphere proximity.
     * @param ub The UB matrix
     * @returns A reflection table.
     */
    af::reflection_table for_ub(const mat3<double> &ub) {

      // Create the reflection table and the local container
      af::reflection_table table;
      stills_prediction_data predictions(table);

      // Create the index generate and loop through the indices. For each index,
      // predict the rays and append to the reflection table
      IndexGenerator indices(unit_cell_, space_group_type_, dmin_);
      for (;;) {
        miller_index h = indices.next();
        if (h.is_zero()) {
          break;
        }

        Ray ray;
        ray = predict_ray_(h, ub);
        double delpsi = std::abs(predict_ray_.get_delpsi());
        if(delpsi < 0.0015)
          append_for_index(predictions, ub, h);
      }

      // Return the reflection table
      return table;
    }

    /**
     * Predict the reflections with given HKL.
     * @param h The miller index
     * @returns The reflection list
     */
    af::reflection_table operator()(
        const af::const_ref< miller_index > &h) {
      af::reflection_table table;
      stills_prediction_data predictions(table);
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predictions, ub_, h[i]);
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
        std::size_t panel) {
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
        const af::const_ref<std::size_t> &panel) {
      DIALS_ASSERT(h.size() == panel.size());
      af::reflection_table table;
      stills_prediction_data predictions(table);
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predictions, ub_, h[i], (int)panel[i]);
      }
      return table;
    }

    /**
     * Predict reflections for UB
     * @param h The array of miller indices
     * @param panel The array of panels
     * @returns A reflection table.
     */
    af::reflection_table for_hkl_with_individual_ub(
        const af::const_ref<miller_index> &h,
        const af::const_ref<std::size_t> &panel,
        const af::const_ref< mat3<double> >&ub) {
      DIALS_ASSERT(ub.size() == h.size());
      DIALS_ASSERT(ub.size() == panel.size());
      af::reflection_table table;
      af::shared<double> column;
      table["delpsical.rad"] = column;
      stills_prediction_data predictions(table);
      for (std::size_t i = 0; i < h.size(); ++i) {
        append_for_index(predictions, ub[i], h[i], panel[i]);
      }
      DIALS_ASSERT(table.nrows() == h.size());
      return table;
    }

    /**
     * Predict reflections to the entries in the table.
     * @param table The reflection table
     * @param ub The ub matrix
     */
    void for_reflection_table(
        af::reflection_table table,
        const mat3<double> &ub) {
      af::shared< mat3<double> > uba(table.nrows(), ub);
      for_reflection_table_with_individual_ub(table, uba.const_ref());
    }

    /**
     * Predict reflections and add to the entries in the table.
     * @param table The reflection table
     */
    void for_reflection_table_with_individual_ub(
        af::reflection_table table,
        const af::const_ref< mat3<double> > &ub) {
      DIALS_ASSERT(ub.size() == table.nrows());
      af::reflection_table new_table = for_hkl_with_individual_ub(
        table["miller_index"],
        table["panel"],
        ub);
      DIALS_ASSERT(new_table.nrows() == table.nrows());
      table["miller_index"] = new_table["miller_index"];
      table["panel"] = new_table["panel"];
      table["s1"] = new_table["s1"];
      table["xyzcal.px"] = new_table["xyzcal.px"];
      table["xyzcal.mm"] = new_table["xyzcal.mm"];

      // Add "delpsical.rad" key to table if it is not there already
      //if (table.count("delpsical.rad") != 1) {
      //  af::shared<double> column(table.nrows());
      //  table["delpsical.rad"] = column;
      //}
      table["delpsical.rad"] = new_table["delpsical.rad"];
      af::shared<std::size_t> flags = table["flags"];
      af::shared<std::size_t> new_flags = new_table["flags"];
      for (std::size_t i = 0; i < flags.size(); ++i) {
        flags[i] &= ~af::Predicted;
        flags[i] |= new_flags[i];
      }
      DIALS_ASSERT(table.is_consistent());
    }

  protected:

    /**
     * Predict reflections for the given HKL.
     * @param p The reflection data
     * @param h The miller index
     */
    virtual void append_for_index(stills_prediction_data &p,
        const mat3<double> ub,
        const miller_index &h, int panel=-1) {
      //af::small<Ray, 2> rays = predict_rays_(h, ub_);
      Ray ray;
      ray = predict_ray_(h, ub);
      double delpsi = predict_ray_.get_delpsi();
      append_for_ray(p, h, ray, panel, delpsi);
    }

    /**
     * Predict the reflection for the given ray data.
     * @param p The reflection data
     * @param h The miller index
     * @param ray The ray data
     * @param panel The panel number
     * @param delpsi The calculated minimum rotation to Ewald sphere (delPsi) value
     */
    void append_for_ray(stills_prediction_data &p,
        const miller_index &h, const Ray &ray, int panel, double delpsi) const {
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
        p.flags.push_back(af::Predicted);
        p.delpsi.push_back(delpsi);

      } catch(dxtbx::error) {
        // do nothing
      }
    }

  private:

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

  protected:

    boost::shared_ptr<BeamBase> beam_;
    Detector detector_;
    mat3<double> ub_;
    cctbx::uctbx::unit_cell unit_cell_;
    cctbx::sgtbx::space_group_type space_group_type_;
    const double dmin_;
    StillsRayPredictor predict_ray_;
  };

  class NaveStillsReflectionPredictor : public StillsDeltaPsiReflectionPredictor {

    typedef cctbx::miller::index<> miller_index;
    /**
      * Initialise the predictor
      */
  public:

    NaveStillsReflectionPredictor(
        const boost::shared_ptr<BeamBase> beam,
        const Detector &detector,
        mat3<double> ub,
        const cctbx::uctbx::unit_cell &unit_cell,
        const cctbx::sgtbx::space_group_type &space_group_type,
        const double &dmin,
        const double &ML_half_mosaicity_deg,
        const double &ML_domain_size_ang)
      : StillsDeltaPsiReflectionPredictor(beam,
        detector, ub, unit_cell, space_group_type, dmin),
        ML_half_mosaicity_deg_(ML_half_mosaicity_deg),
        ML_domain_size_ang_(ML_domain_size_ang) {}

    /**
     * Predict reflections for UB, using the Nave models from Nave 2014, JSR
     * as implemented in Sauter 2014, Akta D, equations 16-17
     * @param ub The UB matrix
     * @returns A reflection table.
     */
    af::reflection_table for_ub(const mat3<double> &ub) {

      // Create the reflection table and the local container
      af::reflection_table table;
      stills_prediction_data predictions(table);

      // Create the index generate and loop through the indices. For each index,
      // predict the rays and append to the reflection table
      IndexGenerator indices(unit_cell_, space_group_type_, dmin_);
      for (;;) {
        miller_index h = indices.next();
        if (h.is_zero()) {
          break;
        }
        double d = unit_cell_.d(h);
        double deltapsi_model = (d/ML_domain_size_ang_) + (ML_half_mosaicity_deg_*pi_180/2); // equation 16

        Ray ray;
        ray = predict_ray_(h, ub);
        double delpsi = std::abs(predict_ray_.get_delpsi());
        if(delpsi < deltapsi_model) // equation 17
          append_for_index(predictions, ub, h);
      }

      // Return the reflection table
      return table;
    }

  private:
    const double ML_half_mosaicity_deg_;
    const double ML_domain_size_ang_;
  };

  /**
   * Stills prediction using the SphericalRelpStillsRayPredictor.
   */
  class SphericalRelpStillsReflectionPredictor : public StillsDeltaPsiReflectionPredictor {

  typedef cctbx::miller::index<> miller_index;

  public:

    /**
      * Initialise the predictor
      */
    SphericalRelpStillsReflectionPredictor(
        const boost::shared_ptr<BeamBase> beam,
        const Detector &detector,
        mat3<double> ub,
        const cctbx::uctbx::unit_cell &unit_cell,
        const cctbx::sgtbx::space_group_type &space_group_type,
        const double &dmin)
      : StillsDeltaPsiReflectionPredictor(beam,
        detector, ub, unit_cell, space_group_type, dmin),
        spherical_relp_predict_ray_(beam->get_s0()) {}

    /**
     * Predict reflections for UB. Also filters based on ewald sphere proximity.
     * Override uses SphericalRelpStillsRayPredictor.
     * @param ub The UB matrix
     * @returns A reflection table.
     */
    af::reflection_table for_ub(const mat3<double> &ub) {

      // Create the reflection table and the local container
      af::reflection_table table;
      stills_prediction_data predictions(table);

      // Create the index generate and loop through the indices. For each index,
      // predict the rays and append to the reflection table
      IndexGenerator indices(unit_cell_, space_group_type_, dmin_);
      for (;;) {
        miller_index h = indices.next();
        if (h.is_zero()) {
          break;
        }

        Ray ray;
        ray = spherical_relp_predict_ray_(h, ub);
        double delpsi = std::abs(spherical_relp_predict_ray_.get_delpsi());
        if(delpsi < 0.0015)
          append_for_index(predictions, ub, h);
      }

      // Return the reflection table
      return table;
    }

  protected:

    /**
     * Predict reflections for the given HKL.
     * Override uses SphericalRelpStillsRayPredictor.
     * @param p The reflection data
     * @param h The miller index
     */
    virtual void append_for_index(stills_prediction_data &p,
        const mat3<double> ub,
        const miller_index &h, int panel=-1) {
      //af::small<Ray, 2> rays = predict_rays_(h, ub_);
      Ray ray;
      ray = spherical_relp_predict_ray_(h, ub);
      double delpsi = spherical_relp_predict_ray_.get_delpsi();
      append_for_ray(p, h, ray, panel, delpsi);
    }

    SphericalRelpStillsRayPredictor spherical_relp_predict_ray_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_REFLECTION_PREDICTOR_H
