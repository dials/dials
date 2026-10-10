/*
 * bbox_calculator.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_BBOX_CALCULATOR_H
#define DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_BBOX_CALCULATOR_H

#include <memory>
#include <cmath>
#include <scitbx/constants.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/profile_model/gaussian_rs/coordinate_system.h>

namespace dials { namespace algorithms { namespace profile_model {
  namespace gaussian_rs {

    // Use a load of stuff from other namespaces
    using dxtbx::model::BeamBase;
    using dxtbx::model::Detector;
    using dxtbx::model::Goniometer;
    using dxtbx::model::Scan;
    using scitbx::vec2;
    using scitbx::vec3;
    using scitbx::af::double2;
    using scitbx::af::double3;
    using scitbx::af::double4;
    using scitbx::af::int6;
    using scitbx::af::max;
    using scitbx::af::min;
    using std::ceil;
    using std::floor;

    /**
     * Interface for bounding box calculator.
     */
    class BBoxCalculatorIface {
    public:
      virtual ~BBoxCalculatorIface() {}

      virtual int6 single(vec3<double> s1, double frame, std::size_t panel) const = 0;

      virtual af::shared<int6> array(const af::const_ref<vec3<double> >& s1,
                                     const af::const_ref<double>& frame,
                                     const af::const_ref<std::size_t>& panel) const = 0;
    };

    /** Calculate the bounding box for each reflection */
    class BBoxCalculator3D : public BBoxCalculatorIface {
    public:
      /**
       * Initialise the bounding box calculation.
       * @param beam The beam parameters
       * @param detector The detector parameters
       * @param goniometer The goniometer parameters
       * @param delta_divergence The xds delta_divergence parameter
       * @param delta_mosaicity The xds delta_mosaicity parameter
       */
      BBoxCalculator3D(const BeamBase& beam,
                       const Detector& detector,
                       const Goniometer& gonio,
                       const Scan& scan,
                       double delta_divergence,
                       double delta_mosaicity)
          : s0_(beam.get_s0()),
            m2_(gonio.get_rotation_axis()),
            detector_(detector),
            scan_(scan),
            delta_divergence_(1, delta_divergence),
            delta_mosaicity_(1, delta_mosaicity) {
        DIALS_ASSERT(delta_divergence > 0.0);
        DIALS_ASSERT(delta_mosaicity > 0.0);
      }

      /**
       * Initialise the bounding box calculation.
       * @param beam The beam parameters
       * @param detector The detector parameters
       * @param goniometer The goniometer parameters
       * @param delta_divergence The xds delta_divergence parameter
       * @param delta_mosaicity The xds delta_mosaicity parameter
       */
      BBoxCalculator3D(const BeamBase& beam,
                       const Detector& detector,
                       const Goniometer& gonio,
                       const Scan& scan,
                       const af::const_ref<double>& delta_divergence,
                       const af::const_ref<double>& delta_mosaicity)
          : s0_(beam.get_s0()),
            m2_(gonio.get_rotation_axis()),
            detector_(detector),
            scan_(scan),
            delta_divergence_(delta_divergence.begin(), delta_divergence.end()),
            delta_mosaicity_(delta_mosaicity.begin(), delta_mosaicity.end()) {
        DIALS_ASSERT(delta_divergence.all_gt(0.0));
        DIALS_ASSERT(delta_mosaicity.all_gt(0.0));
        DIALS_ASSERT(delta_divergence_.size() == delta_mosaicity_.size());
        DIALS_ASSERT(delta_divergence_.size() == scan.get_num_images());
        DIALS_ASSERT(delta_divergence_.size() > 0);
      }

      /**
       * Calculate the bbox on the detector image volume for the reflection.
       *
       * The roi is calculated using the parameters delta_divergence and
       * delta_mosaicity. The reflection mask comprises all pixels where:
       *  |e1| <= delta_d, |e2| <= delta_d, |e3| <= delta_m
       *
       * We transform the coordinates of the box
       *   (-delta_d, -delta_d, 0)
       *   (+delta_d, -delta_d, 0)
       *   (-delta_d, +delta_d, 0)
       *   (+delta_d, +delta_d, 0)
       *
       * to the detector image volume and return the minimum and maximum values
       * for the x, y, z image volume coordinates.
       *
       * @param s1 The diffracted beam vector
       * @param frame The predicted frame number
       * @returns A 6 element array: (minx, maxx, miny, maxy, minz, maxz)
       */
      virtual int6 single(vec3<double> s1, double frame, std::size_t panel) const {
        // Ensure our values are ok
        DIALS_ASSERT(s1.length_sq() > 0);

        // Get the rotation angle
        double phi = scan_.get_angle_from_array_index(frame);

        // Create the coordinate system for the reflection
        CoordinateSystem xcs(m2_, s0_, s1, phi);

        // Get the divergence and mosaicity for this point
        double delta_d = 0.0;
        double delta_m = 0.0;
        if (delta_divergence_.size() == 1) {
          delta_d = delta_divergence_[0];
          delta_m = delta_mosaicity_[0];
        } else {
          int frame0 = scan_.get_array_range()[0];
          int index = (int)std::floor(frame) - frame0;
          if (index < 0) {
            delta_d = delta_divergence_.front();
            delta_m = delta_mosaicity_.front();
          } else if (index >= delta_divergence_.size()) {
            delta_d = delta_divergence_.back();
            delta_m = delta_mosaicity_.back();
          } else {
            delta_d = delta_divergence_[index];
            delta_m = delta_mosaicity_[index];
          }
        }

        // Calculate the beam vectors at the following xds coordinates:
        //   (-delta_d, -delta_d, 0)
        //   (+delta_d, -delta_d, 0)
        //   (-delta_d, +delta_d, 0)
        //   (+delta_d, +delta_d, 0)
        double point = delta_d;
        double3 sdash1 = xcs.to_beam_vector(double2(-point, -point));
        double3 sdash2 = xcs.to_beam_vector(double2(+point, -point));
        double3 sdash3 = xcs.to_beam_vector(double2(-point, +point));
        double3 sdash4 = xcs.to_beam_vector(double2(+point, +point));

        // Get the detector coordinates (px) at the ray intersections
        double2 xy1 = detector_[panel].get_ray_intersection_px(sdash1);
        double2 xy2 = detector_[panel].get_ray_intersection_px(sdash2);
        double2 xy3 = detector_[panel].get_ray_intersection_px(sdash3);
        double2 xy4 = detector_[panel].get_ray_intersection_px(sdash4);

        /// Calculate the rotation angles at the following XDS
        // e3 coordinates: -delta_m, +delta_m
        double phi1 = xcs.to_rotation_angle_fast(-delta_m);
        double phi2 = xcs.to_rotation_angle_fast(+delta_m);

        // Get the array indices at the rotation angles
        double z1 = scan_.get_array_index_from_angle(phi1);
        double z2 = scan_.get_array_index_from_angle(phi2);

        // Return the roi in the following form:
        // (minx, maxx, miny, maxy, minz, maxz)
        // Min's are rounded down to the nearest integer, Max's are rounded up
        double4 x(xy1[0], xy2[0], xy3[0], xy4[0]);
        double4 y(xy1[1], xy2[1], xy3[1], xy4[1]);
        double2 z(z1, z2);
        int6 bbox((int)floor(min(x)),
                  (int)ceil(max(x)),
                  (int)floor(min(y)),
                  (int)ceil(max(y)),
                  (int)floor(min(z)),
                  (int)ceil(max(z)));

        vec2<int> array_range = scan_.get_array_range();
        DIALS_ASSERT(bbox[4] <= frame && frame < bbox[5]);
        bbox[4] = std::max(bbox[4], array_range[0]);
        bbox[4] = std::min(bbox[4], array_range[1] - 1);
        bbox[5] = std::min(bbox[5], array_range[1]);
        bbox[5] = std::max(bbox[5], array_range[0] + 1);
        DIALS_ASSERT(bbox[1] > bbox[0]);
        DIALS_ASSERT(bbox[3] > bbox[2]);
        DIALS_ASSERT(bbox[5] > bbox[4]);
        return bbox;
      }

      /**
       * Calculate the rois for an array of reflections given by the array of
       * diffracted beam vectors and rotation angles.
       * @param s1 The array of diffracted beam vectors
       * @param frame The array of frame numbers.
       */
      virtual af::shared<int6> array(const af::const_ref<vec3<double> >& s1,
                                     const af::const_ref<double>& frame,
                                     const af::const_ref<std::size_t>& panel) const {
        DIALS_ASSERT(s1.size() == frame.size());
        DIALS_ASSERT(s1.size() == panel.size());
        af::shared<int6> result(s1.size(), af::init_functor_null<int6>());
        for (std::size_t i = 0; i < s1.size(); ++i) {
          result[i] = single(s1[i], frame[i], panel[i]);
        }
        return result;
      }

    private:
      vec3<double> s0_;
      vec3<double> m2_;
      Detector detector_;
      Scan scan_;
      af::shared<double> delta_divergence_;
      af::shared<double> delta_mosaicity_;
    };

    /** Calculate the bounding box for each reflection */
    class BBoxCalculator2D : public BBoxCalculatorIface {
    public:
      /**
       * Initialise the bounding box calculation.
       * @param beam The beam parameters
       * @param detector The detector parameters
       * @param delta_divergence The xds delta_divergence parameter
       * @param delta_mosaicity The xds delta_mosaicity parameter
       */
      BBoxCalculator2D(const BeamBase& beam,
                       const Detector& detector,
                       double delta_divergence,
                       double delta_mosaicity)
          : s0_(beam.get_s0()),
            detector_(detector),
            delta_divergence_(delta_divergence) {
        DIALS_ASSERT(delta_divergence > 0.0);
        DIALS_ASSERT(delta_mosaicity >= 0.0);
      }

      /**
       * Calculate the bbox on the detector image volume for the reflection.
       *
       * The roi is calculated using the parameters delta_divergence and
       * delta_mosaicity. The reflection mask comprises all pixels where:
       *  |e1| <= delta_d, |e2| <= delta_d, |e3| <= delta_m
       *
       * We transform the coordinates of the box
       *   (-delta_d, -delta_d, 0)
       *   (+delta_d, -delta_d, 0)
       *   (-delta_d, +delta_d, 0)
       *   (+delta_d, +delta_d, 0)
       *
       * to the detector image volume and return the minimum and maximum values
       * for the x, y, z image volume coordinates.
       *
       * @param s1 The diffracted beam vector
       * @param frame The predicted frame number
       * @returns A 6 element array: (minx, maxx, miny, maxy, minz, maxz)
       */
      virtual int6 single(vec3<double> s1, double frame, std::size_t panel) const {
        // Ensure our values are ok
        DIALS_ASSERT(s1.length_sq() > 0);

        // Create the coordinate system for the reflection
        CoordinateSystem2d xcs(s0_, s1);

        // Get the divergence and mosaicity for this point
        double delta_d = delta_divergence_;

        // Calculate the beam vectors at the following xds coordinates:
        //   (-delta_d, -delta_d, 0)
        //   (+delta_d, -delta_d, 0)
        //   (-delta_d, +delta_d, 0)
        //   (+delta_d, +delta_d, 0)
        double point = delta_d;
        double3 sdash1 = xcs.to_beam_vector(double2(-point, -point));
        double3 sdash2 = xcs.to_beam_vector(double2(+point, -point));
        double3 sdash3 = xcs.to_beam_vector(double2(-point, +point));
        double3 sdash4 = xcs.to_beam_vector(double2(+point, +point));

        // Get the detector coordinates (px) at the ray intersections
        double2 xy1 = detector_[panel].get_ray_intersection_px(sdash1);
        double2 xy2 = detector_[panel].get_ray_intersection_px(sdash2);
        double2 xy3 = detector_[panel].get_ray_intersection_px(sdash3);
        double2 xy4 = detector_[panel].get_ray_intersection_px(sdash4);

        // Return the roi in the following form:
        // (minx, maxx, miny, maxy, minz, maxz)
        // Min's are rounded down to the nearest integer, Max's are rounded up
        double4 x(xy1[0], xy2[0], xy3[0], xy4[0]);
        double4 y(xy1[1], xy2[1], xy3[1], xy4[1]);
        int6 bbox((int)floor(min(x)),
                  (int)ceil(max(x)),
                  (int)floor(min(y)),
                  (int)ceil(max(y)),
                  (int)floor(frame),
                  (int)floor(frame) + 1);
        DIALS_ASSERT(bbox[1] > bbox[0]);
        DIALS_ASSERT(bbox[3] > bbox[2]);
        DIALS_ASSERT(bbox[5] > bbox[4]);
        return bbox;
      }

      /**
       * Calculate the rois for an array of reflections given by the array of
       * diffracted beam vectors and rotation angles.
       * @param s1 The array of diffracted beam vectors
       * @param phi The array of rotation angles.
       */
      virtual af::shared<int6> array(const af::const_ref<vec3<double> >& s1,
                                     const af::const_ref<double>& frame,
                                     const af::const_ref<std::size_t>& panel) const {
        DIALS_ASSERT(s1.size() == frame.size());
        DIALS_ASSERT(s1.size() == panel.size());
        af::shared<int6> result(s1.size(), af::init_functor_null<int6>());
        for (std::size_t i = 0; i < s1.size(); ++i) {
          result[i] = single(s1[i], frame[i], panel[i]);
        }
        return result;
      }

    private:
      vec3<double> s0_;
      Detector detector_;
      double delta_divergence_;
    };

    /**
     * Class to help compute bbox for multiple experiments.
     */
    class BBoxMultiCalculator {
    public:
      /**
       * Add a bbox calculator to the list.
       */
      void push_back(std::shared_ptr<BBoxCalculatorIface> obj) {
        compute_.push_back(obj);
      }

      /**
       * Get the number of calculators.
       */
      std::size_t size() const {
        return compute_.size();
      }

      /**
       * Calculate the rois for an array of reflections given by the array of
       * diffracted beam vectors and rotation angles.
       * @param id The experiment id
       * @param s1 The array of diffracted beam vectors
       * @param phi The array of rotation angles.
       * @param panel The panel number
       */
      af::shared<int6> operator()(const af::const_ref<std::size_t>& id,
                                  const af::const_ref<vec3<double> >& s1,
                                  const af::const_ref<double>& phi,
                                  const af::const_ref<std::size_t>& panel) const {
        DIALS_ASSERT(s1.size() == id.size());
        DIALS_ASSERT(s1.size() == phi.size());
        DIALS_ASSERT(s1.size() == panel.size());
        af::shared<int6> result(s1.size(), af::init_functor_null<int6>());
        for (std::size_t i = 0; i < s1.size(); ++i) {
          DIALS_ASSERT(id[i] < size());
          result[i] = compute_[id[i]]->single(s1[i], phi[i], panel[i]);
        }
        return result;
      }

    private:
      std::vector<std::shared_ptr<BBoxCalculatorIface> > compute_;
    };

    /**
     * Inflate the x and y extents of an array of bounding boxes by a margin,
     * targeting a minimum added total pixel volume per box, which will all
     * be used as background pixels.
     *
     * For each bounding box the margin is chosen as the smallest integer in
     * [margin_min, margin_max] such that the resulting added pixel volume is at
     * least min_added_volume.  If even margin_max is insufficient to reach
     * min_added_volume, the margin is capped at margin_max anyway (the margin
     * bounds always win).
     *
     * Only the x and y extents are inflated; the z (frame) extent is unchanged.
     *
     * @param bbox              Input bounding boxes (minx,maxx,miny,maxy,minz,maxz)
     * @param margin_min        Minimum number of pixels to add to each side (default 3)
     * @param margin_max        Maximum number of pixels to add to each side (default
     * 16)
     * @param min_added_volume  Target minimum total pixel volume (default 4096)
     * @returns                 New array of inflated bounding boxes
     */
    af::shared<int6> inflate_bboxes(const af::const_ref<int6>& bbox,
                                    int margin_min = 3,
                                    int margin_max = 16,
                                    int min_added_volume = 3600) {
      DIALS_ASSERT(margin_min >= 0);
      DIALS_ASSERT(margin_max >= margin_min);
      DIALS_ASSERT(min_added_volume > 0);

      af::shared<int6> result(bbox.size(), af::init_functor_null<int6>());

      for (std::size_t i = 0; i < bbox.size(); ++i) {
        const int6& b = bbox[i];

        // Current extents
        int nx = b[1] - b[0];  // width  in x
        int ny = b[3] - b[2];  // height in y
        int nz = b[5] - b[4];  // depth  in z (unchanged)

        DIALS_ASSERT(nx > 0 && ny > 0 && nz > 0);

        // Find the smallest margin in [margin_min, margin_max] that satisfies
        // the target volume.
        int target_volume = nx * ny * nz + min_added_volume;
        int margin = margin_min;
        for (int m = margin_min; m <= margin_max; ++m) {
          int inflated_volume = (nx + 2 * m) * (ny + 2 * m) * nz;
          if (inflated_volume >= target_volume) {
            margin = m;
            break;
          }
          // If we exhaust the range without hitting target_volume, margin_max is used.
          margin = m;  // keeps updating so we end on margin_max if needed
        }

        result[i] =
          int6(b[0] - margin, b[1] + margin, b[2] - margin, b[3] + margin, b[4], b[5]);

        DIALS_ASSERT(result[i][1] > result[i][0]);
        DIALS_ASSERT(result[i][3] > result[i][2]);
        DIALS_ASSERT(result[i][5] > result[i][4]);
      }

      return result;
    }

}}}}  // namespace dials::algorithms::profile_model::gaussian_rs

#endif  // DIALS_ALGORITHMS_PROFILE_MODEL_GAUSSIAN_RS_BBOX_CALCULATOR_H
