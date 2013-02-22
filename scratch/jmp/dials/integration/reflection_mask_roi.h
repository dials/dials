
#ifndef DIALS_INTEGRATION_REFLECTION_MASK_ROI_CALCULATOR_H
#define DIALS_INTEGRATION_REFLECTION_MASK_ROI_CALCULATOR_H

#include <cmath>
#include <scitbx/constants.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/ref_reductions.h>
#include "../equipment/beam.h"
#include "../equipment/goniometer.h"
#include "../equipment/detector.h"
#include "../geometry/xds_coordinate_system.h"
#include "../geometry/transform/from_xds_to_detector.h"
#include "../geometry/transform/from_xds_e3_to_phi.h"
#include "../array_family/array_types.h"
#include "../reflection/reflection.h"

namespace dials { namespace integration {

/** Calculate the reflection mask roi for each reflection */
class ReflectionMaskRoi {

public:

    /** Default constructor */
    ReflectionMaskRoi() {}

    /**
     * Initialise the reflection mask roi calculation.
     * @param beam The beam parameters
     * @param detector The detector parameters
     * @param goniometer The goniometer parameters
     * @param delta_divergence The xds delta_divergence parameter
     * @param delta_mosaicity The xds delta_mosaicity parameter
     */
    ReflectionMaskRoi(const equipment::Beam &beam,
                      const equipment::Detector &detector,
                      const equipment::Goniometer &goniometer,
                      double delta_divergence,
                      double delta_mosaicity)
        : beam_(beam),
          detector_(detector),
          goniometer_(goniometer),
          delta_divergence_(delta_divergence),
          delta_mosaicity_(delta_mosaicity) {}

    /**
     * Calculate the roi on the detector image volume for the reflection.
     *
     * The roi is calculated using the parameters delta_divergence and
     * delta_mosaicity. The reflection mask comprises all pixels where:
     *  |e1| <= delta_d / 2, |e2| <= delta_d < 2, |e3| <= delta_m / 2
     *
     * We transform the coordinates of the box
     *   (-delta_d/2, -delta_d/2, 0)
     *   (+delta_d/2, -delta_d/2, 0)
     *   (-delta_d/2, +delta_d/2, 0)
     *   (+delta_d/2, +delta_d/2, 0)
     *
     * to the detector image volume and return the minimum and maximum values
     * for the x, y, z image volume coordinates.
     *
     * @todo Figure out why scitbx::af::tiny <int, 6> doesn't convert to a
     * python type but scitbx::af::tiny <double, 6> does.
     *
     * @param s1 The diffracted beam vector
     * @param phi The rotation angle
     * @returns A 6 element array: (minx, maxx, miny, maxy, minz, maxz)
     */
    scitbx::af::tiny <int, 6> calculate(scitbx::vec3 <double> s1, double phi) {

        // Create the coordinate system for the reflection
        geometry::XdsCoordinateSystem xcs(beam_.get_direction(),
                                          s1,
                                          goniometer_.get_rotation_axis(),
                                          phi);

        // Create the transformer from the xds coordinate system to the detector
        geometry::transform::FromXdsToDetector from_xds_to_detector(
                                                xcs, s1, detector_);

        // Create the transformer from Xds E3 to rotation angle
        geometry::transform::FromXdsE3ToPhi from_xds_e3_to_phi(
            xcs.get_zeta(), phi);

        // Find the detector coordinates at the following xds coordinates:
        //   (-delta_d/2, -delta_d/2, 0)
        //   (+delta_d/2, -delta_d/2, 0)
        //   (-delta_d/2, +delta_d/2, 0)
        //   (+delta_d/2, +delta_d/2, 0)
        double point = delta_divergence_ / 2.0;
        scitbx::vec2 <double> xy1 = from_xds_to_detector.apply(
            scitbx::vec3 <double> (-point, -point, 0.0));
        scitbx::vec2 <double> xy2 = from_xds_to_detector.apply(
            scitbx::vec3 <double> (+point, -point, 0.0));
        scitbx::vec2 <double> xy3 = from_xds_to_detector.apply(
            scitbx::vec3 <double> (-point, +point, 0.0));
        scitbx::vec2 <double> xy4 = from_xds_to_detector.apply(
            scitbx::vec3 <double> (+point, +point, 0.0));

        /// Get the image volume z coordinates (zero based) at the following XDS
        // e3 coordinates: -delta_m/2, +delta_m/2
        double z1 = goniometer_.get_zero_based_frame_from_angle(
                        from_xds_e3_to_phi.apply(-delta_mosaicity_ / 2.0),
                        false);
        double z2 = goniometer_.get_zero_based_frame_from_angle(
                        from_xds_e3_to_phi.apply(+delta_mosaicity_ / 2.0),
                        false);

        // Return the roi in the following form:
        // (minx, maxx, miny, maxy, minz, maxz)
        // Min's are rounded down to the nearest integer, Max's are rounded up
        scitbx::af::tiny <double, 4> x(xy1[0], xy2[0], xy3[0], xy4[0]);
        scitbx::af::tiny <double, 4> y(xy1[1], xy2[1], xy3[1], xy4[1]);
        scitbx::af::tiny <double, 2> z(z1, z2);
        return scitbx::af::tiny <int, 6> (
                (int)std::floor(scitbx::af::min(x)),
                (int)std::ceil (scitbx::af::max(x)),
                (int)std::floor(scitbx::af::min(y)),
                (int)std::ceil (scitbx::af::max(y)),
                (int)std::floor(scitbx::af::min(z)),
                (int)std::ceil (scitbx::af::max(z)));
    }

    /**
     * Calculate the rois for an array of reflections given by the array of
     * diffracted beam vectors and rotation angles.
     * @param s1 The array of diffracted beam vectors
     * @param phi The array of rotation angles.
     */
    af::flex_tiny6_int calculate(const af::flex_vec3_double &s1,
                                 const scitbx::af::flex_double &phi) {
        DIALS_ASSERT(s1.size() == phi.size());
        af::flex_tiny6_int result(s1.size());
        for (int i = 0; i < s1.size(); ++i) {
            result[i] = calculate(s1[i], phi[i]);
        }
        return result;
    }

    /**
     * Calculate the rois for a reflection given by reflection struct
     * @param reflection The reflection data
     */
    void calculate(Reflection &reflection) {
        reflection.set_region_of_interest(
            calculate(reflection.get_beam_vector(),
                      reflection.get_rotation_angle()));
    }

    /**
     * Calculate the rois for an array of reflections given by the array of
     * reflections
     * @param reflections The list of reflections
     */
    void calculate(ReflectionList &reflections) {
        for (int i = 0; i < reflections.size(); ++i) {
            calculate(reflections[i]);
        }
    }

private:

    equipment::Beam beam_;
    equipment::Detector detector_;
    equipment::Goniometer goniometer_;
    double delta_divergence_;
    double delta_mosaicity_;
};

}} // namespace dials::integration

#endif // DIALS_INTEGRATION_REFLECTION_MASK_ROI_CALCULATOR_H
