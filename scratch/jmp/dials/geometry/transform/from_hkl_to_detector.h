
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_DETECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_DETECTOR_H

#include "from_hkl_to_beam_vector.h"
#include "from_beam_vector_to_detector.h"

namespace dials { namespace geometry { namespace transform {

/** A transform from miller indices to detector coordinates */
class FromHklToDetector {

public:

    /** Default constructor */
    FromHklToDetector() {}

    /**
     * Initialise the transform using component transform objects.
     * @param hkl_to_s1 The hkl to beam vector transform
     * @param s1_to_xy The beam vector to detector transform
     */
    FromHklToDetector(FromHklToBeamVector hkl_to_s1,
                      FromBeamVectorToDetector s1_to_xy)
        : _hkl_to_s1(hkl_to_s1),
          _s1_to_xy(s1_to_xy) {}
     
    /**
     * Initialise the transform using component transform objects.
     * @param rlcs The reciprocal lattice coordinate system class
     * @param s0 The incident beam vector
     * @param m2 The rotation axis
     * @param dcs The detector coordinate system
     * @param origin The origin of the detector coordinate system
     * @param distance The distance from the detector to the crystal     
     */        
    FromHklToDetector(ReciprocalLatticeCoordinateSystem rlcs,
                      scitbx::vec3 <double> s0,
                      scitbx::vec3 <double> m2,
                      DetectorCoordinateSystem dcs,
                      scitbx::vec2 <double> origin,
                      double distance)
        : _hkl_to_s1(rlcs, s0, m2),
          _s1_to_xy(dcs, origin, distance) {}

public:

    /**
     * Apply the transform to the miller indices and rotation angle.
     * @param hkl The miller indices
     * @param phi The rotation angle
     * @returns The detector coordinates
     */
    scitbx::vec2 <double> apply(scitbx::vec3 <int> hkl, double phi) {
        return _s1_to_xy.apply(_hkl_to_s1.apply(hkl, phi));
    }

private:

    FromHklToBeamVector _hkl_to_s1;
    FromBeamVectorToDetector _s1_to_xy;
};

}}} // namespace = dials::geometry::transform

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_DETECTOR_H
