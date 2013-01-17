
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_XDS_TO_DETECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_XDS_TO_DETECTOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include "from_xds_to_beam_vector.h"
#include "from_beam_vector_to_detector.h"

namespace dials { namespace geometry { namespace transform {

/** A class representing transform from XDS to detector coordinate */
class FromXdsToDetector {

public:

    /** The default constructor */
    FromXdsToDetector() {}

    /**
     * Initialise the transform
     * @param from_xds_to_beam_vector The XDS -> Beam vector transform
     * @param from_beam_vector_to_detector The Beam vector to Detector transform
     */
    FromXdsToDetector(FromXdsToBeamVector from_xds_to_beam_vector,
                      FromBeamVectorToDetector from_beam_vector_to_detector)
        : from_xds_to_beam_vector_(from_xds_to_beam_vector),
          from_beam_vector_to_detector_(from_beam_vector_to_detector) {}
          
    /**
     * Initialise the transform
     * @param xcs The XDS coordinate system
     * @param s1 The diffracted beam vector
     * @param dcs The detector coordinate system
     * @param pixel_size The detector pixel size in mm
     * @param origin The detector origin
     * @param distance The detector distance
     */
    FromXdsToDetector(XdsCoordinateSystem xcs,
                      scitbx::vec3 <double> s1,
                      DetectorCoordinateSystem dcs,
                      scitbx::vec2 <double> pixel_size,                      
                      scitbx::vec2 <double> origin,
                      double distance)
        : from_xds_to_beam_vector_(xcs, s1),
          from_beam_vector_to_detector_(dcs, pixel_size, origin, distance) {}

    /**
     * Apply the transform
     * @param c The xds coordinate
     * @returns The detector coordinate
     */
    scitbx::vec2 <double> apply(scitbx::vec3 <double> c) {
        return from_beam_vector_to_detector_.apply(
                    from_xds_to_beam_vector_.apply(c));
    }

private:

    FromXdsToBeamVector from_xds_to_beam_vector_;
    FromBeamVectorToDetector from_beam_vector_to_detector_;
};

}}}

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_XDS_TO_DETECTOR_H
