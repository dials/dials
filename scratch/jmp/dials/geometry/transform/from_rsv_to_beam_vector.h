
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_RSV_TO_BEAM_VECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_RSV_TO_BEAM_VECTOR_H

#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include "../../equipment/beam.h"
#include "../../equipment/goniometer.h"
#include "../../error.h"

namespace dials { namespace geometry { namespace transform {

typedef scitbx::af::flex <scitbx::vec3 <double> >::type flex_vec3_double;

/** Transform from reciprocal space vector to diffracted beam vector */
class FromRsvToBeamVector {

public:

    /** Default constructor */
    FromRsvToBeamVector() {}

    /** 
     * Initialise the transform.
     * @param beam The beam parameters
     * @param goniometer The goniometer parameters
     */
    FromRsvToBeamVector(equipment::Beam beam,
                        equipment::Goniometer goniometer)
        : incident_beam_vector_(beam.get_direction().normalize() * 
                                beam.get_wavelength()),
          rotation_axis_(goniometer.get_rotation_axis().normalize()) {}

    /**
     * Apply the transform
     * @param reciprocal_space_vector The reciprocal space vector
     * @param rotation_angle The diffracting rotation angle
     */
    scitbx::vec3 <double> apply(scitbx::vec3 <double> reciprocal_space_vector, 
                                double rotation_angle) {
        return reciprocal_space_vector.unit_rotate_around_origin(
                rotation_axis_, rotation_angle) + incident_beam_vector_;
    }

    /**
     * Apply the transform
     * @param reciprocal_space_vectors The reciprocal space vectors
     * @param rotation_angles The diffracting rotation angles
     */
    flex_vec3_double apply(flex_vec3_double reciprocal_space_vectors,
                           scitbx::af::flex_double rotation_angles) {
        DIALS_ASSERT(reciprocal_space_vectors.size() == rotation_angles.size());
        flex_vec3_double result(reciprocal_space_vectors.size());
        for (int i = 0; i < reciprocal_space_vectors.size(); ++i) {
            result[i] = apply(reciprocal_space_vectors[i],
                              rotation_angles[i]);
        }
        return result;
    }

private:

    scitbx::vec3 <double> incident_beam_vector_;
    scitbx::vec3 <double> rotation_axis_;
};

}}} // namespace dials::geometry:: transform

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_RSV_TO_BEAM_VECTOR_H