
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_DETECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_DETECTOR_H

#include "from_hkl_to_beam_vector.h"
#include "from_beam_vector_to_detector.h"

namespace dials { namespace geometry { namespace transform {

class from_hkl_to_detector {

public:

    from_hkl_to_detector() {}

    from_hkl_to_detector(from_hkl_to_beam_vector hkl_to_s1,
                         from_beam_vector_to_detector s1_to_xy)
        : _hkl_to_s1(hkl_to_s1),
          _s1_to_xy(s1_to_xy) {}

public:

    scitbx::vec2 <double> apply(scitbx::vec3 <int> hkl, double phi) {
        return _s1_to_xy.apply(_hkl_to_s1.apply(hkl, phi));
    }

private:

    from_hkl_to_beam_vector _hkl_to_s1;
    from_beam_vector_to_detector _s1_to_xy;
};

}}}

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_BEAM_VECTOR_H
