
#ifndef DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_BEAM_VECTOR_H
#define DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_BEAM_VECTOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>

namespace dials { namespace geometry { namespace transform {

class from_hkl_to_beam_vector {

public:

    from_hkl_to_beam_vector() {}

    from_hkl_to_beam_vector(scitbx::vec3 <double> b1_star,
                            scitbx::vec3 <double> b2_star,
                            scitbx::vec3 <double> b3_star,
                            scitbx::vec3 <double> s0,
                            scitbx::vec3 <double> m2) 
        : _b1_star(b1_star),
          _b2_star(b2_star),
          _b3_star(b3_star),
          _s0(s0),
          _m2(m2.normalize()) {}

public:

    scitbx::vec3 <double> apply(scitbx::vec3 <int> hkl, double phi) {
        return _s0 + (double(hkl[0]) * _b1_star + 
                      double(hkl[1]) * _b2_star + 
                      double(hkl[2]) * _b3_star)
                        .unit_rotate_around_origin(_m2, phi);
    }

private:

    scitbx::vec3 <double> _b1_star;
    scitbx::vec3 <double> _b2_star;
    scitbx::vec3 <double> _b3_star;
    scitbx::vec3 <double> _s0;
    scitbx::vec3 <double> _m2;
};

}}}

#endif // DIALS_GEOMETRY_TRANSFORM_FROM_HKL_TO_BEAM_VECTOR_H
