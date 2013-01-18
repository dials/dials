
#ifndef DIALS_SPOT_PREDICTION_ANGLE_FILTER_H
#define DIALS_SPOT_PREDICTION_ANGLE_FILTER_H

#include <scitbx/constants.h>
#include <scitbx/vec2.h>

namespace dials { namespace spot_prediction { 

/** Convert the angle mod 2PI */
inline 
double mod_2pi(double angle) {
    return angle - scitbx::constants::two_pi * 
        std::floor(angle / scitbx::constants::two_pi); 
}

/** Convert a pair of angles to mod 2PI */
inline 
scitbx::vec2 <double> mod_2pi(scitbx::vec2 <double> angles) {
	return scitbx::vec2 <double> (mod_2pi(angles[0]), mod_2pi(angles[1]));
}

/** Convert and angle (degrees or radians) to mod 2PI radians */
inline 
double mod_2pi_radians(double angle, bool deg) {
    return mod_2pi(deg ? scitbx::deg_as_rad(angle) : angle);
}

/** Convert and angle (degrees or radians) to mod 2PI radians */
inline 
scitbx::vec2 <double> mod_2pi_radians(scitbx::vec2 <double> angle, bool deg) {
    return scitbx::vec2 <double> (
        mod_2pi_radians(angle[0], deg), 
        mod_2pi_radians(angle[1], deg));
}

/**
    * Check if the angle is within the filter range
    * @param angle The angle to check
    * @param range The angular range
    * @param def In degrees (True/False)
    * @returns True/False the angle is in the filter range
    */
inline 
bool angle_filter(double angle, scitbx::vec2 <double> range, bool deg = false) {
    return (mod_2pi_radians(angle - range[1], deg) >= 
            mod_2pi_radians(angle - range[0], deg) || 
            mod_2pi_radians(angle - range[1], deg) == 0 || 
            mod_2pi_radians(range[0] - range[1], deg) == 0);
}

}} // namespace dials::spot_prediction

#endif // DIALS_SPOT_PREDICTION_ANGLE_FILTER_H