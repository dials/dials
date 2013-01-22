#ifndef DIALS_SPOT_PREDICTION_COORDINATE_FILTER_H
#define DIALS_SPOT_PREDICTION_COORDINATE_FILTER_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include "../error.h"

namespace dials { namespace spot_prediction {

typedef scitbx::af::flex <scitbx::vec2 <double> >::type flex_vec2_double;
typedef scitbx::af::flex <scitbx::vec3 <double> >::type flex_vec3_double;

scitbx::af::flex_bool in_range(const scitbx::af::flex_double &x, 
                               scitbx::vec2 <double> range) {
    scitbx::af::flex_bool result(x.size());
    for (int i = 0; i < x.size(); ++i) {
        result[i] = (range[0] <= x[i] && x[i] <= range[1]) ? 1 : 0;
    }
    return result;
}

scitbx::af::flex_bool in_rect(const flex_vec2_double &xy, 
                              scitbx::vec2 <double> xrange, 
                              scitbx::vec2 <double> yrange) {
    scitbx::af::flex_bool result(xy.size());
    for (int i = 0; i < xy.size(); ++i) {
        result[i] = ((xrange[0] <= xy[i][0] && xy[i][0] <= xrange[1]) &&
                     (yrange[0] <= xy[i][1] && xy[i][1] <= yrange[1])) ? 1 : 0;
    }
    return result;
}

scitbx::af::flex_bool in_volume(const flex_vec3_double &xyz, 
                                scitbx::vec2 <double> xrange, 
                                scitbx::vec2 <double> yrange,
                                scitbx::vec2 <double> zrange) {
    scitbx::af::flex_bool result(xyz.size());
    for (int i = 0; i < xyz.size(); ++i) {
        result[i] = ((xrange[0] <= xyz[i][0] && xyz[i][0] <= xrange[1]) &&
                     (yrange[0] <= xyz[i][1] && xyz[i][1] <= yrange[1]) &&
                     (zrange[0] <= xyz[i][2] && xyz[i][2] <= zrange[1])) ? 1 : 0;
    }
    return result;
}

template <typename ArrayType>
ArrayType remove_if(const ArrayType &in, const scitbx::af::flex_bool &yesno) {
    DIALS_ASSERT(in.size() == yesno.size());
    int size = 0;
    for (int i = 0; i < yesno.size(); ++i) {
        if (!yesno[i]) size++;
    }
    ArrayType result(size);
    for (int i = 0, j = 0; i < in.size(); ++i) {
        if (!yesno[i]) {
            result[j++] = in[i];
        }
    }
    return result;
}

template <typename ArrayType>
ArrayType remove_if_not(const ArrayType &in, const scitbx::af::flex_bool &yesno) {
    DIALS_ASSERT(in.size() == yesno.size());
    int size = 0;
    for (int i = 0; i < yesno.size(); ++i) {
        if (yesno[i]) size++;
    }
    ArrayType result(size);
    for (int i = 0, j = 0; i < in.size(); ++i) {
        if (yesno[i]) {
            result[j++] = in[i];
        }
    }
    return result;
}

}} // namespace dials::spot_prediction

#endif // DIALS_SPOT_PREDICTION_COORDINATE_FILTER_H
