#ifndef DIALS_ALGORITHMS_INTEGRATE_LUI_HELPER_H
#define DIALS_ALGORITHMS_INTEGRATE_LUI_HELPER_H
#include <iostream>
#include <cmath>
const float pi=3.14159265358;
namespace dials { namespace algorithms {

  void get_polar (float& ang, float& dst, float x, float y){
    ang = atan2(x, y);
    if(ang<0){
      ang+=pi*2.0;
    }
    dst = sqrt(x * x + y * y);
  }

}}

#endif
