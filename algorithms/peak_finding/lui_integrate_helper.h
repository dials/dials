#ifndef DIALS_ALGORITHMS_INTEGRATE_LUI_HELPER_H
#define DIALS_ALGORITHMS_INTEGRATE_LUI_HELPER_H
#include <iostream>
#include <cmath>
#include <scitbx/array_family/flex_types.h>
const float pi=3.141595358;
namespace dials { namespace algorithms {
  using scitbx::af::flex_int;
  void rotate (float& x, float& y, float delta_ang){
    float ang,dist;
    ang = atan2(x, y);
    dist = sqrt(x * x + y * y);
    ang = ang + delta_ang * pi;
    x = dist * sin(ang);
    y = dist * cos(ang);
  }

  flex_int ref_2d(flex_int & data2d, float a, float b,
                       float delta_ang, float imax) {
    int ncol=data2d.accessor().all()[1];
    int nrow=data2d.accessor().all()[0];
    float dx,dy,dd,xc,yc,l;
    float mw = 0.5;
    xc=ncol/2;
    yc=nrow/2;
    flex_int curv3d(data2d.accessor(),0);
    std::cout << "size(x) =" << ncol <<"\n";
    std::cout << "size(y) =" << nrow <<"\n";
    for (int row = 0; row<nrow; row++) {
      for (int col = 0; col<ncol; col++) {
        dx = col - xc;
        dy = row - yc;
        rotate(dx, dy, delta_ang);
        dx = dx / a;
        dy = dy / b;
        dd = sqrt(dx*dx + dy*dy);
        if (dd < 0.0000000000000001){
          dd=0.0000000000000001;
        }
        l=(1.0/pi) * ( mw / ( dd*dd + mw*mw ) );
        curv3d(row, col) = data2d(row, col) + int(l*imax);
      }
    }
    return curv3d;
  }

}}

#endif
