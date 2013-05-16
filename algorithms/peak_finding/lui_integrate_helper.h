#ifndef DIALS_ALGORITHMS_INTEGRATE_LUI_HELPER_H
#define DIALS_ALGORITHMS_INTEGRATE_LUI_HELPER_H
#include <iostream>
#include <cmath>
#include <scitbx/array_family/flex_types.h>
const float pi=3.14159265358;
namespace dials { namespace algorithms {
  using scitbx::af::flex_int;
  using scitbx::af::flex_grid;
  void rotate (float& x, float& y, float delta_ang){
    float ang,dist;
    ang = atan2(x, y);
    dist = sqrt(x * x + y * y);
    ang = ang + delta_ang * pi;
    x = dist * sin(ang);
    y = dist * cos(ang);
  }

  flex_int model_2d(int nrow, int ncol, float a, float b,
                float delta_ang, float imax, float asp) {
    float dx,dy,dd,xc,yc;
    float mw = 0.5;
    float cntnt, gss, lrz, i_tt;
    // float tot=0;
    xc=float(ncol)/2;
    yc=float(nrow)/2;
    std::cout << "\n" << "xc =" << xc << ", yc =" << yc << "\n";
    // flex_int curv3d(data2d.accessor(),0);
    flex_int curv3d(flex_grid<>(nrow, ncol),0);
    cntnt=1.0 / sqrt(2.0 * pi);
    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        dx = float(col) - xc;
        dy = float(row) - yc;
        rotate(dx, dy, delta_ang);
        dx = dx / a;
        dy = dy / b;
        dd = sqrt(dx*dx + dy*dy);
        if (dd < 0.000000000000000000001){
          dd=0.000000000000000000001;
        }
        gss = cntnt * exp(-mw * (dd*dd));
        lrz = 1.0 /(pi * (1.0 + dd*dd ));
        i_tt=gss*asp+lrz*(1.0 - asp);
        // curv3d(row, col) = data2d(row, col) + int(i_tt*imax);
        curv3d(row, col) = int(i_tt*imax);
        // tot+=i_tt*imax;
      }
    }
    // std::cout <<"\n"<<"tot ="<<tot<<"\n";
    return curv3d;
  }

  float measure_2d_angl(flex_int & data2d, float xpos, float ypos) {
    int ncol=data2d.accessor().all()[1];
    int nrow=data2d.accessor().all()[0];
    float dx, dy, xfl, yfl;
    float r_ang, r_avg_ang, r_ang_tot = 0, r_tot = 0;
    float l_ang, l_avg_ang, l_ang_tot = 0, l_tot = 0;
    float avg_ang;
    // std::cout <<"\n" << "x,y =" << xpos << ", " << ypos << "\n";
    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        xfl = float(col);// + 0.5;
        yfl = float(row);// + 0.5;
        dx = xfl - xpos;
        dy = yfl - ypos;
        if( dx < 0 ){
          r_ang = atan2(dx, dy);
          r_ang_tot += r_ang * float(data2d(row, col));
          r_tot += data2d(row, col);
        }else{
          l_ang = atan2(dx, dy);
          l_ang_tot += l_ang * float(data2d(row, col));
          l_tot += data2d(row, col);
        }
      }
    }
    r_avg_ang = (r_ang_tot/r_tot)/pi;
    l_avg_ang = (l_ang_tot/l_tot)/pi;
    avg_ang=(r_avg_ang+(l_avg_ang-pi))/2.0;
    return avg_ang;
  }

}}

#endif
