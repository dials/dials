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

  float measure_2d_angl(flex_int & data2d, flex_int & mask2d,float xpos, float ypos) {
    int ncol=data2d.accessor().all()[1];
    int nrow=data2d.accessor().all()[0];
    std::cout <<"\n" << "ncol =" << ncol <<"\n" << "nrow =" << nrow <<"\n";
    float dx, dy;
    int pie_size;
    // std::cout <<"\n" << "x,y =" << xpos << ", " << ypos << "\n";
    int contr = 0;
    float pi_sqrt = 1.77245385091;

    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        if( mask2d(row, col)==1 ){
          contr++;
        }
      }
    }

    std::cout <<"\n" << "cont =" << contr <<"\n";
    float angl_tmp[int(contr)];
    float dist_tmp[int(contr)];
    contr = 0;

    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        if( mask2d(row, col)==1 ){
          contr++;
          dx = float(col) - xpos;
          dy = float(row) - ypos;
          angl_tmp[contr] = atan2(dx, dy);
          dist_tmp[contr] = sqrt(dx * dx + dy * dy);

        }
      }
    }
    int area=contr;
    pie_size=int(2.0 * pi_sqrt * sqrt(float(area)));
    std::cout << "pie_size =" << pie_size << "\n";
    /*float dist[int(pie_size)];
    float angl[int(pie_size)];*/
    float ang_min, ang_max;
    for (int pos=0; pos < pie_size; pos++){
      ang_min = 2.0 * pi * float(pos) / float(pie_size);
      ang_max = 2.0 * pi * float(pos + 1) / float(pie_size);
   /*   if (ang_max > 2.0*pi){
        ang_max = 2.0*pi;
      }
      if (ang_max < 0.0){
        ang_max = 0.0;
      } */

      std::cout << "pos =" << pos << "\n";

      std::cout << "\n" << "ang min =" << ang_min << "\n" << "ang max =" << ang_max << "\n";
    }
    std::cout << "\n" << "2 * pi =" << 2.0 * pi << "\n";

    return pie_size;
  }

}}

#endif
