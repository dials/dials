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

  void get_polar (float& ang, float& dst, float x, float y){
    ang = atan2(x, y);
    if(ang<0){
      ang+=pi*2.0;
    }
    dst = sqrt(x * x + y * y);
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

    //std::cout <<"\n" << "cont =" << contr <<"\n";
    float angl_tmp[int(contr)];
    float dist_tmp[int(contr)];
    int   xcol_tmp[int(contr)];
    int   yrow_tmp[int(contr)];
    int area=contr;
    contr = 0;
    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        if( mask2d(row, col)==1 ){
          contr++;
          dx = float(col) - xpos;
          dy = float(row) - ypos;
          get_polar (angl_tmp[contr], dist_tmp[contr], dx, dy);
          xcol_tmp[contr]=col;
          yrow_tmp[contr]=row;
        }
      }
    }
/*
    for (int i=0;i<contr;i++){
      std::cout << "\n" << "ang ="<<angl_tmp[i];
    }
*/

    pie_size=int(2.0 * pi_sqrt * sqrt(float(area)));
    //std::cout << "pie_size =" << pie_size << "\n";
    float dist[int(pie_size)];
    float angl[int(pie_size)];
    int   xcol[int(pie_size)];
    int   yrow[int(pie_size)];

    float ang_min, ang_max, max_dst;
    for (int pos=0; pos < pie_size; pos++){
      ang_min = 2.0 * pi * float(pos) / float(pie_size);
      ang_max = 2.0 * pi * float(pos + 1) / float(pie_size);
      if (ang_max > 2.0*pi){
        ang_max = 2.0*pi;
      }
      if (ang_max < 0.0){
        ang_max = 0.0;
      }

      dist[pos] = -1.0;
      angl[pos] = -1.0;
      max_dst=0;
      for (int locl_pos = 0; locl_pos < area; locl_pos++){
        if(angl_tmp[locl_pos] > ang_min and angl_tmp[locl_pos] < ang_max){
          if(dist_tmp[locl_pos] > max_dst){
            max_dst = dist_tmp[locl_pos];
            dist[pos] = dist_tmp[locl_pos];
            angl[pos] = angl_tmp[locl_pos];
            xcol[pos] = xcol_tmp[locl_pos];
            yrow[pos] = yrow_tmp[locl_pos];
          }
        }
      }
    }

    int prev_pos, next_pos;
    //float avg_row, avg_col;
    for (int pos=0; pos < pie_size; pos++){
      if (angl[pos] == -1 and dist[pos] == -1){
        prev_pos=pos-1;
        next_pos=pos+1;
        if(prev_pos<0){
          prev_pos+=pie_size;
        }
        if(next_pos>=pie_size){
          next_pos-=pie_size;
        }
        if (angl[next_pos]==-1 and dist[next_pos]==-1
         and angl[prev_pos]==-1 and dist[prev_pos]==-1){
           angl[pos]=float(pos)/pie_size;
        }
        else if (angl[next_pos]!=-1 and dist[next_pos]!=-1
            and angl[prev_pos]!=-1 and dist[prev_pos]!=-1){
                                                                                  /* starts block to be fixed
          avg_col=float(xcol_tmp[prev_pos]+xcol_tmp[next_pos])/2.0;
          avg_row=float(yrow_tmp[prev_pos]+yrow_tmp[next_pos])/2.0;

          dx = float(avg_col) - xpos;
          dy = float(avg_row) - ypos;
          get_polar (angl[pos], dist[pos], dx, dy);
                                                                                   ends block to be fixed */
          angl[pos]=(angl[prev_pos]+angl[next_pos])/2.0;
          dist[pos]=(dist[prev_pos]+dist[next_pos])/2.0;

        }
      }
    }
/*
    for (int pos=0; pos < pie_size; pos++){
      std::cout <<"pos ="<<pos<<"\n________________________________________\n";
      std::cout <<"ang ="<<angl[pos]<< "\n"<<"dst ="<<dist[pos]<<"\n";
    }
*/
    int pos2,pos3,pos4, fin_pos = -1;
    float avg_1, avg_2, dd, fin_dd = 0.0;
    for (int pos1=0; pos1 < pie_size / 2; pos1++){
      pos2 = pos1 + pie_size / 4;
      pos3 = pos1 + pie_size / 2;
      pos4 = pos1 + pie_size * 3 / 4;
      if (pos3 > pie_size){
        pos3 -= pie_size;
      }
      if (pos4 > pie_size){
        pos4 -= pie_size;
      }
      avg_1=(dist[pos1]+dist[pos3])/2.0;
      avg_2=(dist[pos2]+dist[pos4])/2.0;
      dd=avg_1-avg_2;
      if (dd > fin_dd){
        fin_dd = dd;
        fin_pos = pos1;
      }
    }
    float fin_ang = angl[fin_pos];

    return fin_ang;
  }

}}

#endif
