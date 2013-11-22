#ifndef DIALS_SCRATCH_LUI_HELPER_H
#define DIALS_SCRATCH_LUI_HELPER_H
#include <iostream>
#include <scitbx/vec2.h>
#include <scitbx/array_family/flex_types.h>

#include <scitbx/array_family/versa_matrix.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <scitbx/matrix/inversion.h>


#include <cstdlib>
#include <scitbx/array_family/accessors/mat_grid.h>


#include <cmath>
#include <stdio.h>
#include <cstdlib>

const float pi=3.14159265358;

namespace dials { namespace scratch {
  using scitbx::af::flex_int;
  using scitbx::af::flex_double;
  using scitbx::af::flex_grid;
  using scitbx::vec2;
  int hello_tst() {
    std::cout << "Hello world \n";
    int a=5;
    return a;
  }

  af::versa< double, af::c_grid<2> >  tst_ref_prod(
         af::const_ref<double, af::c_grid<2> > matr01,
         af::ref<double, af::c_grid<2> > matr02)
  {
    af::versa< double, af::c_grid<2> > prod(matr01.accessor());
    std::cout << "Hello from  tst_ref_prod\n";

    matr02(0,0) = 555;
    matr02(1,1) = 666;
    matr02(2,2) = 777;

    prod = af::matrix_multiply(matr01, matr02);

    return prod;
  }

  af::versa< double, af::c_grid<2> >  tst_prod(
         af::const_ref<double, af::c_grid<2> > matr01,
         af::const_ref<double, af::c_grid<2> > matr02)
  {
    af::versa< double, af::c_grid<2> > prod(matr01.accessor());
    std::cout << "Hello from prod_tst \n";
//    af::versa< double, af::c_grid<2> > prod(matr01.accessor());
  //  af::versa< double, af::c_grid<2> > dat_in(matr01.accessor());


    prod = af::matrix_multiply(matr01, matr02);

    return prod;
  }
/*
    // example code that adds the contents of two versa arrays
    for (std::size_t i = 0; i < prod.size(); ++i) {
      prod[i] = matr01[i] + matr02[i];
    }

    // attempt to do inversion of matrix, this is not working

      af::versa< double, af::c_grid<2> >  tst_prod(
         af::const_ref<double, af::c_grid<2> > matr01,
         af::const_ref<double, af::c_grid<2> > matr02)
  {
    std::cout << "Hello from prod_tst \n";
    af::versa< double, af::c_grid<2> > prod(matr01.accessor());
    af::versa< double, af::c_grid<2> > dat_in(matr01.accessor());
    // copying the content of a versa array
        for (std::size_t i = 0; i < dat_in.size(); ++i) {
          dat_in[i] = matr01[i];
        }

    af::ref<double, af::c_grid<2> > b(0, af::c_grid<2>(0,0));
    scitbx::matrix::inversion_in_place(
        dat_in.ref().begin(),
        static_cast<std::size_t>(dat_in.accessor()[0]),
        b.begin(),
        static_cast<std::size_t>(b.accessor()[0]));

    //prod = af::matrix_multiply(matr01, matr02);

    return prod;
  }

*/

  int write_2d(flex_double & data2d) {
    int ncol=data2d.accessor().all()[1];
    int nrow=data2d.accessor().all()[0];
    int num=0;

    for (int row = 0; row<=nrow-1;row++) {
      if (row==0){
        std::cout << "\n  [ [ ";
      } else {
        std::cout << "\n    [ ";
      }

      for (int col = 0; col<=ncol-1;col++) {
        //printf(" %3d ", int(data2d(row,col)));

        //[ 6 ] = minimum width (no maximum given)
        //[ 2 ] = precision after the period
        printf("%6.2f ", data2d(row,col));

        num++;

        //std::cout << int(matx2d[row][col]) << " ,   ";
      }
      //fflush(stdout);
      std::cout << "  ]";
    }
    std::cout << " ] \n";
  return num;
  }

  flex_double add_2d(flex_double descriptor, flex_double data2d, flex_double total) {
    flex_double data2dreturn(total);
    int ncol_in = data2d.accessor().all()[1];
    int nrow_in = data2d.accessor().all()[0];
    int ncol_tot = total.accessor().all()[1];
    int nrow_tot = total.accessor().all()[0];
    double centr_col = descriptor(0,0);
    double centr_row = descriptor(0,1);
    double scale = descriptor(0,2);
    double tot_row, tot_col;
    double x_contrib, y_contrib, x_pix_pos, y_pix_pos;
    int xpos_ex, ypos_ex;
    /*
    std::cout <<"\n ncol_tot ="<< ncol_tot <<"\n";
    std::cout <<"\n nrow_tot ="<< nrow_tot <<"\n";
    std::cout <<"\n centr_col ="<< centr_col <<"\n";
    std::cout <<"\n centr_row ="<< centr_row <<"\n";
     */

    int tot_row_centr=int(nrow_tot / 2);
    int tot_col_centr=int(ncol_tot / 2);

    for (int row = 0; row <= nrow_in - 1;row++) {
      for (int col = 0; col <= ncol_in - 1;col++) {
        tot_row = row + tot_row_centr - centr_row + 1;
        tot_col = col + tot_col_centr - centr_col + 1;
        if (tot_row >= 0 and tot_col >= 0 and tot_row < nrow_tot and tot_col < ncol_tot) {

          x_pix_pos=tot_col - int(tot_col);
          if (x_pix_pos < 0.5) {
            x_contrib = 0.5 + x_pix_pos;
            xpos_ex = tot_col - 1;
          } else if ( x_pix_pos > 0.5 ) {
            x_contrib = 1.5 - x_pix_pos;
            xpos_ex = tot_col + 1;
          } else {
            x_contrib = 1;
            xpos_ex = tot_col;
          }

          y_pix_pos=tot_row - int(tot_row);
          if (y_pix_pos < 0.5) {
            y_contrib = 0.5 + y_pix_pos;
            ypos_ex = tot_row - 1;
          } else if ( y_pix_pos > 0.5 ) {
            y_contrib = 1.5 - y_pix_pos;
            ypos_ex = tot_row + 1;
          } else {
            y_contrib = 1;
            ypos_ex = tot_row;
          }
          /*
          std::cout << "\n tot_col =" << tot_col << "\n";
          std::cout << "\n double(int(tot_row)) =" << double(int(tot_col)) << "\n";

          std::cout << "\n y_contrib =" << y_contrib;
          std::cout << "\n x_contrib =" << x_contrib << "\n____________________________________________________________";
           */
          data2dreturn(tot_row, tot_col)=total(tot_row, tot_col) + data2d(row, col) * scale * x_contrib * y_contrib;
          if( xpos_ex != tot_col or ypos_ex != tot_row ){
            if( xpos_ex != tot_col ){
              data2dreturn(tot_row, xpos_ex)=total(tot_row, xpos_ex) + data2d(row,col) * scale * (1 - x_contrib) * y_contrib;
            }
            if( ypos_ex != tot_row ){
              data2dreturn(ypos_ex, tot_col)=total(ypos_ex, tot_col) + data2d(row, col) * scale * x_contrib * (1 - y_contrib);
            }
            if( xpos_ex != tot_col and ypos_ex != tot_row ){
              data2dreturn(ypos_ex, xpos_ex)=total(ypos_ex, xpos_ex) + data2d(row, col)* scale * (1 - x_contrib) * (1 - y_contrib);
            }
          }


        } else {
          std::cout << "\n ERROR Not fitting in the area to be added \n";
          std::cout <<"  tot_row =" <<tot_row << "  tot_col =" << tot_col <<
              "  nrow_tot =" << nrow_tot << "  ncol_tot =" << ncol_tot << "\n";
          std::cout <<" centr_col =" <<  centr_col << " centr_row =" << centr_row << "\n";
        }
      }
    }
    return data2dreturn;
  }


  vec2<double> fitting_2d(flex_double descriptor, flex_double data2d, flex_double backg2d, flex_double profile2d) {
      vec2<double> integr_data(0,1);
      int ncol = profile2d.accessor().all()[1];
      int nrow = profile2d.accessor().all()[0];
      int counter = 0;
      flex_double data2dmov(profile2d.accessor(), 0);
      flex_double backg2dmov(profile2d.accessor(), 0);
      descriptor(0,2) = 1;

      data2dmov = add_2d(descriptor, data2d, data2dmov);

      backg2dmov = add_2d(descriptor, backg2d, backg2dmov);

      // Counting how many pixels are useful so far
      for (int row = 0; row <= nrow - 1; row++) {
        for (int col = 0; col <= ncol - 1; col++) {
          if (data2dmov(row,col) != backg2dmov(row,col) ) {
            counter++ ;
          }
        }
      }

      // Building a set 1D lists with the useful pixels
      double iexpr_lst[counter];
      double imodl_lst[counter];
      double iback_lst[counter];
      double modl_scal_lst[counter];
      double scale = 0, sum = 0, i_var;
      //double bkg_var;
      double avg_i_scale, diff, df_sqr;
      counter = 0;
      for (int row = 0; row <= nrow - 1; row++) {
        for (int col = 0; col <= ncol - 1; col++) {
          if (data2dmov(row,col) != backg2dmov(row,col) and profile2d(row,col) > 0 ) {
            iexpr_lst[counter] = data2dmov(row,col) - backg2dmov(row,col);
            imodl_lst[counter] = profile2d(row,col);// * conv_scale;
            iback_lst[counter] = backg2dmov(row,col);
            counter++ ;
          }
        }
      }

      // finding the scale needed to fit profile list to experiment list
      for (int i = 0; i < counter; i++){
        scale = iexpr_lst[i] / imodl_lst[i];
        sum += scale * imodl_lst[i];
      }

      avg_i_scale = sum;
      for (int i = 0; i < counter; i++){
        modl_scal_lst[i] =imodl_lst[i] * avg_i_scale;
      }

      sum = 0;
      for (int i = 0; i < counter; i++){
        diff = (modl_scal_lst[i] - iexpr_lst[i]);
        df_sqr = diff * diff;
        sum += df_sqr;
      }

      i_var = sum;

      sum = 0;
        for (int i = 0; i < counter; i++){
          df_sqr = iback_lst[i] * iback_lst[i];
          sum += df_sqr;
        }
      //bkg_var = sum;

      integr_data[0] = avg_i_scale;//                            // intensity
      integr_data[1] = i_var;// + bkg_var ;          // intensity variance

      return integr_data;;
    }

    flex_double subtrac_bkg_2d(flex_double data2d, flex_double backg2d) {
      int ncol = data2d.accessor().all()[1];
      int nrow = data2d.accessor().all()[0];
      flex_double data2dreturn(data2d.accessor(), 0);
      for (int row = 0; row <= nrow - 1; row++) {
        for (int col = 0; col <= ncol - 1; col++) {
          if (data2d(row,col) > backg2d(row,col) ) {
            data2dreturn(row,col) = data2d(row,col) - backg2d(row,col);
          } else {
            data2dreturn(row,col) = 0.0;
          }
        }
      }
      return data2dreturn;
    }

  void rotate (float& x, float& y, float delta_ang){
    float ang,dist;
    ang = atan2(x, y);
    dist = sqrt(x * x + y * y);
    ang = ang + delta_ang * pi;
    x = dist * sin(ang);
    y = dist * cos(ang);
  }


  flex_double model_2d(int nrow, int ncol, float a, float b,
                float delta_ang, float imax, float asp) {
    float dx,dy,dd,xc,yc;
    float mw = 0.5;
    float cntnt, gss, lrz, i_tt;
    // float tot=0;
    xc=float(ncol)/2;
    yc=float(nrow)/2;
    // std::cout << "\n" << "xc =" << xc << ", yc =" << yc << "\n";
    // flex_int curv3d(data2d.accessor(),0);
    flex_double curv3d(flex_grid<>(nrow, ncol),0);
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
        //curv3d(row, col) = int(i_tt*imax);
        curv3d(row, col) = i_tt*imax;

        curv3d(row, col) += rand() % 10 - 5;

        // tot+=i_tt*imax;
      }
    }
    // std::cout <<"\n"<<"tot ="<<tot<<"\n";
    return curv3d;
  }


  void get_polar (float& ang, float& dst, float x, float y){
    ang = atan2(x, y);
    if(ang<0){
      ang+=pi*2.0;
    }
    dst = sqrt(x * x + y * y);
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
    //int   xcol_tmp[int(contr)];
    //int   yrow_tmp[int(contr)];
    int area=contr;
    contr = 0;
    for (int row = 0; row < nrow; row++) {
      for (int col = 0; col < ncol; col++) {
        if( mask2d(row, col)==1 ){
          contr++;
          dx = float(col) - xpos;
          dy = float(row) - ypos;
          get_polar (angl_tmp[contr], dist_tmp[contr], dx, dy);
          //xcol_tmp[contr]=col;
          //yrow_tmp[contr]=row;
        }
      }
    }

//    for (int i=0;i<contr;i++){
//      std::cout << "\n" << "ang ="<<angl_tmp[i];
//    }


    pie_size=int(2.0 * pi_sqrt * sqrt(float(area)));
    //std::cout << "pie_size =" << pie_size << "\n";
    float dist[int(pie_size)];
    float angl[int(pie_size)];
    //int   xcol[int(pie_size)];
    //int   yrow[int(pie_size)];

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
            //xcol[pos] = xcol_tmp[locl_pos];
            //yrow[pos] = yrow_tmp[locl_pos];
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
           std::cout <<"___________________________________________________________ exeption";
        }
        else if (angl[next_pos]!=-1 and dist[next_pos]!=-1
            and angl[prev_pos]!=-1 and dist[prev_pos]!=-1){
//                                                                                 //  starts block to be fixed
//          avg_col=float(xcol_tmp[prev_pos]+xcol_tmp[next_pos])/2.0;
//          avg_row=float(yrow_tmp[prev_pos]+yrow_tmp[next_pos])/2.0;
//
//          dx = float(avg_col) - xpos;
//          dy = float(avg_row) - ypos;
//          get_polar (angl[pos], dist[pos], dx, dy);
//                                                                                 //  ends block to be fixed
          angl[pos]=(angl[prev_pos]+angl[next_pos])/2.0;
          dist[pos]=(dist[prev_pos]+dist[next_pos])/2.0;

        }
      }
    }

//
//    for (int pos=0; pos < pie_size; pos++){
//      std::cout <<"________________________________________\n"<<"pos ="<<pos<<"\n";
//      std::cout <<"ang ="<<angl[pos]<< "\n"<<"dst ="<<dist[pos]<<"\n";
//    }
//
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
        std::cout <<"____"<< pos1 <<"____here";
      }
    }
    std::cout << "__fn=__" << fin_pos << "\n";
    if (fin_pos == -1) {
      std::cout <<"____________________________________________________ -1 ang\n";
    }

    float fin_ang = angl[fin_pos];

    return fin_ang;
  }


  flex_double test_compare_2d(flex_double descriptor, flex_double data2d, flex_double total) {
    flex_double data2dreturn(total);
    int ncol_in = data2d.accessor().all()[1];
    int nrow_in = data2d.accessor().all()[0];
    int ncol_tot = total.accessor().all()[1];
    int nrow_tot = total.accessor().all()[0];


  /*
  af::versa< double, af::c_grid<2> > test_compare_2d(
    const af::const_ref< double, af::c_grid<2> > &descriptor,
    const af::const_ref< double, af::c_grid<2> > &data2d,
    const af::const_ref< double, af::c_grid<2> > &total) {
    af::versa< double, af::c_grid<2> > data2dreturn(total.accessor(),0);
    int ncol_in=data2d.accessor()[1];
    int nrow_in=data2d.accessor()[0];
    int ncol_tot=total.accessor()[1];
    int nrow_tot=total.accessor()[0];
  */

    for (int row = 0; row<nrow_in;row++) {
      if (row==0){
        std::cout << "\n  [ [ ";
      } else {
        std::cout << "\n    [ ";
      }
      for (int col = 0; col<ncol_in;col++) {
        //[ 6 ] = minimum width (no maximum given)
        //[ 2 ] = precision after the period
        printf("%6.3f ", data2d(row,col));
        //std::cout << int(matx2d[row][col]) << " ,   ";
      }
      //fflush(stdout);
      std::cout << "  ]";
    }
    std::cout << " ] \n";

    for (int row = 0; row<nrow_tot;row++) {
      if (row==0){
        std::cout << "\n  [ [ ";
      } else {
        std::cout << "\n    [ ";
      }

      for (int col = 0; col<ncol_tot;col++) {

        //[ 6 ] = minimum width (no maximum given)
        //[ 2 ] = precision after the period
        printf("%6.3f ", total(row,col));

        //std::cout << int(matx2d[row][col]) << " ,   ";
      }
      //fflush(stdout);
      std::cout << "  ]";
    }
    std::cout << " ] \n";

    return data2dreturn;

  }




}}


#endif
