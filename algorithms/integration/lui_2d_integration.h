/*
 * lui_2d_integration.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: Luis Fuentes-Montero
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_LUI_INTEGRATION_2D_H
#define DIALS_ALGORITHMS_LUI_INTEGRATION_2D_H

#include <iostream>
#include <scitbx/vec2.h>
#include <cmath>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/shoebox.h>

namespace dials { namespace algorithms {

  using scitbx::vec2;
  // using dials::model::Valid;
  using dials::model::Foreground;
  using dials::model::Background;
  vec2<double> raw_2d_cut(
      const af::const_ref< double, af::c_grid<2> > &data2d,
      const af::const_ref< int, af::c_grid<2> > &mask2d,
      const af::const_ref< double, af::c_grid<2> > &background2d) {
        double i_s = 0, i_bg = 0, rho_j = 0;
        double n = 0, m = 0;
        int cont = 0;
        double var_i;
        std::size_t ncol=data2d.accessor()[1];
        std::size_t nrow=data2d.accessor()[0];
        vec2<double> integr_data(0,1);

        /*
      for(int row = 0; row < ny; row++){
        for( int col = 0; col < nx; col++ ){
        */
        for (int row = 0; row<nrow;row++) {
          for (int col = 0; col<ncol;col++) {
            if ( mask2d(row,col) & Foreground ) {
              i_s  += data2d(row,col) - background2d(row,col);
              i_bg += background2d(row,col);
              m++;
            } else if(mask2d(row,col) & Background) {
              rho_j += background2d(row,col);
              n++;
            }
            cont++;
          }
        }
        if( i_bg>0 && cont>0 && m>0 && n>0 ){
          var_i = i_s + i_bg + (m / n) * ( m / n) * rho_j;
          //std::cout << "\n m =" << m << "\n n =" << n;
          //std::cout << "\n i_s =" << i_s << "\n i_bg =" << i_bg << "\n (m/n)**2 =" << double((m / n) * ( m / n)) << "\n rho_j =" << rho_j;
        } else {
          var_i = i_s;
        }

        integr_data[0]=i_s;            // intensity sumation
        integr_data[1]=var_i;          // intensity variance
        return integr_data;

    }

} }

#endif
