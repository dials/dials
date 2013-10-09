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
        double i_tot = 0, tot_bkgr = 0;
        int npix_bkgr = 0, npix_mask = 0, cont = 0;
        double bkgr, sig_i;
        std::size_t ncol=data2d.accessor()[1];
        std::size_t nrow=data2d.accessor()[0];
        vec2<double> integr_data(0,1);

        for (int row = 0; row<=nrow-1;row++) {
          for (int col = 0; col<=ncol-1;col++) {
            if ( mask2d(row,col) & Foreground ){
              i_tot = i_tot + data2d(row,col) - background2d(row,col);
              npix_mask++;
            } else if(mask2d(row,col) & Background) {
              npix_bkgr++;
            }

            cont++;
            tot_bkgr+=background2d(row,col);

          }
        }
        if( tot_bkgr>0 && cont>0 && npix_mask>0 && npix_bkgr>0 ){
          bkgr = tot_bkgr / cont;
          sig_i = std::sqrt(i_tot + (1.0 + (npix_mask) / (npix_bkgr)) * (npix_mask * bkgr));

        } else {
          bkgr = 0;
          sig_i = std::sqrt(i_tot);
        }

        integr_data[0]=i_tot;        // intensity
        integr_data[1]=sig_i;          // intensity sigma

        return integr_data;

    }

} }

#endif
