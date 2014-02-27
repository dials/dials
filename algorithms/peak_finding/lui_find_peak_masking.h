#ifndef DIALS_ALGORITHMS_PEAK_FINDING_LUI_HELPER_H
#define DIALS_ALGORITHMS_PEAK_FINDING_LUI_HELPER_H
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <iostream>
#include <cmath>
namespace dials { namespace algorithms {

  af::versa< double, af::c_grid<2> > find_mask_2d(
        const af::const_ref< double, af::c_grid<2> > &data2d
      , const af::const_ref< double, af::c_grid<2> > &data2dsmoth
      , int tot_times) {
        af::versa< double, af::c_grid<2> > mask2d(data2d.accessor(),0);
        std::size_t ncol=data2d.accessor()[1];
        std::size_t nrow=data2d.accessor()[0];
        for (int row = 1; row<nrow - 1;row++) {
          for (int col = 1; col<ncol - 1;col++){
            if (data2d(row,col) > data2dsmoth(row,col) + 15.0){
              mask2d(row,col) = 1;
            }else{
              mask2d(row,col) = 0;
            }
          }
        }

        int nblock = 10;
        int col_block_size = int(double(ncol)/double(nblock));
        int row_block_size = int(double(nrow)/double(nblock));

        int col_ini;
        int col_end;
        int col_block_num;

        int row_ini;
        int row_end;
        int row_block_num;

        double dif, tot, c, var;

        for(
        col_ini = 0,
        col_end = col_block_size,
        col_block_num = 0
        ;
        col_block_num < nblock
        ;
        col_block_num++,
        col_ini = col_block_num * col_block_size,
        col_end = (col_block_num + 1) * col_block_size
        ){

          for(
          row_ini = 0,
          row_end = row_block_size,
          row_block_num = 0
          ;
          row_block_num < nblock
          ;
          row_block_num++,
          row_ini = row_block_num * row_block_size,
          row_end = (row_block_num + 1) * row_block_size

          ){
            std::cout << "\ncol_ini =" << col_ini <<
            "\ncol_end =" << col_end << "\ncol_block_num =" << col_block_num;
            std::cout << "\nrow_ini =" << row_ini << "\nrow_end =" << row_end <<
             "\n" << "row_block_num =" << row_block_num;

            tot = 0;
            c = 0;
            for (int row = row_ini + 1; row < row_end - 1;row++) {
              for (int col = col_ini + 1; col < col_end - 1;col++) {
                dif =  data2dsmoth(row,col) - data2d(row,col);
                tot += abs(dif);
                c++;
              }
            }
            var = tot/c;
            std::cout << "\n" << "tot = "<< tot << "\n" << "var =" << var
             << "\n____________________\n";
          }


        }


    return mask2d;
  }

}}

#endif
