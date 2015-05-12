/*
 * rgb_2d.h
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: Luis Fuentes-Montero (Luiso)
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

//#ifndef DIALS_RGB_IMG_BUILDER_H
//#define DIALS_RGB_IMG_BUILDER_H
#include <iostream>
#include <string>
#include <scitbx/array_family/flex_types.h>
#include <cmath>

//namespace dials { namespace viewer { namespace boost_python {
  using scitbx::af::flex_double;
  using scitbx::af::flex_int;
  using scitbx::af::flex_grid;


  int get_digits( double nm, int (&dgt_num)[12]){

    int err_cod = 0;

    char asc_str[12];

    std::string str;
    std::cout << "nm =" << nm << "\n";
    sprintf (asc_str, "           ");

    sprintf (asc_str, "%g", nm);
    std::cout << "asc_str = <<" << asc_str << ">>\n\n";
    str = asc_str;

    for (int i = 0; i < 12; i++){
      str = asc_str[i];
      dgt_num[i] = asc_str[i] - 48;

      if( dgt_num[i] < 0 or dgt_num[i] > 9 ){
        if( str == "." ){
          dgt_num[i] = 10;
        }else if( str == "e" ){
          dgt_num[i] = 11;
        }else if( str == "+" ){
          dgt_num[i] = 12;
        }else if( str == "-" ){
          dgt_num[i] = 13;
        }else if( str == " " ){
          dgt_num[i] = 14;
        }else if( asc_str[i] == 0 ){
          dgt_num[i] = 15;
        }else{
          std::cout << "\nfound \'" << str << "\' and not converted";
          err_cod = asc_str[i];
          std::cout << "\n char(" << err_cod << ") found\n";
          err_cod = 1;
        }
      }

    std::cout << " {" << asc_str[i] << "} " << "["<< dgt_num[i] <<"],";
    }

    std::cout << "\nerr_cod =" << err_cod << "\n";

    return err_cod;
  }

//}}}

//#endif
