#ifndef DIALS_ALGORITHMS_LUI_2D_BACKGROUND_H
#define DIALS_ALGORITHMS_LUI_2D_BACKGROUND_H
#include <iostream>
#include <scitbx/array_family/flex_types.h>

namespace dials { namespace algorithms {
  using scitbx::af::flex_int;

  flex_int background_subtract_2d(flex_int & data2d) {
        // std::cout << "length =" << data2d.size() << " \n";
        //std::size_t ncol=data2d.accessor().all()[1];
        //std::size_t nrow=data2d.accessor().all()[0];

        std::cout << "Hi";
        flex_int data2dsmoth(data2d.accessor(),0);

    return data2dsmoth;
  }


}}

#endif
