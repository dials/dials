/*
 * transform_ext.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/reflexion_basis/rebin_pixels.h>

namespace dials { namespace algorithms { namespace reflexion_basis {
  namespace boost_python {

  using namespace boost::python;
  using scitbx::af::int2;
  using scitbx::af::flex_grid;
  
  inline
  flex_double rebin_pixels_wrapper(const flex_double &input, 
      const flex_vec2_double &inputxy, int2 size) {
    flex_double output(flex_grid<>(size[0], size[1]));
    rebin_pixels(output, input, inputxy);
    return output;
  }

  void export_rebin_pixels() 
  {
    def("rebin_pixels", &rebin_pixels_wrapper, (
      arg("input"), arg("inputxy")));
  }

  BOOST_PYTHON_MODULE(dials_algorithms_reflexion_basis_transform_ext)
  {
    export_rebin_pixels();
  }

}}}} // namespace = dials::algorithms::reflexion_basis::boost_python
