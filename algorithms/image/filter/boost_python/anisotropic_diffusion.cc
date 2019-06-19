/*
 * anisotropic_diffusion.cc
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
#include <dials/algorithms/image/filter/anisotropic_diffusion.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_anisotropic_diffusion() {
    def("anisotropic_diffusion",
        &anisotropic_diffusion,
        (arg("data"), arg("niter") = 1, arg("kappa") = 50, arg("gamma") = 0.1));

    def("anisotropic_diffusion",
        &masked_anisotropic_diffusion,
        (arg("data"),
         arg("mask"),
         arg("niter") = 1,
         arg("kappa") = 50,
         arg("gamma") = 0.1));
  }

}}}  // namespace dials::algorithms::boost_python
