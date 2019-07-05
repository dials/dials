/*
 * filter_ext.cc
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

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_summed_area();
  void export_mean_and_variance();
  void export_index_of_dispersion_filter();
  void export_convolve();
  void export_median();
  void export_distance();
  void export_anisotropic_diffusion();

  BOOST_PYTHON_MODULE(dials_algorithms_image_filter_ext) {
    export_summed_area();
    export_mean_and_variance();
    export_index_of_dispersion_filter();
    export_convolve();
    export_median();
    export_distance();
    export_anisotropic_diffusion();
  }

}}}  // namespace dials::algorithms::boost_python
