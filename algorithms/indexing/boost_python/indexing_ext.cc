/*
 * indexing_ext.cc
 *
 *  Copyright (C) 2014 Diamond Light Source
 *
 *  Author: Richard Gildea
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/indexing/index.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_fft3d();

  void export_assign_indices() {

    typedef AssignIndices w_t;

    class_<w_t>(
        "AssignIndices", no_init)
      .def(init<af::const_ref<scitbx::vec3<double> > const &,
                af::const_ref<scitbx::mat3<double> > const &,
                double>((
        arg("reciprocal_space_points"),
        arg("UB_matrices"),
        arg("tolerance") = 10)))
      .def("miller_indices", &w_t::miller_indices)
      .def("crystal_ids", &w_t::crystal_ids)
      .def("n_rejects", &w_t::n_rejects);

  }

  BOOST_PYTHON_MODULE(dials_algorithms_indexing_ext)
  {
    export_fft3d();
    export_assign_indices();
  }

}}} // namespace = dials::algorithms::boost_python
