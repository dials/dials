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

    class_<w_t>("AssignIndices", no_init)
      .def(init<af::const_ref<scitbx::vec3<double> > const &,
                af::const_ref<double> const &,
                af::const_ref<scitbx::mat3<double> > const &,
                double>((arg("reciprocal_space_points"),
                         arg("phi"),
                         arg("UB_matrices"),
                         arg("tolerance") = 0.3)))
      .def("miller_indices", &w_t::miller_indices)
      .def("crystal_ids", &w_t::crystal_ids);
  }

  void export_assign_indices_local() {
    typedef AssignIndicesLocal w_t;
    class_<w_t>("AssignIndicesLocal", no_init)
      .def(init<af::const_ref<scitbx::vec3<double> > const &,
                af::const_ref<double> const &,
                af::const_ref<scitbx::mat3<double> > const &,
                const double,
                const double,
                const double,
                const int>((arg("reciprocal_space_points"),
                            arg("phi"),
                            arg("UB_matrices"),
                            arg("epsilon") = 0.05,
                            arg("delta") = 8,
                            arg("l_min") = 0.8,
                            arg("nearest_neighbours") = 20)))
      .def("miller_indices", &w_t::miller_indices)
      .def("crystal_ids", &w_t::crystal_ids);
  }

  BOOST_PYTHON_MODULE(dials_algorithms_indexing_ext) {
    export_fft3d();
    export_assign_indices();
    export_assign_indices_local();
  }

}}}  // namespace dials::algorithms::boost_python
