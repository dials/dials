/*
 * spatial_interpolation_ext.cc
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
#include <scitbx/vec2.h>
#include <dials/algorithms/polygon/spatial_interpolation.h>
#include <dials/array_family/scitbx_shared_and_versa.h>

namespace dials { namespace algorithms { namespace polygon {
  namespace spatial_interpolation { namespace boost_python {

    using namespace boost::python;

    void export_spatial_interpolation() {
      class_<Match>("Match")
        .def(init<int, int, double>())
        .def_readwrite("in", &Match::in)
        .def_readwrite("out", &Match::out)
        .def_readwrite("fraction", &Match::fraction);

      def("irregular_grid_to_grid", &irregular_grid_to_grid);
      def("grid_to_irregular_grid", &grid_to_irregular_grid);

      def("regrid_irregular_grid_to_grid", &regrid_irregular_grid_to_grid);
      def("regrid_grid_to_irregular_grid", &regrid_grid_to_irregular_grid);
    }

    BOOST_PYTHON_MODULE(dials_algorithms_polygon_spatial_interpolation_ext) {
      export_spatial_interpolation();
    }

}}}}}  // namespace dials::algorithms::polygon::spatial_interpolation::boost_python
