/*
 * fft3d.cc
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
#include <dials/algorithms/indexing/fft3d.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_fft3d() {
    def("sampling_volume_map",
        &sampling_volume_map,
        (arg("data"),
         arg("angle_ranges"),
         arg("s0"),
         arg("m2"),
         arg("rl_grid_spacing"),
         arg("d_min"),
         arg("b_iso")));

    def("clean_3d",
        &clean_3d,
        (arg("dirty_beam"), arg("dirty_map"), arg("n_peaks"), arg("gamma") = 1));

    def("map_centroids_to_reciprocal_space_grid",
        &map_centroids_to_reciprocal_space_grid,
        (arg("grid"),
         arg("reciprocal_space_vectors"),
         arg("selection"),
         arg("d_min"),
         arg("b_iso") = 0));
  }

}}}  // namespace dials::algorithms::boost_python
