/*
 * clipping_ext.cc
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
#include <boost_adaptbx/std_pair_conversion.h>
#include <dials/algorithms/polygon/clip/clip.h>

namespace dials { namespace algorithms { namespace polygon { namespace clip {
  namespace boost_python {

    using namespace boost::python;
    using boost_adaptbx::std_pair_conversions::to_tuple;

    void export_clip() {
      to_tuple<vert2, bool>();

      def("simple_with_convex", &simple_with_convex, (arg("subject"), arg("target")));
      def("simple_with_rect", &simple_with_rect, (arg("poly"), arg("rect")));
      def("triangle_with_triangle",
          &triangle_with_triangle,
          (arg("subject"), arg("target")));
      def("triangle_with_convex_quad",
          &triangle_with_convex_quad,
          (arg("subject"), arg("target")));
      def("quad_with_triangle", &quad_with_triangle, (arg("subject"), arg("target")));
      def("quad_with_convex_quad",
          &quad_with_convex_quad,
          (arg("subject"), arg("target")));
      def("line_with_rect", &line_with_rect, (arg("line"), arg("rect")));
    }

    BOOST_PYTHON_MODULE(dials_algorithms_polygon_clip_ext) {
      export_clip();
    }

}}}}}  // namespace dials::algorithms::polygon::clip::boost_python
