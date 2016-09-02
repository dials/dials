/*
 * FIXME add a header
 */

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/util/scale_down_array.h>
#include <dials/util/masking.h>

namespace dials { namespace util { namespace boost_python {

  af::shared<bool> is_inside_polygon_a(
      const af::const_ref< vec2<double> > &poly,
      const af::const_ref< vec2<double> > &xy) {
    af::shared<bool> inside(xy.size(), false);
    for (std::size_t i = 0; i < xy.size(); i++) {
      inside[i] = is_inside_polygon(poly, xy[i][0], xy[i][1]);
    }
    return inside;
  }

  using namespace boost::python;
  BOOST_PYTHON_MODULE(dials_util_ext) {

    def("scale_down_array", &scale_down_array,
        (arg("image"), arg("scale_factor")));

    def("mask_untrusted_rectangle",
        &mask_untrusted_rectangle);

    def("mask_untrusted_circle",
        &mask_untrusted_circle);

    def("mask_untrusted_resolution_range",
        &mask_untrusted_resolution_range);

    def("mask_untrusted_polygon",
        &mask_untrusted_polygon);

    def("is_inside_polygon",
        &is_inside_polygon);

    def("is_inside_polygon",
        &is_inside_polygon_a);

    class_<ResolutionMaskGenerator>("ResolutionMaskGenerator", no_init)
      .def(init<const Beam&,const Panel&>())
      .def("apply", &ResolutionMaskGenerator::apply)
      ;
  }
}}}
