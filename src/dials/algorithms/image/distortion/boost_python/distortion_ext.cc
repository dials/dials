#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/image/distortion/plane.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_create_plane_distortion_maps() {
    class_<PlaneLinearTransformationMaps>("PlaneLinearTransformationMaps", no_init)
      .def(init<const Panel&,
                scitbx::mat2<double>,
                scitbx::vec3<double>,
                scitbx::vec3<double>,
                scitbx::vec3<double>>(
        (arg("panel"), arg("ellipse_matrix"), arg("fast"), arg("slow"), arg("mid"))))
      .def("get_dx", &PlaneLinearTransformationMaps::get_dx)
      .def("get_dy", &PlaneLinearTransformationMaps::get_dy);
  }

  BOOST_PYTHON_MODULE(dials_algorithms_image_distortion_ext) {
    export_create_plane_distortion_maps();
  }

}}}  // namespace dials::algorithms::boost_python
