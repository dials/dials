#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../gaussian_smoother_3D.h"

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  void export_gaussian_smoother_3D() {
    class_<GaussianSmoother3D>("GaussianSmoother3D", no_init)
      .def(init<vec2<double>,
                std::size_t,
                vec2<double>,
                std::size_t,
                vec2<double>,
                std::size_t>((arg("x_range"),
                              arg("num_x_intervals"),
                              arg("y_range"),
                              arg("num_y_intervals"),
                              arg("z_range"),
                              arg("num_z_intervals"))))
      .def("set_smoothing", &GaussianSmoother3D::set_smoothing)
      .def("num_x_values", &GaussianSmoother3D::num_x_values)
      .def("num_y_values", &GaussianSmoother3D::num_y_values)
      .def("num_z_values", &GaussianSmoother3D::num_z_values)
      .def("num_x_samples", &GaussianSmoother3D::num_x_samples)
      .def("num_y_samples", &GaussianSmoother3D::num_y_samples)
      .def("num_z_samples", &GaussianSmoother3D::num_z_samples)
      .def("num_x_average", &GaussianSmoother3D::num_x_average)
      .def("num_y_average", &GaussianSmoother3D::num_y_average)
      .def("num_z_average", &GaussianSmoother3D::num_z_average)
      .def("sigma", &GaussianSmoother3D::sigma)
      .def("x_spacing", &GaussianSmoother3D::x_spacing)
      .def("y_spacing", &GaussianSmoother3D::y_spacing)
      .def("z_spacing", &GaussianSmoother3D::z_spacing)
      .def("x_positions", &GaussianSmoother3D::x_positions)
      .def("y_positions", &GaussianSmoother3D::y_positions)
      .def("z_positions", &GaussianSmoother3D::z_positions)
      .def("value_weight", &GaussianSmoother3D::value_weight)
      .def("multi_value_weight", &GaussianSmoother3D::multi_value_weight);
  }

}}}  // namespace dials::refinement::boost_python
