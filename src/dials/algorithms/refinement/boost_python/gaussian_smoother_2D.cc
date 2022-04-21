#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../gaussian_smoother_2D.h"

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  void export_gaussian_smoother_2D() {
    class_<GaussianSmoother2D>("GaussianSmoother2D", no_init)
      .def(init<vec2<double>, std::size_t, vec2<double>, std::size_t>(
        (arg("x_range"),
         arg("num_x_intervals"),
         arg("y_range"),
         arg("num_y_intervals"))))
      .def("set_smoothing", &GaussianSmoother2D::set_smoothing)
      .def("num_x_values", &GaussianSmoother2D::num_x_values)
      .def("num_y_values", &GaussianSmoother2D::num_y_values)
      .def("num_x_samples", &GaussianSmoother2D::num_x_samples)
      .def("num_y_samples", &GaussianSmoother2D::num_y_samples)
      .def("num_x_average", &GaussianSmoother2D::num_x_average)
      .def("num_y_average", &GaussianSmoother2D::num_y_average)
      .def("sigma", &GaussianSmoother2D::sigma)
      .def("x_spacing", &GaussianSmoother2D::x_spacing)
      .def("y_spacing", &GaussianSmoother2D::y_spacing)
      .def("x_positions", &GaussianSmoother2D::x_positions)
      .def("y_positions", &GaussianSmoother2D::y_positions)
      .def("value_weight", &GaussianSmoother2D::value_weight)
      .def("multi_value_weight", &GaussianSmoother2D::multi_value_weight);
  }

}}}  // namespace dials::refinement::boost_python
