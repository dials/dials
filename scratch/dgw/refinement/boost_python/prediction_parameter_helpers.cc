
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include "../prediction_parameter_helpers.h"

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  void export_prediction_parameter_helpers()
  {
    scitbx::af::boost_python::flex_wrapper <mat3 <double> >::plain("flex_mat3_double");

    vec3 <double> (*detector_pv_derivative_single)(
      mat3 <double>, mat3 <double>, vec3 <double>) =
        &detector_pv_derivative;

    flex_vec3_double (*detector_pv_derivative_array)(
      mat3 <double>, const flex_mat3_double&, vec3 <double>) =
        &detector_pv_derivative;

    double (*source_phi_derivative_single)(
      vec3 <double>, vec3 <double>, double) =
        &source_phi_derivative;

    flex_double (*source_phi_derivative_array)(
      vec3 <double>, const flex_vec3_double&, double) =
        &source_phi_derivative;

    vec3 <double> (*source_pv_derivative_single)(
      mat3 <double>, vec3 <double>, double, vec3 <double>) =
        &source_pv_derivative;

    flex_vec3_double (*source_pv_derivative_array)(
      mat3 <double>, vec3 <double>, const flex_double&,
      const flex_vec3_double&) =  &source_pv_derivative;

    vec3 <double> (*crystal_orientation_r_derivative_single)(
      mat3 <double>, mat3 <double>, mat3 <double>, miller_index) =
        &crystal_orientation_r_derivative;

    flex_vec3_double (*crystal_orientation_r_derivative_array)(
      mat3 <double>, const flex_mat3_double&, mat3 <double>,
      miller_index) = &crystal_orientation_r_derivative;

    double (*crystal_orientation_phi_derivative_single)(
      vec3 <double>, vec3 <double>, double) =
        &crystal_orientation_phi_derivative;

    flex_double (*crystal_orientation_phi_derivative_array)(
      const flex_vec3_double&, vec3 <double>, double) =
        &crystal_orientation_phi_derivative;

    vec3 <double> (*crystal_orientation_pv_derivative_single)(
      mat3 <double>, vec3 <double>, vec3 <double>, double) =
        &crystal_orientation_pv_derivative;

    flex_vec3_double (*crystal_orientation_pv_derivative_array)(
      mat3 <double>, const flex_vec3_double&, vec3 <double>,
      const flex_double&) = &crystal_orientation_pv_derivative;

    vec3 <double> (*crystal_cell_r_derivative_single)(
      mat3 <double>, mat3 <double>, mat3 <double>, miller_index) =
        &crystal_cell_r_derivative;

    flex_vec3_double (*crystal_cell_r_derivative_array)(
      mat3 <double>, mat3 <double>, const flex_mat3_double&, miller_index) =
        &crystal_cell_r_derivative;

    double (*crystal_cell_phi_derivative_single)(
      vec3 <double>, vec3 <double>, double) =
        &crystal_cell_phi_derivative;

    flex_double (*crystal_cell_phi_derivative_array)(
      const flex_vec3_double&, vec3 <double>, double) =
        &crystal_cell_phi_derivative;

    vec3 <double> (*crystal_cell_pv_derivative_single)(
      mat3 <double>, vec3 <double>, vec3 <double>, double) =
        &crystal_cell_pv_derivative;

    flex_vec3_double (*crystal_cell_pv_derivative_array)(
      mat3 <double>, const flex_vec3_double&, vec3 <double>, const flex_double&) =
        &crystal_cell_pv_derivative;

    def("detector_pv_derivative",
      detector_pv_derivative_single);
    def("detector_pv_derivative",
      detector_pv_derivative_array);
    def("source_phi_derivative",
      source_phi_derivative_single);
    def("source_phi_derivative",
      source_phi_derivative_array);
    def("source_pv_derivative",
      source_pv_derivative_single);
    def("source_pv_derivative",
      source_pv_derivative_array);
    def("crystal_orientation_r_derivative",
      crystal_orientation_r_derivative_single);
    def("crystal_orientation_r_derivative",
      crystal_orientation_r_derivative_array);
    def("crystal_orientation_phi_derivative",
      crystal_orientation_phi_derivative_single);
    def("crystal_orientation_phi_derivative",
      crystal_orientation_phi_derivative_array);
    def("crystal_orientation_pv_derivative",
      crystal_orientation_pv_derivative_single);
    def("crystal_orientation_pv_derivative",
      crystal_orientation_pv_derivative_array);
    def("crystal_cell_r_derivative",
      crystal_cell_r_derivative_single);
    def("crystal_cell_r_derivative",
      crystal_cell_r_derivative_array);
    def("crystal_cell_phi_derivative",
      crystal_cell_phi_derivative_single);
    def("crystal_cell_phi_derivative",
      crystal_cell_phi_derivative_array);
    def("crystal_cell_pv_derivative",
      crystal_cell_pv_derivative_single);
    def("crystal_cell_pv_derivative",
      crystal_cell_pv_derivative_array);
  }

}}} // namespace dials::refinement::boost_python
