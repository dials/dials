#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/scaling/scaling_helper.h>

namespace dials_scaling { namespace boost_python {

  using namespace boost::python;

  using scitbx::sparse::matrix;

  void export_determine_outlier_indices() {
    def("determine_outlier_indices",
        &determine_outlier_indices,
        (arg("h_index_matrix"), arg("z_scores"), arg("zmax")));
  }

  void export_elementwise_square() {
    def("elementwise_square", &elementwise_square, (arg("m")));
  }

  void export_calc_dIh_by_dpi() {
    def("calc_dIh_by_dpi",
        &calculate_dIh_by_dpi,
        (arg("a"), arg("sumgsq"), arg("h_index_mat"), arg("derivatives")));
  }

  void export_calc_jacobian() {
    def("calc_jacobian",
        &calc_jacobian,
        (arg("derivatives"),
         arg("h_index_mat"),
         arg("Ih"),
         arg("g"),
         arg("dIh"),
         arg("sumgsq")));
  }

  void export_sph_harm_table() {
    def("create_sph_harm_table",
        &create_sph_harm_table,
        (arg("s0_theta_phi"), arg("s1_theta_phi"), arg("lmax")));
  }

  void export_rotate_vectors_about_axis() {
    def("rotate_vectors_about_axis",
        &rotate_vectors_about_axis,
        (arg("rot_axis"), arg("vectors"), arg("angles")));
  }

  void export_calc_theta_phi() {
    def("calc_theta_phi", &calc_theta_phi, (arg("xyz")));
  }

  void export_calc_sigmasq() {
    def("calc_sigmasq", &calc_sigmasq, (arg("jacobian_transpose"), arg("var_cov")));
  }

  void export_row_multiply() {
    def("row_multiply", &row_multiply, (arg("m"), arg("v")));
  }

  void export_calc_lookup_index() {
    def("calc_lookup_index",
        &calc_lookup_index,
        (arg("thetaphi"), arg("points_per_degree")));
  }

  void export_create_sph_harm_lookup_table() {
    def("create_sph_harm_lookup_table",
        &create_sph_harm_lookup_table,
        (arg("lmax"), arg("points_per_degree")));
  }

  void export_calculate_harmonic_tables_from_selections() {
    def("calculate_harmonic_tables_from_selections",
        &calculate_harmonic_tables_from_selections,
        (arg("s0_selection"), arg("s1_selection"), arg("coefficients_list")));
  }

}}  // namespace dials_scaling::boost_python
