#include <boost/python.hpp>
#include <boost/python/def.hpp>

using namespace boost::python;

namespace dials_scaling { namespace boost_python {

  void export_elementwise_square();
  void export_sph_harm_table();
  void export_rotate_vectors_about_axis();
  void export_calc_theta_phi();
  void export_calc_sigmasq();
  void export_row_multiply();
  void export_determine_outlier_indices();
  void export_calc_dIh_by_dpi();
  void export_calc_jacobian();
  void export_calculate_harmonic_tables_from_selections();
  void export_calc_lookup_index();
  void export_create_sph_harm_lookup_table();
  void export_gaussian_smoother_first_fixed();
  void export_limit_outlier_weights();
  void export_split_unmerged();
  void export_mean_sample_variance();
  void export_average_intensity_variance();
  void export_mean_sample_variance_unweighted();
  void export_average_intensity_variance_unweighted();
  void export_split_into_hemispheres();

  BOOST_PYTHON_MODULE(dials_scaling_ext) {
    export_elementwise_square();
    export_sph_harm_table();
    export_rotate_vectors_about_axis();
    export_calc_theta_phi();
    export_calc_sigmasq();
    export_row_multiply();
    export_determine_outlier_indices();
    export_calc_dIh_by_dpi();
    export_calc_jacobian();
    export_calculate_harmonic_tables_from_selections();
    export_calc_lookup_index();
    export_create_sph_harm_lookup_table();
    export_gaussian_smoother_first_fixed();
    export_limit_outlier_weights();
    export_split_unmerged();
    export_mean_sample_variance();
    export_average_intensity_variance();
    export_mean_sample_variance_unweighted();
    export_average_intensity_variance_unweighted();
    export_split_into_hemispheres();
  }

}}  // namespace dials_scaling::boost_python
