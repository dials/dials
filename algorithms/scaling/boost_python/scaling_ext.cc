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

  BOOST_PYTHON_MODULE(dials_scaling_ext)
  {
    export_elementwise_square();
    export_sph_harm_table();
    export_rotate_vectors_about_axis();
    export_calc_theta_phi();
    export_calc_sigmasq();
    export_row_multiply();
    export_determine_outlier_indices();
  }

}} // namespace dials_scaling::boost_python
