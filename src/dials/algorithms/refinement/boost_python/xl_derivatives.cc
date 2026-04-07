// xl_derivatives.cc
// boost.python wrapper for compute_xl_derivative_one.
// Registered into dials_refinement_helpers_ext.

#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/args.hpp>
#include "../parameterisation/xl_derivatives.h"

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  // Returns Python tuple (dpv_flex_vec3, dphi_flex_double).
  // Flex arrays convert automatically to af::const_ref<T> via the standard
  // cctbx boost.python converters (already used in parameterisation_helpers.cc).
  // The returned af::shared<vec3<double>> and af::shared<double> convert back
  // to flex.vec3_double / flex.double automatically.
  boost::python::tuple compute_xl_derivative_one_py(
    scitbx::mat3<double> const& der,
    scitbx::af::const_ref<scitbx::mat3<double> > const& fixed_rotation,
    scitbx::af::const_ref<scitbx::mat3<double> > const& setting_rotation,
    scitbx::af::const_ref<scitbx::mat3<double> > const& B_or_U,
    scitbx::af::const_ref<scitbx::vec3<double> > const& h,
    scitbx::af::const_ref<scitbx::vec3<double> > const& axis,
    scitbx::af::const_ref<double> const& phi_calc,
    scitbx::af::const_ref<scitbx::vec3<double> > const& s1,
    scitbx::af::const_ref<scitbx::vec3<double> > const& e_X_r,
    scitbx::af::const_ref<double> const& e_r_s0,
    scitbx::af::const_ref<scitbx::mat3<double> > const& D,
    bool b_matrix) {
    XlDerivResult r = compute_xl_derivative_one(der,
                                                fixed_rotation,
                                                setting_rotation,
                                                B_or_U,
                                                h,
                                                axis,
                                                phi_calc,
                                                s1,
                                                e_X_r,
                                                e_r_s0,
                                                D,
                                                b_matrix);
    return boost::python::make_tuple(r.dpv, r.dphi);
  }

  void export_xl_derivatives() {
    def("compute_xl_derivative_one",
        &compute_xl_derivative_one_py,
        (arg("der"),
         arg("fixed_rotation"),
         arg("setting_rotation"),
         arg("B_or_U"),
         arg("h"),
         arg("axis"),
         arg("phi_calc"),
         arg("s1"),
         arg("e_X_r"),
         arg("e_r_s0"),
         arg("D"),
         arg("b_matrix")));
  }

}}}  // namespace dials::refinement::boost_python
