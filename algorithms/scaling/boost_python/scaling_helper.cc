#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials_scratch/jbe/scaling_code/scaling_helper.h>

namespace dials_scratch { namespace scaling { namespace boost_python {

  using namespace boost::python;
    
  using scitbx::sparse::matrix;

    void export_elementwise_square()
  {
    def("elementwise_square", &elementwise_square, (
      arg("m")));
    }

    void export_sph_harm_table()
    {
    def("create_sph_harm_table", &create_sph_harm_table, (
      arg("s0_theta"),
      arg("s0_phi"),
      arg("s1_theta"),
      arg("s1_phi"),
      arg("lmax")));
    }

}}} // namespace = dials_scratch::scaling::boost_python
