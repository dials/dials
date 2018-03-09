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

}}} // namespace = dials_scratch::scaling::boost_python
