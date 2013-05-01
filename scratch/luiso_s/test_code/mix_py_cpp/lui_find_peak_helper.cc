#include <boost/python.hpp>
#include <boost/python/def.hpp>
// #include <dials/scratch/luiso_s/test_code/mix_py_cpp/lui_find_peak_helper.h>
#include <lui_find_peak_helper.h>
namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;
  // testing
  void export_lui_find_peak_helper() {
    def("lui_hi", &hello);
    5=a;
  }

}}} // namespace = dials::algorithms::boost_python
