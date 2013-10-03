//++ b/trunk/scratch/luiso_s/boost_python/lui_example_helper.cc
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/shoebox/mask_builder.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;
  // testing
  void export_mask_builder() {
    def("hello_tst", &hello_tst);
  }

}}}}
