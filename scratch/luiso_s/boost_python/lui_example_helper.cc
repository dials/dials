#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/scratch/luiso_s/lui_example_helper.h>

namespace dials { namespace scratch { namespace boost_python {

  using namespace boost::python;
  // testing
  void luiso_s_scratch_ext() {
    def("hello_tst", &hello_tst);
    def("tst_01", &tst_01, (arg("data2d")));
    def("raw_2d_cut", &raw_2d_cut, (arg("data2d")), (arg("mask2d")),
        (arg("background2d")));


  }

}}}
