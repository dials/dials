#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/scratch/luiso_s/lui_example_helper.h>

namespace dials { namespace scratch { namespace boost_python {

  using namespace boost::python;
  // testing
  void luiso_s_scratch_ext() {
    def("hello_tst", &hello_tst);
    def("write_2d", &write_2d, (arg("data2d")));
    def("add_2d", &add_2d, arg("descriptor"), arg("data2d"), arg("total"));
    def("subtrac_bkg_2d", &subtrac_bkg_2d, (arg("data2d"), arg("background2d")));
    def("fitting_2d", &fitting_2d, (arg("descriptor"), arg("data2d"), arg("background2d"), arg("profile2d") ));

    def("model_2d", &model_2d, (arg("nrow")=100, arg("ncol")=100,
           arg("a") = 10, arg("b") = 20, arg("delta_ang") = 1,
           arg("imax") = 50 ,   arg("asp") = 0.5 ) );

    def("measure_2d_angl", &measure_2d_angl, ( arg("data2d"),
         arg("mask2d"), arg("xpos")=1, arg("ypos")=20 ) );
  }

}}}
