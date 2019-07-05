#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/panel.h>
#include <math.h>
#include <vector>

namespace kapton {

using dxtbx::model::Detector;
using scitbx::vec2;
using scitbx::vec3;

/**
 * Implementing a c++ version of the get_kapton_path function. Significant speedups
 * observed
 */

scitbx::af::shared<double> get_kapton_path_cpp(
  boost::python::list kapton_faces,
  scitbx::af::const_ref<vec3<double> > s1_flex) {
  scitbx::af::shared<double> kapton_path_mm;
  int nfaces = boost::python::len(kapton_faces);
  std::vector<Detector> all_kapton_faces;
  // Store all the detectors ahead of time to avoid making expensive extract calls every
  // time
  for (std::size_t i = 0; i < nfaces; ++i) {
    all_kapton_faces.push_back(boost::python::extract<Detector>(kapton_faces[i]));
  }
  for (scitbx::af::const_ref<vec3<double> >::const_iterator it = s1_flex.begin();
       it != s1_flex.end();
       ++it) {
    double max_d2 = -999.9;
    scitbx::af::shared<vec3<double> > intersection_xy_list;
    // Find out intersection points of s1 vector with all faces of kapton volume
    for (std::size_t j = 0; j < nfaces; ++j) {
      Detector detector = all_kapton_faces.at(j);
      try {
        vec2<double> px = detector[0].get_ray_intersection_px(*it);
        if (px[0] > 0.0 && px[1] > 0.0 && px[0] < detector[0].get_image_size()[0]
            && px[1] < detector[0].get_image_size()[1]) {
          intersection_xy_list.push_back(
            detector[0].get_lab_coord(detector[0].get_ray_intersection(*it)));
        }
      } catch (dxtbx::error) {
        // Do nothing
      }
    }  // kapton faces
    /*
     * Now find out the maximum path length through the kapton. That is the path
     * traversed by s1 vector if no intersection or intersects with only one face then
     * path length is set to 0.0
     */
    if (intersection_xy_list.size() == 0) {
      kapton_path_mm.push_back(0.0);
    } else if (intersection_xy_list.size() == 1) {
      kapton_path_mm.push_back(0.0);
    } else {
      double d2 = 0.0;
      for (std::size_t k1 = 0; k1 < intersection_xy_list.size() - 1; ++k1) {
        for (std::size_t k2 = k1 + 1; k2 < intersection_xy_list.size(); ++k2) {
          vec3<double> pt1 = intersection_xy_list[k1];
          vec3<double> pt2 = intersection_xy_list[k2];
          d2 = (pt1[0] - pt2[0]) * (pt1[0] - pt2[0])
               + (pt1[1] - pt2[1]) * (pt1[1] - pt2[1])
               + (pt1[2] - pt2[2]) * (pt1[2] - pt2[2]);
          if (d2 > max_d2) {
            max_d2 = d2;
          }
        }
      }
      double d = sqrt(max_d2);
      kapton_path_mm.push_back(d);
    }
  }  // s1
  return kapton_path_mm;
}
}  // namespace kapton

using namespace boost::python;
namespace kapton { namespace boost_python { namespace {

  void kapton_init_module() {
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

    def("get_kapton_path_cpp", &kapton::get_kapton_path_cpp);
  }

}}}  // namespace kapton::boost_python::

BOOST_PYTHON_MODULE(dials_algorithms_integration_kapton_ext) {
  kapton::boost_python::kapton_init_module();
}
