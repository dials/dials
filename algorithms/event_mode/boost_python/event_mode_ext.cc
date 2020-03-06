#include <boost/python.hpp>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/panel.h>
#include <dxtbx/model/beam.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <cctype>
#include <boost/random.hpp>

namespace dials { 
  namespace algorithms { 
    namespace event_mode {

      // Create event list

      static boost::python::tuple event_list(int image_n,
					     scitbx::af::shared<size_t> idx, scitbx::af::shared<int> counts)
      {
	scitbx::af::shared<int> position, time;
	boost::random::mt19937 rng;
	boost::random::uniform_real_distribution<double> dist(0, 1);
	for (size_t j = 0; j < idx.size(); j++) {  
	  for (size_t c = 0; c < counts[j]; c++) {
	    position.push_back(idx[j]);
            time.push_back(1000*(image_n + dist(rng)));
	  }
	}
	return boost::python::make_tuple(position, time);
      }

      // Get image from event list

      static scitbx::af::shared<double> image_coord(int n, int expT,
						    scitbx::af::shared<int> pos, scitbx::af::shared<int> t)
      {
	scitbx::af::shared<double> coord;
	int t_i = expT*n;
	int t_f = expT*(n+1);
	//scitbx::af::shared<bool> deltaT;
	//scitbx::af::shared<size_t> idx;
	for (size_t j = 0; j < t.size(); j++) {
          if ((t[j] >= t_i) && (t[j] < t_f)) {
	    coord.push_back(pos[j]);
          }
	}
                
	return coord;
      }

      void init_module()
      {
	using namespace boost::python;
	def("event_list", event_list, (arg("image_n"), arg("idx"), arg("counts")));
	def("image_coord", image_coord, (arg("n"), arg("expT"), arg("pos"), arg("t")));
      }
    }
  }
} // namespace dials::algorithms::event_mode

BOOST_PYTHON_MODULE(dials_algorithms_event_mode_ext)
{
  dials::algorithms::event_mode::init_module();
}
