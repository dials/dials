#include <boost/python.hpp>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <cctype>

#include <x2tbx.h>

namespace x2tbx { 
  namespace ext {

    static merged_isig merge(observation_list)
    {
      /* fixme implement */
      merged_isig result;
      result.I = 0.0;
      result.sigI = 0.0;
      return result;
    }

    static float 
    isig(scitbx::af::const_ref<float> const & i_data,
	 scitbx::af::const_ref<float> const & sigi_data)
    {
      float result = 0.0;

      CCTBX_ASSERT(i_data.size() == sigi_data.size());

      for (size_t i = 0; i < i_data.size(); i++) {
	result += i_data[i] / sigi_data[i];
      }

      return result / i_data.size();
    }

    void init_module()
    {
      using namespace boost::python;
      def("isig", isig, (arg("i_data"), arg("sigi_data")));
    }

  }
} // namespace x2tbx::ext

BOOST_PYTHON_MODULE(x2tbx_ext)
{
  x2tbx::ext::init_module();
}
