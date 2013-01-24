#include <x2tbx.h>

namespace x2tbx {
  namespace ext {

    static merged_isig merge(observation_list ol)
    {
      float sum_wi = 0.0;
      float sum_w = 0.0;

      for (size_t j = 0; j < ol.size(); j ++) {
        float i = ol[j].I;
        float w = 1.0 / (ol[j].sigI * ol[j].sigI);
        sum_w += w;
        sum_wi += w * i;
      }

      merged_isig result;
      result.I = sum_wi / sum_w;
      result.sigI = sqrt(1.0 / sum_w);
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

    static float
    isig_proper(scitbx::af::const_ref<cmil::index<int> > const & indices,
                scitbx::af::const_ref<float> const & i_data,
                scitbx::af::const_ref<float> const & sigi_data)
    {
      float result = 0.0;

      unmerged_reflections ur;
      unmerged_reflections_iterator uri;
      observation o;

      CCTBX_ASSERT(indices.size() == i_data.size());
      CCTBX_ASSERT(i_data.size() == sigi_data.size());

      for (size_t i = 0; i < i_data.size(); i++) {
        o.I = i_data[i];
        o.sigI = sigi_data[i];
        o.property = 0.0;
        o.flag = 0;

        ur[indices[i]].push_back(o);
      }

      int unique = 0;

      for (uri = ur.begin(); uri != ur.end(); ++uri) {
        merged_isig mi = merge(uri->second);
        result += mi.I / mi.sigI;
        unique += 1;
      }

      return result / unique;
    }

    void init_module()
    {
      using namespace boost::python;
      def("isig", isig, (arg("i_data"), arg("sigi_data")));
      def("isig_proper", isig_proper,
          (arg("indices"), arg("i_data"), arg("sigi_data")));
    }
  }

  scitbx::af::shared<cmil::index<int> >
  resolutionizer::sorted_indices(void)
  {
    scitbx::af::shared<cmil::index<int> > result;
    unmerged_reflections_iterator uri;
    for (uri = ur.begin(); uri != ur.end(); ++uri) {
      result.push_back(uri->first);
    }

    std::sort(result.begin(), result.end(), sorter_by_resolution(unit_cell));

    return result;
  }

  bool
  resolutionizer::compare_resolution(cmil::index<int> const & a,
                                     cmil::index<int> const & b)
  {
    return unit_cell.d(a) < unit_cell.d(b);
  }

  void
  resolutionizer::set_unit_cell(scitbx::af::tiny<double, 6>
                                unit_cell_parameters)
  {
    unit_cell = cuc::unit_cell(unit_cell_parameters);
  }

  void
  resolutionizer::setup(scitbx::af::const_ref<cmil::index<int> > const & indices,
                        scitbx::af::const_ref<float> const & i_data,
                        scitbx::af::const_ref<float> const & sigi_data)
  {
    CCTBX_ASSERT(indices.size() == i_data.size());
    CCTBX_ASSERT(i_data.size() == sigi_data.size());

    observation o;

    for (size_t i = 0; i < i_data.size(); i++) {
      o.I = i_data[i];
      o.sigI = sigi_data[i];
      o.property = 0.0;
      o.flag = 0;

      ur[indices[i]].push_back(o);
    }

    s_indices = sorted_indices();
  }

  void
  resolutionizer::setup_shells(size_t nshells)
  {
    shells.clear();

    size_t n_per_shell = s_indices.size() / nshells;

    shell s;

    for (size_t j = 0; j < (nshells - 1); j ++) {
      s = shell(s_indices.begin() + j * n_per_shell,
		s_indices.begin() + (j + 1) * n_per_shell);
      shells.push_back(s);
    }

    s = shell(s_indices.begin() + (nshells - 1) * n_per_shell, s_indices.end());
    shells.push_back(s);

    size_t n_tot = 0;
    for(size_t i = 0; i < shells.size(); i++) {
      n_tot += shells[i].size();
    }

    CCTBX_ASSERT(s_indices.size() == n_tot);

  }

  float
  resolutionizer::isig(void)
  {
    unmerged_reflections_iterator uri;
    float result = 0.0;
    int unique = 0;

    for (uri = ur.begin(); uri != ur.end(); ++uri) {
      merged_isig mi = ext::merge(uri->second);
      result += mi.I / mi.sigI;
      unique += 1;
    }

    return result / unique;
  }

  scitbx::af::shared<float>
  resolutionizer::isig_shells(void)
  {
    scitbx::af::shared<float> result;

    CCTBX_ASSERT(shells.size() != 0);

    for (size_t i = 0; i < shells.size(); i++) {
      float total = 0.0;
      int n = 0;
      for (size_t j = 0; j < shells[i].size(); j ++) {
        merged_isig mi = ext::merge(ur[shells[i][j]]);
        total += mi.I / mi.sigI;
        n += 1;
      }
      result.push_back(total / n);
    }

    return result;
  }



} // namespace x2tbx::ext

BOOST_PYTHON_MODULE(x2tbx_ext)
{
  x2tbx::ext::init_module();
  boost::python::class_<x2tbx::resolutionizer>("resolutionizer")
    .def("set_unit_cell", & x2tbx::resolutionizer::set_unit_cell)
    .def("compare_resolution", & x2tbx::resolutionizer::compare_resolution)
    .def("setup", & x2tbx::resolutionizer::setup)
    .def("isig", & x2tbx::resolutionizer::isig)
    .def("sorted_indices", & x2tbx::resolutionizer::sorted_indices)
    .def("setup_shells", & x2tbx::resolutionizer::setup_shells)
    .def("isig_shells", & x2tbx::resolutionizer::isig_shells);
}
