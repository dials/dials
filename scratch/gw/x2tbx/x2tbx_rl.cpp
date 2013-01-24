#include <x2tbx.h>

namespace x2tbx {
  ReflectionList::ReflectionList(void) { };
  ReflectionList::~ReflectionList(void) { };

  void
  ReflectionList::setup(miller_index_list_type indices,
                        float_value_list_type i_data,
                        float_value_list_type sigi_data)
  {
    CCTBX_ASSERT(indices.size() == i_data.size());
    CCTBX_ASSERT(i_data.size() == sigi_data.size());

    i_sig_type o;

    for (size_t i = 0; i < i_data.size(); i++) {
      o = i_sig_type(i_data[i], sigi_data[i]);
      reflections[indices[i]].add(o);
    }

    // gather unique reflection indices

    std::map<miller_index_type, ObservationList>::iterator ri;

    for (ri = reflections.begin(); ri != reflections.end(); ++ri) {
      unique_indices.push_back(ri->first);
    }

  }

  miller_index_list_type
  ReflectionList::get_indices(void)
  {
    return unique_indices;
  }

  void
  ReflectionList::set_unit_cell(scitbx::af::tiny<double, 6> unit_cell_params)
  {
    unit_cell = cuc::unit_cell(unit_cell_params);
  }

  void
  ReflectionList::setup_resolution_shells(size_t n_shells)
  {
    shells.clear();
    high_limits.clear();
    low_limits.clear();

    std::sort(unique_indices.begin(), unique_indices.end(),
              sorter_by_resolution(unit_cell));

    size_t n_per_shell = unique_indices.size() / n_shells;

    miller_index_list_type s;
    miller_index_list_type::iterator shell_start, shell_end;

    for (size_t j = 0; j < (n_shells - 1); j ++) {
      shell_start = unique_indices.begin() + j * n_per_shell;
      shell_end = unique_indices.begin() + (j + 1) * n_per_shell;
      s = miller_index_list_type(shell_start, shell_end);
      shells.push_back(s);

      // and update the resolution shells
      high_limits.push_back(unit_cell.d(s[0]));
      low_limits.push_back(unit_cell.d(s[s.size() - 1]));
    }

    s = miller_index_list_type(unique_indices.begin() +
                               (n_shells - 1) * n_per_shell,
                               unique_indices.end());
    shells.push_back(s);
    high_limits.push_back(unit_cell.d(s[0]));
    low_limits.push_back(unit_cell.d(s[s.size() - 1]));

    size_t n_tot = 0;
    for(size_t i = 0; i < shells.size(); i++) {
      n_tot += shells[i].size();
    }

    CCTBX_ASSERT(unique_indices.size() == n_tot);

  }

  miller_index_list_type
  ReflectionList::get_shell(size_t shell)
  {
    CCTBX_ASSERT(shell < shells.size());
    return shells[shell];
  }

  scitbx::af::shared<float>
  ReflectionList::shell_high_limits(void)
  {
    return high_limits;
  }

  scitbx::af::shared<float>
  ReflectionList::shell_low_limits(void)
  {
    return low_limits;
  }

  void
  ReflectionList::merge(void)
  {
    CCTBX_ASSERT(reflections.size() > 0);

    std::map<miller_index_type, ObservationList>::iterator r;

    for(r = reflections.begin(); r != reflections.end(); ++r) {
      (r->second).merge();
    }
  }

  float
  ReflectionList::i_sigma(void)
  {
    CCTBX_ASSERT(reflections.size() > 0);

    i_sig_type i_s;
    float result = 0.0;
    std::map<miller_index_type, ObservationList>::iterator r;

    for(r = reflections.begin(); r != reflections.end(); ++r) {
      i_s = (r->second).get_i_sigma();
      result += i_s[0] / i_s[1];
    }

    return result / reflections.size();
  }

  float
  ReflectionList::rmerge(void)
  {
    CCTBX_ASSERT(reflections.size() > 0);

    float r_sum = 0.0, i_sum = 0;
    std::map<miller_index_type, ObservationList>::iterator r;

    for(r = reflections.begin(); r != reflections.end(); ++r) {
      ObservationList o = r->second;

      r_sum += o.rmerge();
      i_sum += o.get_i_sigma()[0] * o.multiplicity();
    }

    return r_sum / i_sum;
  }

  scitbx::af::shared<float>
  ReflectionList::i_sigma_shells(void)
  {
    CCTBX_ASSERT(reflections.size() > 0);
    CCTBX_ASSERT(shells.size() > 0);

    scitbx::af::shared<float> result;

    for (size_t shell = 0; shell < shells.size(); shell ++) {
      i_sig_type i_s;
      float shell_result = 0.0;
      miller_index_list_type indices = shells[shell];
      for (size_t i = 0; i < indices.size(); i++) {
        i_s = reflections[indices[i]].get_i_sigma();
        shell_result += i_s[0] / i_s[1];
      }

      result.push_back(shell_result / indices.size());
    }

    return result;
  }

  scitbx::af::shared<float>
  ReflectionList::rmerge_shells(void)
  {
    CCTBX_ASSERT(reflections.size() > 0);
    CCTBX_ASSERT(shells.size() > 0);

    scitbx::af::shared<float> result;

    for (size_t shell = 0; shell < shells.size(); shell ++) {
      float r_sum = 0.0, i_sum = 0;
      miller_index_list_type indices = shells[shell];
      for (size_t i = 0; i < indices.size(); i++) {
        ObservationList o = reflections[indices[i]];
        r_sum += o.rmerge();
        i_sum += o.get_i_sigma()[0] * o.multiplicity();
      }

      result.push_back(r_sum / i_sum);
    }

    return result;
  }
}
