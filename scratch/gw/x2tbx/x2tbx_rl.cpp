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
    std::sort(unique_indices.begin(), unique_indices.end(),
              sorter_by_resolution(unit_cell));

    size_t n_per_shell = unique_indices.size() / n_shells;

    miller_index_list_type s;

    for (size_t j = 0; j < (n_shells - 1); j ++) {
      s = shell(unique_indices.begin() + j * n_per_shell,
                unique_indices.begin() + (j + 1) * n_per_shell);
      shells.push_back(s);
    }

    s = shell(unique_indices.begin() + (n_shells - 1) * n_per_shell,
              unique_indices.end());
    shells.push_back(s);

    size_t n_tot = 0;
    for(size_t i = 0; i < shells.size(); i++) {
      n_tot += shells[i].size();
    }

    CCTBX_ASSERT(unique_indices.size() == n_tot);

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

    i_sig_type i_s;
    float r_sum = 0.0, i_sum = 0;
    std::map<miller_index_type, ObservationList>::iterator r;

    for(r = reflections.begin(); r != reflections.end(); ++r) {
      ObservationList o = r->second;

      r_sum += o.rmerge();
      i_sum += o.get_i_sigma()[0] * o.multiplicity();
    }

    return r_sum / i_sum;
  }
}
