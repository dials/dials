#include <x2tbx.h>

namespace x2tbx {
  ReflectionList::ReflectionList(void) { };

  /**
   * Set up the reflection list from arrays pulled (perhaps) from an MTZ file
   * - this will do all of the sorting of reflections into the hash table.
   * @param indices Miller indices, as extracted from reflection file (not uniq)
   * @param i_data the measured intensities
   * @param sigi_data the corresponding error estimates
   */

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
      reflections_[indices[i]].add(o);
    }

    // gather unique reflection indices

    std::map<miller_index_type, ObservationList>::iterator ri;

    for (ri = reflections_.begin(); ri != reflections_.end(); ++ri) {
      unique_indices_.push_back(ri->first);
    }

  }

  /**
   * Get the unique Miller indices
   * @returns flex array of (integer) Miller indices
   */

  miller_index_list_type
  ReflectionList::get_indices()
  {
    return unique_indices_;
  }

  /**
   * Set the unit cell parameters, used to sort the indices into resolution
   * order
   * @param unit cell parameters as tuple / scitbx::af::tiny<double, 6>
   */

  void
  ReflectionList::set_unit_cell(scitbx::af::tiny<double, 6> unit_cell_params)
  {
    unit_cell_ = cctbx::uctbx::unit_cell(unit_cell_params);
  }

  /**
   * Set up the resolution shells for analysis / binning - needed before
   * asking for per-shell statistics. Divides unique reflections into
   * equally sized bins, resolution ordered.
   * @param number of shells to use
   */

  void
  ReflectionList::setup_resolution_shells(size_t n_shells)
  {
    shells_.clear();
    high_limits_.clear();
    low_limits_.clear();

    std::sort(unique_indices_.begin(), unique_indices_.end(),
              sorter_by_resolution(unit_cell_));

    size_t n_per_shell = unique_indices_.size() / n_shells;

    miller_index_list_type s;
    miller_index_list_type::iterator shell_start, shell_end;

    for (size_t j = 0; j < (n_shells - 1); j ++) {
      shell_start = unique_indices_.begin() + j * n_per_shell;
      shell_end = unique_indices_.begin() + (j + 1) * n_per_shell;
      s = miller_index_list_type(shell_start, shell_end);
      shells_.push_back(s);

      // and update the resolution shells
      high_limits_.push_back(unit_cell_.d(s[0]));
      low_limits_.push_back(unit_cell_.d(s[s.size() - 1]));
    }

    s = miller_index_list_type(unique_indices_.begin() +
                               (n_shells - 1) * n_per_shell,
                               unique_indices_.end());
    shells_.push_back(s);
    high_limits_.push_back(unit_cell_.d(s[0]));
    low_limits_.push_back(unit_cell_.d(s[s.size() - 1]));

    size_t n_tot = 0;
    for(size_t i = 0; i < shells_.size(); i++) {
      n_tot += shells_[i].size();
    }

    CCTBX_ASSERT(unique_indices_.size() == n_tot);

  }

  /**
   * Get the unique Miller indices for a given shell
   * @returns flex array of (integer) Miller indices
   */

  miller_index_list_type
  ReflectionList::get_shell(size_t shell)
  {
    CCTBX_ASSERT(shell < shells_.size());
    return shells_[shell];
  }

  /**
   * Get the high resolution limits for the resolution shells. Computed from
   * Miller index of highest resolution reflection in shell.
   * @returns scitbx array of high limits in Angstroms
   */

  scitbx::af::shared<float>
  ReflectionList::shell_high_limits()
  {
    return high_limits_;
  }

  /**
   * Get the low resolution limits for the resolution shells. Computed from
   * Miller index of lowest resolution reflection in shell.
   * @returns scitbx array of low limits in Angstroms
   */

  scitbx::af::shared<float>
  ReflectionList::shell_low_limits()
  {
    return low_limits_;
  }

  /**
   * Merge all of the observations for each unique reflection, required before
   * asking for any statistics, not needed to be repeated unless reflection
   * list modified.
   */

  void
  ReflectionList::merge()
  {
    CCTBX_ASSERT(reflections_.size() > 0);

    std::map<miller_index_type, ObservationList>::iterator r;

    for(r = reflections_.begin(); r != reflections_.end(); ++r) {
      (r->second).merge();
    }
  }

  /**
   * Get the overall merged I/sigma for the reflection list.
   * @returns Mn(I/sigma)
   */

  float
  ReflectionList::i_sigma()
  {
    CCTBX_ASSERT(reflections_.size() > 0);

    i_sig_type i_s;
    float result = 0.0;
    std::map<miller_index_type, ObservationList>::iterator r;

    for(r = reflections_.begin(); r != reflections_.end(); ++r) {
      i_s = (r->second).i_sigma();
      result += i_s[0] / i_s[1];
    }

    return result / reflections_.size();
  }

  /**
   * Get the overall unmerged I/sigma for the reflection list.
   * @returns I/sigma
   */

  float
  ReflectionList::total_i_sigma()
  {
    CCTBX_ASSERT(reflections_.size() > 0);

    i_sig_type i_s;
    float i_tot = 0.0;
    size_t n_tot = 0;
    std::map<miller_index_type, ObservationList>::iterator r;
    ObservationList o;

    for(r = reflections_.begin(); r != reflections_.end(); ++r) {
      o = (r->second);
      i_tot += o.total_i_sigma();
      n_tot += o.multiplicity();
    }

    return i_tot / n_tot;
  }

  /**
   * Get the overall Rmerge for the reflection list.
   * @returns rmerge
   */

  float
  ReflectionList::rmerge()
  {
    CCTBX_ASSERT(reflections_.size() > 0);

    float r_sum = 0.0, i_sum = 0;
    std::map<miller_index_type, ObservationList>::iterator r;

    for(r = reflections_.begin(); r != reflections_.end(); ++r) {
      ObservationList o = r->second;

      r_sum += o.rmerge();
      i_sum += o.i_sigma()[0] * o.multiplicity();
    }

    return r_sum / i_sum;
  }

  /**
   * Get the overall merged I/sigma for the each resolution shell in the
   * reflection list.
   * @returns flex array of Mn(I/sigma)
   */

  scitbx::af::shared<float>
  ReflectionList::i_sigma_shells()
  {
    CCTBX_ASSERT(reflections_.size() > 0);
    CCTBX_ASSERT(shells_.size() > 0);

    scitbx::af::shared<float> result;

    for (size_t shell = 0; shell < shells_.size(); shell ++) {
      i_sig_type i_s;
      float shell_result = 0.0;
      miller_index_list_type indices = shells_[shell];
      for (size_t i = 0; i < indices.size(); i++) {
        i_s = reflections_[indices[i]].i_sigma();
        shell_result += i_s[0] / i_s[1];
      }

      result.push_back(shell_result / indices.size());
    }

    return result;
  }

  /**
   * Get the overall unmerged I/sigma for the each resolution shell in the
   * reflection list.
   * @returns flex array of I/sigma
   */

  scitbx::af::shared<float>
  ReflectionList::total_i_sigma_shells()
  {
    CCTBX_ASSERT(reflections_.size() > 0);
    CCTBX_ASSERT(shells_.size() > 0);

    scitbx::af::shared<float> result;
    ObservationList o;

    for (size_t shell = 0; shell < shells_.size(); shell ++) {
      float shell_i_tot = 0.0;
      size_t shell_n_tot = 0;
      miller_index_list_type indices = shells_[shell];
      for (size_t i = 0; i < indices.size(); i++) {
        o = reflections_[indices[i]];
        shell_i_tot += o.total_i_sigma();
        shell_n_tot += o.multiplicity();
      }

      result.push_back(shell_i_tot / shell_n_tot);
    }

    return result;
  }

  /**
   * Get the overall Rmerge for the each resolution shell in the reflection
   * list.
   * @returns flex array of rmerge values
   */

  scitbx::af::shared<float>
  ReflectionList::rmerge_shells()
  {
    CCTBX_ASSERT(reflections_.size() > 0);
    CCTBX_ASSERT(shells_.size() > 0);

    scitbx::af::shared<float> result;

    for (size_t shell = 0; shell < shells_.size(); shell ++) {
      float r_sum = 0.0, i_sum = 0;
      miller_index_list_type indices = shells_[shell];
      for (size_t i = 0; i < indices.size(); i++) {
        ObservationList o = reflections_[indices[i]];
        r_sum += o.rmerge();
        i_sum += o.i_sigma()[0] * o.multiplicity();
      }

      result.push_back(r_sum / i_sum);
    }

    return result;
  }
}
