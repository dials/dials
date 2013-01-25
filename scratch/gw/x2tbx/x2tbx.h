/* x2tbx
 *
 * a toolbox to gracefully handle unmerged reflections for (in the first
 * instance) calculations in PyChef and resolution limits. N.B. will have
 * fundamental data structures:
 *
 * observation - float I, float sigI (at the moment) as a tiny array
 *
 * unmerged_reflections - map(cctbx::miller::index, std::vector<observation>)
 *
 * though ideally want these nicely available from Python too (though that I
 * can live without.)
 *
 * first task try just implementing one calculation for e.g. resolution limits.
 *
 */

#ifndef X2TBX_X2TBX_H
#define X2TBX_X2TBX_H

#include <boost/python.hpp>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <miller.h>
#include <uctbx.h>
#include <cctype>
#include <map>
#include <algorithm>
#include <vector>

namespace cmil = cctbx::miller;
namespace cuc = cctbx::uctbx;

namespace x2tbx {

  void init_module();

  struct merged_isig {
    float I;
    float sigI;
  };

  typedef scitbx::af::tiny<float, 2> i_sig_type;

  /**
   * ObservationList class - to store multiple observations of a given reflection
   * assuming all belong to one unique miller index.
   */

  class ObservationList {
  public:
    ObservationList();

    void add(i_sig_type);
    void merge();
    i_sig_type i_sigma();
    float total_i_sigma();
    size_t multiplicity();
    float rmerge();

  private:
    scitbx::af::shared<i_sig_type> observations_;
    float imean_, sigimean_, total_i_sigi_;
  };

  typedef cmil::index<int> miller_index_type;
  typedef scitbx::af::shared<miller_index_type> miller_index_list_type;
  typedef scitbx::af::shared<float> float_value_list_type;

  /**
   * ReflectionList class - to store and manipulate an entire data set
   * consisting of unmerged observations of many unique reflections.
   */

  class ReflectionList {
  public:
    ReflectionList();

    void setup(miller_index_list_type,
               float_value_list_type,
               float_value_list_type);
    void set_unit_cell(scitbx::af::tiny<double, 6>);

    miller_index_list_type get_indices();
    void setup_resolution_shells(size_t);

    void merge();
    float i_sigma();
    float total_i_sigma();
    float rmerge();
    scitbx::af::shared<float> rmerge_shells();
    scitbx::af::shared<float> i_sigma_shells();
    scitbx::af::shared<float> total_i_sigma_shells();
    scitbx::af::shared<float> shell_high_limits();
    scitbx::af::shared<float> shell_low_limits();
    miller_index_list_type get_shell(size_t);

  private:
    cuc::unit_cell unit_cell_;
    std::map<miller_index_type, ObservationList> reflections_;
    std::vector<miller_index_list_type> shells_;
    miller_index_list_type unique_indices_;
    scitbx::af::shared<float> high_limits_;
    scitbx::af::shared<float> low_limits_;

  };

  /**
   * Magic tool used to sort reflection indices in terms of their resolution
   * order, from the smallest dspacings to the largest. For use in std::sort.
   * Requires cctbx::uctbx::unit_cell in constructor.
   */

  struct sorter_by_resolution {
    cuc::unit_cell unit_cell;
    sorter_by_resolution(cuc::unit_cell new_unit_cell):
      unit_cell(new_unit_cell) { }
    bool operator() (cmil::index<int> const & a,
                     cmil::index<int> const & b)
    {
      return unit_cell.d(a) < unit_cell.d(b);
    }
  };
}

#endif
