/* x2tbx
 *
 * a toolbox to gracefully handle unmerged reflections for (in the first
 * instance) calculations in PyChef and resolution limits. N.B. will have
 * fundamental data structures:
 *
 * observation - float I, float sigI, float property, int flag
 *
 * unmerged_reflections - map(cctbx::miller::index, std::vector<observation>)
 *
 * though ideally want these nicely available from Python too (though that I
 * can live without.)
 *
 * first task try just implementing one calculation for e.g. resolution limits.
 *
 */

#include <boost/python.hpp>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <miller.h>
#include <uctbx.h>
#include <cctype>
#include <map>
#include <vector>

namespace cmil = cctbx::miller;
namespace cuc = cctbx::uctbx;

namespace x2tbx {

  typedef struct {
    float I;
    float sigI;
    float property;
    int flag;
  } observation;

  typedef struct {
    float I;
    float sigI;
  } merged_isig;

  typedef std::vector<observation> observation_list;

  typedef std::map<cctbx::miller::index<int>, observation_list> \
    unmerged_reflections;
  typedef std::map<cctbx::miller::index<int>, observation_list>::iterator \
    unmerged_reflections_iterator;

  typedef struct {
    unmerged_reflections ur;
    cuc::unit_cell unit_cell;

    void set_unit_cell(scitbx::af::tiny<double, 6> new_unit_cell);
    void setup(scitbx::af::const_ref<cmil::index<int> > const & indices,
               scitbx::af::const_ref<float> const & i_data,
               scitbx::af::const_ref<float> const & sigi_data);
    float isig(void);
  } resolutionizer;
}
