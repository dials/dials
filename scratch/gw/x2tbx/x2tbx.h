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

#include <miller.h>
#include <map>
#include <vector>

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
}
