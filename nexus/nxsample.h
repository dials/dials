


#ifndef DIALS_NEXUS_NXSAMPLE_H
#define DIALS_NEXUS_NXSAMPLE_H

#include <string>
#include <boost/optional.hpp>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <dials/nexus/nxbeam.h>

namespace dials { namespace nexus {

  using scitbx::vec3;
  using scitbx::mat3;

  class NXsample {
  public:

    boost::optional<std::string> name;
    boost::optional<std::string> chemical_formula;
    boost::optional<double> temperature;
    af::shared<std::string> unit_cell_class;
    af::shared<std::string> unit_cell_group;
    boost::optional< vec3<double> > sample_orientation;
    af::shared< mat3<double> > orientation_matrix;
    af::shared< af::tiny<double,6> > unit_cell;
    boost::optional<NXbeam> beam;

  };

  template <>
  struct serialize<NXsample> {

    template <typename Handle>
    static
    NXsample load(const Handle &handle) {
      return NXsample();
    }

    template <typename Handle>
    static
    void dump(const NXsample &obj, Handle &handle) {

    }

  };

}} // namespace dials::nexus

#endif // DIALS_NEXUS_NXSAMPLE_H
