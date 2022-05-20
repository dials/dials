

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

  using scitbx::mat3;
  using scitbx::vec3;

  class NXsample {
  public:
    boost::optional<std::string> name;
    boost::optional<std::string> chemical_formula;
    boost::optional<double> temperature;
    af::shared<std::string> unit_cell_class;
    af::shared<std::string> unit_cell_group;
    boost::optional<vec3<double> > sample_orientation;
    af::shared<mat3<double> > orientation_matrix;
    af::shared<af::tiny<double, 6> > unit_cell;
    boost::optional<NXbeam> beam;
  };

  template <>
  struct serialize<NXsample> {
    template <typename Handle>
    static NXsample load(const Handle &handle) {
      NXsample result;

      // Process the objects in the group
      for (std::size_t i = 0; i < handle.getNumObjs(); ++i) {
        // Get the name of the object
        std::string name = handle.getObjnameByIdx(i);

        switch (handle.getObjTypeByIdx(i)) {
        case H5G_GROUP: {
          H5::Group group = handle.openGroup(name);
          if (is_nx_class(group, "NXbeam")) {
            result.beam.push_back(serialize<NXbeam>::load(group));
          }
        } break;

        case H5G_DATASET: {
          H5::DataSet dset = handle.openDataSet(name);
          if (name == "name") {
            result.name = serialize<std::string>::load(dset);
          } else if (name == "chemical_formula") {
            result.chemical_formula = serialize<std::string>::load(dset);
          } else if (name == "temperature") {
            result.temperature = serialize<double>::load(dset);
          } else if (name == "unit_cell_class") {
            result.unit_cell_class = serialize<af::shared<std::string> >::load(dset);
          } else if (name == "unit_cell_group") {
            result.unit_cell_group = serialize<af::shared<std::string> >::load(dset);
          } else if (name == "sample_orientation") {
            result.sample_orientation = serialize<vec3<double> >::load(dset);
          } else if (name == "orientation_matrix") {
            result.orientation_matrix =
              serialize<af::shared<mat3<double> > >::load(dset);
          } else if (name == "unit_cell") {
            result.unit_cell = serialize<af::shared<af::tiny<double, 6> > >::load(dset);
          }
        } break;

        default:
          break;
        };
      }

      // Return the NXsample object
      return result;
    }

    template <typename Handle>
    static void dump(const NXsample &obj, Handle &handle) {}
  };

}}  // namespace dials::nexus

#endif  // DIALS_NEXUS_NXSAMPLE_H
