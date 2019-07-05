
#ifndef DIALS_NEXUS_NXMX_H
#define DIALS_NEXUS_NXMX_H

#include <string>
#include <boost/optional.hpp>
#include <dials/nexus/serialize.h>
#include <dials/nexus/nxinstrument.h>
#include <dials/nexus/nxsample.h>
#include <dials/nexus/nxdata.h>

namespace dials { namespace nexus {

  class NXmx {
  public:
    boost::optional<std::string> title;
    boost::optional<std::string> start_time;
    boost::optional<std::string> end_time;

    af::shared<NXinstrument> instrument;
    af::shared<NXsample> sample;
  };

  template <>
  struct serialize<NXmx> {
    template <typename Handle>
    static NXmx load(const Handle &handle) {
      NXmx result;

      // Process the objects in the group
      for (std::size_t i = 0; i < handle.getNumObjs(); ++i) {
        // Get the name of the object
        std::string name = handle.getObjnameByIdx(i);

        switch (handle.getObjTypeByIdx(i)) {
        case H5G_GROUP: {
          H5::Group group = handle.openGroup(name);
          if (is_nx_class(group, "NXinstrument")) {
            result.instrument.push_back(serialize<NXinstrument>::load(group));
          } else if (is_nx_class(group, "NXsample")) {
            result.sample.push_back(serialize<NXsample>::load(group));
          }
        } break;

        case H5G_DATASET: {
          H5::DataSet dset = handle.openDataSet(name);
          if (name == "title") {
            result.title = serialize<std::string>::load(dset);
          } else if (name == "start_time") {
            result.start_time = serialize<std::string>::load(dset);
          } else if (name == "end_time") {
            result.end_time = serialize<std::string>::load(dset);
          }
        } break;

        default:
          break;
        };
      }

      // Return the NXmx object
      return result;
    }

    template <typename Handle>
    static void dump(const NXmx &obj, Handle &handle) {}
  };

}}  // namespace dials::nexus

#endif  // DIALS_NEXUS_NXMX_H
