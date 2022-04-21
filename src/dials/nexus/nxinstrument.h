

#ifndef DIALS_NEXUS_NXINSTRUMENT_H
#define DIALS_NEXUS_NXINSTRUMENT_H

#include <dials/nexus/nxattenuator.h>
#include <dials/nexus/nxdetector.h>
#include <dials/nexus/serialize.h>

namespace dials { namespace nexus {

  class NXinstrument {
  public:
    af::shared<NXattenuator> attenuator;
    af::shared<NXdetector> detector;
  };

  template <>
  struct serialize<NXinstrument> {
    template <typename Handle>
    static NXinstrument load(const Handle &handle) {
      NXinstrument result;

      // Process the objects in the group
      for (std::size_t i = 0; i < handle.getNumObjs(); ++i) {
        // Get the name of the object
        std::string name = handle.getObjnameByIdx(i);

        switch (handle.getObjTypeByIdx(i)) {
        case H5G_GROUP: {
          H5::Group group = handle.openGroup(name);
          if (is_nx_class(group, "NXattenuator")) {
            result.attenuator.push_back(serialize<NXattenuator>::load(group));
          } else if (is_nx_class(group, "NXdetector")) {
            result.detector.push_back(serialize<NXdetector>::load(group));
          }
        } break;

        default:
          break;
        };
      }

      // Return the NXinstrument object
      return result;
    }

    template <typename Handle>
    static void dump(const NXinstrument &obj, Handle &handle) {}
  };

}}  // namespace dials::nexus

#endif  // DIALS_NEXUS_NXINSTRUMENT_H
