
#ifndef DIALS_NEXUS_NXATTENUATOR_H
#define DIALS_NEXUS_NXATTENUATOR_H

#include <boost/optional.hpp>

namespace dials { namespace nexus {

  class NXattenuator {
  public:
    boost::optional<double> attenuator_transmission;
  };

  template <>
  struct serialize<NXattenuator> {
    template <typename Handle>
    static NXattenuator load(const Handle &handle) {
      NXattenuator result;

      // Process the objects in the group
      for (std::size_t i = 0; i < handle.getNumObjs(); ++i) {
        // Get the name of the object
        std::string name = handle.getObjnameByIdx(i);

        switch (handle.getObjTypeByIdx(i)) {
        case H5G_DATASET: {
          H5::DataSet dset = handle.openDataSet(name);
          if (name == "attenuator_transmission") {
            result.attenuator_transmission = serialize<double>::load(dset);
          }
        } break;

        default:
          break;
        };
      }

      // Return the NXattenuator object
      return result;
    }

    template <typename Handle>
    static void dump(const NXattenuator &obj, Handle &handle) {}
  };

}}  // namespace dials::nexus

#endif  // DIALS_NEXUS_NXATTENUATOR_H
