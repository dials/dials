

#ifndef DIALS_NEXUS_NXDETECTOR_MODULE_H
#define DIALS_NEXUS_NXDETECTOR_MODULE_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>

namespace dials { namespace nexus {

  using scitbx::vec2;
  using scitbx::vec3;

  class NXdetector_module {
  public:
    vec2<int> data_origin;
    vec2<int> data_size;
    double fast_pixel_size;
    double slow_pixel_size;
    vec3<double> module_offset;
    vec3<double> fast_pixel_direction;
    vec3<double> slow_pixel_direction;
  };

  template <>
  struct serialize<NXdetector_module> {
    template <typename Handle>
    static NXdetector_module load(const Handle &handle) {
      NXdetector_module result;

      // Process the objects in the group
      for (std::size_t i = 0; i < handle.getNumObjs(); ++i) {
        // Get the name of the object
        std::string name = handle.getObjnameByIdx(i);

        switch (handle.getObjTypeByIdx(i)) {
        case H5G_DATASET: {
          H5::DataSet dset = handle.openDataSet(name);
          if (name == "data_origin") {
            result.data_origin = serialize<vec2<int> >::load(dset);
          } else if (name == "data_size") {
            result.data_size = serialize<vec2<int> >::load(dset);
          }
        } break;

        default:
          break;
        };
      }

      // Return the NXdetector_module object
      return result;
    }

    template <typename Handle>
    static void dump(const NXdetector_module &obj, Handle &handle) {}
  };

}}  // namespace dials::nexus

#endif  // DIALS_NEXUS_NXDETECTOR_MODULE_H
