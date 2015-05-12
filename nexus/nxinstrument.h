

#ifndef DIALS_NEXUS_NXINSTRUMENT_H
#define DIALS_NEXUS_NXINSTRUMENT_H

#include <dials/nexus/nxattenuator.h>
#include <dials/nexus/nxdetector.h>
#include <dials/nexus/serialize.h>

namespace dials { namespace nexus {

  class NXinstrument {
  public:

    boost::optional<NXattenuator> attenuator;
    boost::optional<NXdetector> detector;

  };

  template <>
  struct serialize<NXinstrument> {

    template <typename Handle>
    static
    NXinstrument load(const Handle &handle) {
      return NXinstrument();
    }

    template <typename Handle>
    static
    void dump(const NXinstrument &obj, Handle &handle) {

    }

  };

}} // namespace dials::nexus

#endif // DIALS_NEXUS_NXINSTRUMENT_H
