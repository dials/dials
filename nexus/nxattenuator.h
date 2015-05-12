
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
    static
    NXattenuator load(const Handle &handle) {
      return NXattenuator();
    }

    template <typename Handle>
    static
    void dump(const NXattenuator &obj, Handle &handle) {

    }

  };

}} // namespace dials::nexus

#endif // DIALS_NEXUS_NXATTENUATOR_H
