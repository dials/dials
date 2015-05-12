
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
    boost::optional<NXinstrument> instrument;
    boost::optional<NXsample> sample;
    boost::optional<NXdata> data;

  };


  template <>
  class serialize<NXmx> {
  public:

    template <typename Handle>
    static
    NXmx load(const Handle &handle) {
      return NXmx();
    }

    template <typename Handle>
    static
    void dump(const NXmx &obj, Handle &handle) {

    }

  };


}} // namespace dials::nexus

#endif // DIALS_NEXUS_NXMX_H
