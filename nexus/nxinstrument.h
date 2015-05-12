

#ifndef DIALS_NEXUS_NXINSTRUMENT_H
#define DIALS_NEXUS_NXINSTRUMENT_H

#include <dials/nexus/nxattenuator.h>
#include <dials/nexus/nxdetector.h>

namespace dials { namespace nexus {

  class NXinstrument {
  public:

    boost::optional<NXattenuator> attenuator;
    boost::optional<NXdetector> detector;

  };

}} // namespace dials::nexus

#endif // DIALS_NEXUS_NXINSTRUMENT_H
