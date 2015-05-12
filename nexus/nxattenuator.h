
#ifndef DIALS_NEXUS_NXATTENUATOR_H
#define DIALS_NEXUS_NXATTENUATOR_H

#include <boost/optional.hpp>

namespace dials { namespace nexus {

  class NXattenuator {
  public:

    boost::optional<double> attenuator_transmission;

  };

}} // namespace dials::nexus

#endif // DIALS_NEXUS_NXATTENUATOR_H
