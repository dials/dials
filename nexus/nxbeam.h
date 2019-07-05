
#ifndef DIALS_NEXUS_NXBEAM_H
#define DIALS_NEXUS_NXBEAM_H

#include <boost/optional.hpp>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <scitbx/array_family/tiny.h>

namespace dials { namespace nexus {

  class NXbeam {
  public:
    boost::optional<double> incident_wavelength;
    boost::optional<double> flux;
    boost::optional<af::tiny<double, 4> > incident_polarization_stokes;
  };

}}  // namespace dials::nexus

#endif  // DIALS_NEXUS_NXBEAM_H
