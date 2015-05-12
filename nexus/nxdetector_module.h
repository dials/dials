

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
    vec3<double> module_offset;
    vec3<double> fast_pixel_direction;
    vec3<double> slow_pixel_direction;

  };

}} // namespace dials::nexus

#endif // DIALS_NEXUS_NXDETECTOR_MODULE_H
