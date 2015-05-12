
#ifndef DIALS_NEXUS_NXDETECTOR_H
#define DIALS_NEXUS_NXDETECTOR_H

#include <boost/optional.hpp>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/nexus/nxdetector_module.h>

namespace dials { namespace nexus {

  class NXdetector {
  public:

    boost::optional<std::string> description;
    boost::optional<std::string> time_per_channel;
    boost::optional<double> distance;
    boost::optional<double> beam_centre_x;
    boost::optional<double> beam_centre_y;
    boost::optional<double> dead_time;
    boost::optional<double> count_time;
    boost::optional<double> frame_time;
    boost::optional<double> detector_readout_time;
    boost::optional<int> bit_depth_readout;
    boost::optional<int> saturation_value;
    boost::optional<std::string> sensor_material;
    boost::optional<double> sensor_thickness;
    boost::optional<double> threshold_energy;
    boost::optional<std::string> type;
    boost::optional<std::string> gain_setting;
    boost::optional<bool> angular_calibration_applied;
    boost::optional<bool> flatfield_applied;
    boost::optional<bool> pixel_mask_applied;
    boost::optional<bool> countrate_correction_applied_applied;
    af::versa< double, af::c_grid<2> > angular_calibration;
    af::versa< double, af::c_grid<2> > flatfield;
    af::versa< double, af::c_grid<2> > flatfield_error;
    af::versa< int, af::c_grid<2> > pixel_mask;
    af::shared<NXdetector_module> module;

  };

  template <>
  struct serialize<NXdetector> {

    template <typename Handle>
    static
    NXdetector load(const Handle &handle) {
      return NXdetector();
    }

    template <typename Handle>
    static
    void dump(const NXdetector &obj, Handle &handle) {

    }

  };

}} // namespace dials::nexus

#endif // DIALS_NEXUS_NXDETECTOR_H
