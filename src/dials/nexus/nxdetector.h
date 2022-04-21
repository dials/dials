
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
    af::shared<double> dead_time;
    af::shared<double> count_time;
    af::shared<double> frame_time;
    af::shared<double> detector_readout_time;
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
    af::versa<double, af::c_grid<2> > angular_calibration;
    af::versa<double, af::c_grid<2> > flatfield;
    af::versa<double, af::c_grid<2> > flatfield_error;
    af::versa<int, af::c_grid<2> > pixel_mask;
    af::shared<NXdetector_module> module;
  };

  template <>
  struct serialize<NXdetector> {
    template <typename Handle>
    static NXdetector load(const Handle &handle) {
      NXdetector result;

      // Process the objects in the group
      for (std::size_t i = 0; i < handle.getNumObjs(); ++i) {
        // Get the name of the object
        std::string name = handle.getObjnameByIdx(i);

        switch (handle.getObjTypeByIdx(i)) {
        case H5G_GROUP: {
          H5::Group group = handle.openGroup(name);
          if (is_nx_class(group, "NXdetector_module")) {
            result.module.push_back(serialize<NXdetector_module>::load(group));
          }
        } break;

        case H5G_DATASET: {
          H5::DataSet dset = handle.openDataSet(name);
          if (name == "description") {
            result.description = serialize<std::string>::load(dset);
          } else if (name == "time_per_channel") {
            result.time_per_channel = serialize<std::string>::load(dset);
          } else if (name == "distance") {
            result.distance = serialize<double>::load(dset);
          } else if (name == "beam_centre_x") {
            result.beam_centre_x = serialize<double>::load(dset);
          } else if (name == "beam_centre_y") {
            result.beam_centre_y = serialize<double>::load(dset);
          } else if (name == "dead_time") {
            result.dead_time = serialize<af::shared<double> >::load(dset);
          } else if (name == "count_time") {
            result.count_time = serialize<af::shared<double> >::load(dset);
          } else if (name == "frame_time") {
            result.frame_time = serialize<af::shared<double> >::load(dset);
          } else if (name == "detector_readout_time") {
            result.detector_readout_time = serialize<af::shared<double> >::load(dset);
          } else if (name == "bit_depth_readout") {
            result.bit_depth_readout = serialize<int>::load(dset);
          } else if (name == "saturation_value") {
            result.saturation_value = serialize<int>::load(dset);
          } else if (name == "sensor_material") {
            result.sensor_material = serialize<std::string>::load(dset);
          } else if (name == "sensor_thickness") {
            result.sensor_thickness = serialize<double>::load(dset);
          } else if (name == "threshold_energy") {
            result.threshold_energy = serialize<double>::load(dset);
          } else if (name == "type") {
            result.type = serialize<std::string>::load(dset);
          } else if (name == "gain_setting") {
            result.gain_setting = serialize<std::string>::load(dset);
          } else if (name == "angular_calibration_applied") {
            result.angular_calibration_applied = serialize<bool>::load(dset);
          } else if (name == "flatfield_applied") {
            result.flatfield_applied = serialize<bool>::load(dset);
          } else if (name == "pixel_mask_applied") {
            result.pixel_mask_applied = serialize<bool>::load(dset);
          } else if (name == "countrate_correction_applied_applied") {
            result.countrate_correction_applied_applied = serialize<bool>::load(dset);
          } else if (name == "angular_calibration") {
            result.angular_calibration =
              serialize<af::versa<double, af::c_grid<2> > >::load(dset);
          } else if (name == "flatfield") {
            result.flatfield =
              serialize<af::versa<double, af::c_grid<2> > >::load(dset);
          } else if (name == "flatfield_error") {
            result.flatfield_error =
              serialize<af::versa<double, af::c_grid<2> > >::load(dset);
          } else if (name == "pixel_mask") {
            result.pixel_mask = serialize<af::versa<int, af::c_grid<2> > >::load(dset);
          }
        } break;

        default:
          break;
        };
      }

      // Return the NXdetector object
      return result;
    }

    template <typename Handle>
    static void dump(const NXdetector &obj, Handle &handle) {}
  };

}}  // namespace dials::nexus

#endif  // DIALS_NEXUS_NXDETECTOR_H
