#ifndef DX2_MODEL_SCAN_H
#define DX2_MODEL_SCAN_H
#include <nlohmann/json.hpp>
using json = nlohmann::json;

class Scan {
  // A class to represent the physical measurement, consisting of the number of
  // images, starting oscillation and a constant oscillation width between
  // sequential images. This class MUST NOT be modified during processing or
  // used to track additional metadata.
public:
  Scan() = default;
  Scan(std::array<int, 2> image_range, std::array<double, 2> oscillation);
  Scan(json scan_data);
  std::array<int, 2> get_image_range() const;
  std::array<double, 2> get_oscillation() const;
  json to_json() const;

protected:
  std::array<int, 2> image_range_{{0, 0}};
  int num_images_{0};
  double oscillation_width_{0.0};
  double oscillation_start_{0.0};
};

Scan::Scan(std::array<int, 2> image_range, std::array<double, 2> oscillation)
    : image_range_{image_range} {
  num_images_ = image_range_[1] - image_range_[0] + 1;
  oscillation_start_ = oscillation[0];
  oscillation_width_ = oscillation[1];
}

Scan::Scan(json scan_data) {
  // minimal required keys are image range and ["properties"]:"oscillation"
  std::vector<std::string> required_keys = {"image_range", "properties"};
  for (const auto &key : required_keys) {
    if (scan_data.find(key) == scan_data.end()) {
      throw std::invalid_argument("Key " + key +
                                  " is missing from the input scan JSON");
    }
  }
  image_range_ = scan_data["image_range"];
  num_images_ = image_range_[1] - image_range_[0] + 1;
  if (scan_data["properties"].find("oscillation") ==
      scan_data["properties"].end()) {
    throw std::invalid_argument(
        "Key oscillation is missing from the input scan['properties'] JSON");
  }
  std::vector<double> oscillation;
  json osc_array = scan_data["properties"]["oscillation"];
  // Just read the first two oscillation values as that is all that is needed
  if (osc_array.size() < 2) {
    throw std::invalid_argument(
        "scan['properties']['oscillation'] has <2 values in the scan JSON");
  }
  for (json::iterator it = osc_array.begin(); it != osc_array.begin() + 2;
       ++it) {
    oscillation.push_back(*it);
  }
  oscillation_start_ = oscillation[0];
  oscillation_width_ = oscillation[1] - oscillation[0];
}

json Scan::to_json() const {
  json scan_data;
  scan_data["image_range"] = image_range_;
  scan_data["batch_offset"] = 0; // We MUST NOT use batch offsets in dx2,
                                 // written out here for compatibility
  std::vector<double> oscillation_array(num_images_);
  for (int i = 0; i < num_images_; i++) {
    oscillation_array[i] = oscillation_start_ + oscillation_width_ * i;
  }
  scan_data["properties"] = {"oscillation", oscillation_array};
  // FIXME for dials compatibility, do we need exposure times and epochs arrays,
  // also valid_image_ranges empty dict?
  return scan_data;
}

std::array<int, 2> Scan::get_image_range() const { return image_range_; }

std::array<double, 2> Scan::get_oscillation() const {
  return {oscillation_start_, oscillation_width_};
}

#endif // DX2_MODEL_SCAN_H
