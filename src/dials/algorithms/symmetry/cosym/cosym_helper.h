#include <dials/array_family/scitbx_shared_and_versa.h>
#include <cctbx/miller.h>
#include <unordered_map>
#include <tuple>
#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>

#include <functional>
#include <cstddef>

inline void hash_combine(std::size_t& seed, int value) {
  seed ^= std::hash<int>{}(value) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}

struct MillerIndexHash {
  std::size_t operator()(const cctbx::miller::index<>& idx) const noexcept {
    std::size_t seed = 0;
    hash_combine(seed, idx[0]);
    hash_combine(seed, idx[1]);
    hash_combine(seed, idx[2]);
    return seed;
  }
};

class matcher {
public:
  matcher(scitbx::af::const_ref<cctbx::miller::index<> > const& miller_indices) {
    // Create the lookup map for the miller indices of a dataset on initialisation,
    // then future calls can just match against this.
    lookup_map_.reserve(miller_indices.size());
    for (std::size_t i = 0; i < miller_indices.size(); i++) {
      lookup_map_.emplace(miller_indices[i], i);
    }
  }
  boost::python::tuple match(
    scitbx::af::const_ref<cctbx::miller::index<> > const& miller_indices) {
    // Return the two arrays of matching indices (in order) for a dataset
    // against the original dataset used in initialisation.
    scitbx::af::shared<std::size_t> indices_0{};
    scitbx::af::shared<std::size_t> indices_1{};

    indices_0.reserve(miller_indices.size());
    indices_1.reserve(miller_indices.size());

    for (std::size_t i = 0; i < miller_indices.size(); ++i) {
      auto it = lookup_map_.find(miller_indices[i]);
      if (it != lookup_map_.end()) {
        indices_0.push_back(it->second);
        indices_1.push_back(i);
      }
    }
    return boost::python::make_tuple(indices_0, indices_1);
  }

private:
  std::unordered_map<cctbx::miller::index<>, std::size_t, MillerIndexHash>
    lookup_map_{};
};
