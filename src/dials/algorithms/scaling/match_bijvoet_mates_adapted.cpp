#include <dials/algorithms/scaling/match_bijvoet_mates_adapted.h>
#include <cctbx/error.h>
#include <set>

namespace cctbx { namespace miller {

  void match_bijvoet_mates_adapted::match_(
    sgtbx::reciprocal_space::asu const& asu,
    bool /*assert_is_unique_set_under_symmetry*/) {
    typedef std::
      map<cctbx::miller::index<>, scitbx::af::shared<std::size_t>, fast_less_than<>>
        lookup_type;
    lookup_type lookup;

    // ✅ Build groups
    for (std::size_t i = 0; i < miller_indices_.size(); i++) {
      lookup[miller_indices_[i]].push_back(i);
    }

    std::set<cctbx::miller::index<>> processed;

    for (auto const& kv : lookup) {
      cctbx::miller::index<> h = kv.first;
      if (processed.count(h) || processed.count(-h)) continue;

      auto it_plus = lookup.find(h);
      auto it_minus = lookup.find(-h);

      if (it_minus != lookup.end()) {
        group_pair gp;

        // Assign plus/minus using ASU (consistent ordering)
        int which = asu.which(h);

        if (which > 0) {
          gp.plus = it_plus->second;
          gp.minus = it_minus->second;
        } else {
          gp.plus = it_minus->second;
          gp.minus = it_plus->second;
        }

        pairs_.push_back(gp);
      } else {
        // no Friedel mate → singles
        int which = asu.which(h);
        if (which > 0)
          singles_[0].push_back(it_plus->second);
        else
          singles_[1].push_back(it_plus->second);
      }

      processed.insert(h);
      processed.insert(-h);
    }
  }

  scitbx::af::shared<std::size_t> match_bijvoet_mates_adapted::singles(
    char plus_or_minus) const {
    std::size_t j = plus_or_minus_index_(plus_or_minus);
    scitbx::af::shared<std::size_t> flat;

    for (auto const& group : singles_[j]) {
      for (auto i : group)
        flat.push_back(i);
    }
    return flat;
  }

  scitbx::af::shared<std::size_t>
  match_bijvoet_mates_adapted::pairs_hemisphere_selection(char plus_or_minus) const {
    std::size_t j = plus_or_minus_index_(plus_or_minus);
    scitbx::af::shared<std::size_t> result;

    for (auto const& gp : pairs_) {
      const auto& g = (j == 0) ? gp.plus : gp.minus;
      for (auto i : g)
        result.push_back(i);
    }
    return result;
  }

  scitbx::af::shared<cctbx::miller::index<>>
  match_bijvoet_mates_adapted::miller_indices_in_hemisphere(char plus_or_minus) const {
    std::size_t j = plus_or_minus_index_(plus_or_minus);
    scitbx::af::shared<cctbx::miller::index<>> result;

    for (auto const& gp : pairs_) {
      const auto& g = (j == 0) ? gp.plus : gp.minus;

      if (!g.empty()) {
        // ✅ take representative (first index)
        result.push_back(miller_indices_[g[0]]);
      }
    }
    return result;
  }

  std::size_t match_bijvoet_mates_adapted::plus_or_minus_index_(
    char plus_or_minus) const {
    CCTBX_ASSERT(plus_or_minus == '+' || plus_or_minus == '-');
    return (plus_or_minus == '-') ? 1 : 0;
  }

}}  // namespace cctbx::miller
