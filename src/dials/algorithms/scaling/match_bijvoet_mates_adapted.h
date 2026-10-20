#ifndef CCTBX_MILLER_MATCH_BIJVOET_MATES_ADAPTED_H
#define CCTBX_MILLER_MATCH_BIJVOET_MATES_ADAPTED_H

#include <cctbx/miller/match.h>
#include <cctbx/sgtbx/reciprocal_space_asu.h>
#include <map>

namespace cctbx { namespace miller {

  struct group_pair {
    scitbx::af::shared<std::size_t> plus;
    scitbx::af::shared<std::size_t> minus;
  };

  class match_bijvoet_mates_adapted {
  public:
    match_bijvoet_mates_adapted(
      sgtbx::space_group_type const& sg_type,
      scitbx::af::shared<cctbx::miller::index<>> const& miller_indices,
      bool assert_is_unique_set_under_symmetry = true)
        : miller_indices_(miller_indices) {
      match_(sgtbx::reciprocal_space::asu(sg_type),
             assert_is_unique_set_under_symmetry);
    }

    scitbx::af::shared<std::size_t> singles(char plus_or_minus) const;

    scitbx::af::shared<std::size_t> pairs_hemisphere_selection(
      char plus_or_minus) const;

    scitbx::af::shared<cctbx::miller::index<>> miller_indices_in_hemisphere(
      char plus_or_minus) const;

    template <typename NumType>
    scitbx::af::shared<NumType> minus(af::const_ref<NumType> const& data) const {
      scitbx::af::shared<NumType> result;
      for (auto const& gp : pairs_) {
        NumType sum_plus = 0, sum_minus = 0;

        for (auto i : gp.plus)
          sum_plus += data[i];
        for (auto i : gp.minus)
          sum_minus += data[i];

        NumType mean_plus = sum_plus / gp.plus.size();
        NumType mean_minus = sum_minus / gp.minus.size();

        result.push_back(mean_plus - mean_minus);
      }
      return result;
    }

  protected:
    void match_(sgtbx::reciprocal_space::asu const& asu,
                bool assert_is_unique_set_under_symmetry);

    std::size_t plus_or_minus_index_(char plus_or_minus) const;

    scitbx::af::shared<cctbx::miller::index<>> miller_indices_;
    std::vector<group_pair> pairs_;
    af::tiny<std::vector<scitbx::af::shared<std::size_t>>, 2> singles_;
  };

}}  // namespace cctbx::miller

#endif
