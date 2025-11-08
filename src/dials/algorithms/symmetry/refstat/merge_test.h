#pragma once
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/miller/asu.h>
#include <scitbx/array_family/selections.h>
#include <scitbx/array_family/sort.h>
namespace dials { namespace algorithms { namespace symmetry { namespace refstat {
  using namespace cctbx;
  using namespace scitbx;
  using namespace cctbx::sgtbx;

  template <typename FloatType>
  class merge_test_class {
    const af::shared<miller::index<> > indices;
    const af::shared<double> Is;
    const af::shared<double> sigs;

  public:
    merge_test_class(const af::shared<miller::index<> >& indices,
                     const af::shared<FloatType>& Is,
                     const af::shared<FloatType>& sigs)
        : indices(indices), Is(Is), sigs(sigs) {}

    struct merge_stats_result {
      size_t inconsistent_count;
      FloatType r_int;
    };

    merge_stats_result merge_test(const space_group& sg) const {
      af::shared<miller::index<> > sg_indices = indices.deep_copy();
      miller::map_to_asu(sg.type(), sg.is_centric(), sg_indices.ref());
      af::shared<size_t> perm = af::sort_permutation(sg_indices.const_ref());
      dry_shelx_merger sm;
      sm.process(af::select(sg_indices.const_ref(), perm.const_ref()).const_ref(),
                 af::select(Is.const_ref(), perm.const_ref()).const_ref(),
                 af::select(sigs.const_ref(), perm.const_ref()).const_ref());
      merge_stats_result rv;
      rv.inconsistent_count = sm.inconsistent_eq;
      rv.r_int = sm.r_int_num / sm.r_int_den;
      return rv;
    }

    struct sysabs_test_result {
      FloatType weak_I_sum, strong_I_sum, weak_sig_sq_sum, strong_sig_sq_sum;
      size_t weak_count, strong_count;
    };

    sysabs_test_result sysabs_test(const space_group& sg, FloatType scale = 1.0) const {
      sysabs_test_result rv;
      std::memset(&rv, 0, sizeof(rv));
      for (size_t i = 0; i < indices.size(); i++) {
        for (size_t j = 1; j < sg.order_z(); j++) {
          const rt_mx& rt = sg(j);
          if (rt.t().is_zero()) {
            continue;
          }
          if (indices[i] != rt.r() * indices[i]) {
            continue;
          }
          FloatType ps =
            rt.t().as_floating_point<FloatType>(scitbx::type_holder<FloatType>())
            * indices[i];
          if (std::abs(ps - std::round(ps)) < 0.01) {
            rv.strong_count++;
            rv.strong_I_sum += Is[i] * scale;
            rv.strong_sig_sq_sum += fn::pow2(sigs[i] * scale);
          } else {
            rv.weak_count++;
            rv.weak_I_sum += Is[i] * scale;
            rv.weak_sig_sq_sum += fn::pow2(sigs[i] * scale);
          }
        }
      }
      return rv;
    }

  protected:
    struct dry_shelx_merger {
      size_t inconsistent_eq;
      FloatType r_int_num, r_int_den;

      dry_shelx_merger() : inconsistent_eq(0), r_int_num(0), r_int_den(0) {}

      void process(const af::const_ref<miller::index<> >& unmerged_indices,
                   const af::const_ref<FloatType>& unmerged_data,
                   const af::const_ref<FloatType>& unmerged_sigmas) {
        if (unmerged_indices.size() == 0) {
          return;
        }
        size_t group_begin = 0;
        size_t group_end = 1;
        for (; group_end < unmerged_indices.size(); group_end++) {
          if (unmerged_indices[group_end] != unmerged_indices[group_begin]) {
            process_group(group_begin, group_end, unmerged_data, unmerged_sigmas);
            group_begin = group_end;
          }
        }
        process_group(group_begin, group_end, unmerged_data, unmerged_sigmas);
      }

      void process_group(size_t group_begin,
                         size_t group_end,
                         af::const_ref<FloatType> const& unmerged_data,
                         af::const_ref<FloatType> const& unmerged_sigmas) {
        size_t n = group_end - group_begin;
        if (n == 0) {
          return;
        }
        FloatType oss_sum = 0, w_sum = 0, i_wght_sum = 0;
        for (size_t i = 0; i < n; i++) {
          const size_t index = group_begin + i;
          const FloatType s =
            std::max(static_cast<FloatType>(1e-5), unmerged_sigmas[index]);
          const FloatType oss = fn::pow2(1. / s);
          const FloatType val = unmerged_data[index];
          const FloatType w = ((val > 3.0 * s) ? val * oss : 3.0 / s);
          oss_sum += oss;
          w_sum += w;
          i_wght_sum += w * val;
        }
        const FloatType mean = i_wght_sum / w_sum;
        FloatType sum_diff = 0, sum_i = 0, sum_diffs = 0;
        for (size_t i = 0; i < n; i++) {
          const FloatType val = unmerged_data[group_begin + i];
          const FloatType diff = val - mean;
          sum_diff += fn::absolute(diff);
          sum_i += val;
          sum_diffs += fn::pow2(diff);
        }
        FloatType sig = std::sqrt(1. / oss_sum);
        if (n > 1) {
          r_int_num += sum_diff;
          r_int_den += sum_i;
          const FloatType sig_int =
            sum_diff / (n * sqrt(static_cast<FloatType>(n) - 1.0));
          if (sig_int > sig) {
            if (sig_int > 5 * sig) {
              inconsistent_eq++;
            }
          }
        }
      }
    };  // merge_test
  };

}}}}  // namespace dials::algorithms::symmetry::refstat
