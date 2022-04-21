/*
 * pychef.h
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  Author: Richard Gildea
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_PYCHEF_H
#define DIALS_PYCHEF_H

#include <cctbx/miller.h>
#include <cctbx/miller/asu.h>
#include <cctbx/miller/bins.h>
#include <cctbx/miller/sym_equiv.h>
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/reciprocal_space_asu.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/error.h>

namespace dials { namespace pychef {

  using namespace cctbx;

  enum PlusMinusCentricCode {
    PLUS = 1,
    MINUS = -1,
    CENTRIC = 0,
  };

  class ObservationGroup {
  public:
    ObservationGroup() {}

    ObservationGroup(cctbx::miller::index<> miller_index, bool centric)
        : miller_index_(miller_index), centric_(centric) {}

    void add_iplus(std::size_t iplus) {
      iplus_.push_back(iplus);
    }

    void add_iminus(std::size_t iminus) {
      iminus_.push_back(iminus);
    }

    scitbx::af::shared<std::size_t> iplus() {
      return iplus_;
    }

    scitbx::af::shared<std::size_t> iminus() {
      return iminus_;
    }

    cctbx::miller::index<> miller_index() {
      return miller_index_;
    }

    bool is_centric() {
      return centric_;
    }

  private:
    cctbx::miller::index<> miller_index_;
    scitbx::af::shared<std::size_t> iplus_;
    scitbx::af::shared<std::size_t> iminus_;
    bool centric_;
  };

  struct Observations {
    Observations(scitbx::af::const_ref<cctbx::miller::index<> > const &miller_indices,
                 sgtbx::space_group space_group,
                 bool anomalous_flag) {
      scitbx::af::shared<cctbx::miller::index<> > asu_indices(miller_indices.begin(),
                                                              miller_indices.end());
      cctbx::miller::map_to_asu(space_group.type(), anomalous_flag, asu_indices.ref());
      sgtbx::reciprocal_space::asu asu(space_group.type());

      std::size_t n_minus = 0;
      std::size_t n_plus = 0;

      for (std::size_t iref = 0; iref < asu_indices.size(); iref++) {
        cctbx::miller::index<> h_uniq = asu_indices[iref];
        int flag;
        miller::sym_equiv_indices sym_equiv(space_group, h_uniq);
        if (sym_equiv.is_centric()) {
          flag = CENTRIC;
        } else {
          int asu_which = asu.which(h_uniq);
          DIALS_ASSERT((asu_which == 1) || (asu_which == -1));
          if (asu_which == 1) {
            flag = PLUS;
          }
          if (asu_which != 1) {
            for (std::size_t i = 0; i < 3; i++) {
              h_uniq[i] *= -1;
            }
            flag = MINUS;
          }
        }
        ObservationGroup group;
        map_type::iterator it = observation_groups_.find(h_uniq);
        if (it == observation_groups_.end()) {
          group = ObservationGroup(h_uniq, flag == CENTRIC);
          observation_groups_[h_uniq] = group;
        } else {
          group = observation_groups_[h_uniq];
        }
        if (flag == MINUS) {
          n_minus++;
          group.add_iminus(iref);
        } else {
          n_plus++;
          group.add_iplus(iref);
        }
      }
    }

    typedef std::map<cctbx::miller::index<>, ObservationGroup> map_type;

    map_type observation_groups_;

    map_type observation_groups() {
      return observation_groups_;
    }
  };

  namespace accumulator {

    class CompletenessAccumulator {
    public:
      CompletenessAccumulator(af::const_ref<std::size_t> const &dose,
                              af::const_ref<double> const &d_star_sq,
                              cctbx::miller::binner const &binner,
                              int n_steps)
          : finalised_(false),
            dose_(dose.begin(), dose.end()),
            d_star_sq_(d_star_sq.begin(), d_star_sq.end()),
            binner_(binner),
            n_steps_(n_steps),
            iplus_count(af::c_grid<2>(binner.n_bins_used(), n_steps), 0.0),
            iminus_count(af::c_grid<2>(binner.n_bins_used(), n_steps), 0.0),
            ieither_count(af::c_grid<2>(binner.n_bins_used(), n_steps), 0.0),
            iboth_count(af::c_grid<2>(binner.n_bins_used(), n_steps), 0.0),
            iplus_comp_overall(n_steps, 0.0),
            iminus_comp_overall(n_steps, 0.0),
            ieither_comp_overall(n_steps, 0.0),
            iboth_comp_overall(n_steps, 0.0) {}

      void operator()(ObservationGroup group) {
        std::size_t dose_min_iplus = 1e8;
        std::size_t dose_min_iminus = 1e8;
        std::size_t i_bin;
        if (group.iplus().size()) {
          i_bin = binner_.get_i_bin(d_star_sq_[group.iplus()[0]]);
        } else {
          i_bin = binner_.get_i_bin(d_star_sq_[group.iminus()[0]]);
        }

        if (i_bin == 0) {
          // outside "used" bins
          return;
        }

        i_bin -= 1;

        for (std::size_t i = 0; i < group.iplus().size(); i++) {
          std::size_t dose_i = dose_[group.iplus()[i]];
          dose_min_iplus = std::min(dose_i, dose_min_iplus);
          if (group.is_centric()) {
            dose_min_iminus = std::min(dose_i, dose_min_iminus);
          }
        }

        for (std::size_t i = 0; i < group.iminus().size(); i++) {
          for (std::size_t j = 0; j < group.iplus().size(); j++) {
            DIALS_ASSERT(group.iminus()[i] != group.iplus()[j]);
          }
          std::size_t dose_i = dose_[group.iminus()[i]];
          dose_min_iminus = std::min(dose_i, dose_min_iminus);
        }

        if (dose_min_iplus < n_steps_) {
          iplus_count(i_bin, dose_min_iplus) += 1.0;
        }
        if (dose_min_iminus < n_steps_) {
          iminus_count(i_bin, dose_min_iminus) += 1.0;
        }
        if (std::min(dose_min_iplus, dose_min_iminus) < n_steps_) {
          ieither_count(i_bin, std::min(dose_min_iplus, dose_min_iminus)) += 1.0;
        }
        if (std::max(dose_min_iplus, dose_min_iminus) < n_steps_) {
          iboth_count(i_bin, std::max(dose_min_iplus, dose_min_iminus)) += 1.0;
        }
      }

      void finalise(af::const_ref<std::size_t> counts_complete) {
        DIALS_ASSERT(!finalised_);

        DIALS_ASSERT(counts_complete.size() == binner_.n_bins_used());

        // accumulate as a function of time

        for (std::size_t i_bin = 0; i_bin < binner_.n_bins_used(); i_bin++) {
          for (std::size_t i_step = 1; i_step < n_steps_; i_step++) {
            iplus_count(i_bin, i_step) += iplus_count(i_bin, i_step - 1);
            iminus_count(i_bin, i_step) += iminus_count(i_bin, i_step - 1);
            ieither_count(i_bin, i_step) += ieither_count(i_bin, i_step - 1);
            iboth_count(i_bin, i_step) += iboth_count(i_bin, i_step - 1);
          }
        }

        // accumulate as a function of dose and resolution

        std::size_t tot_complete = 0;

        for (std::size_t i_bin = 0; i_bin < binner_.n_bins_used(); i_bin++) {
          std::size_t n_complete = counts_complete[i_bin];
          tot_complete += n_complete;
          double one_over_n_complete = 1.0 / static_cast<double>(n_complete);
          for (std::size_t i_step = 0; i_step < n_steps_; i_step++) {
            iplus_comp_overall[i_step] += iplus_count(i_bin, i_step);
            iminus_comp_overall[i_step] += iminus_count(i_bin, i_step);
            ieither_comp_overall[i_step] += ieither_count(i_bin, i_step);
            iboth_comp_overall[i_step] += iboth_count(i_bin, i_step);

            iplus_count(i_bin, i_step) *= one_over_n_complete;
            iminus_count(i_bin, i_step) *= one_over_n_complete;
            ieither_count(i_bin, i_step) *= one_over_n_complete;
            iboth_count(i_bin, i_step) *= one_over_n_complete;
          }
        }

        for (std::size_t i_step = 0; i_step < n_steps_; i_step++) {
          iplus_comp_overall[i_step] /= tot_complete;
          iminus_comp_overall[i_step] /= tot_complete;
          ieither_comp_overall[i_step] /= tot_complete;
          iboth_comp_overall[i_step] /= tot_complete;
        }

        finalised_ = true;
      }

      af::versa<double, af::c_grid<2> > iplus_completeness_bins() {
        DIALS_ASSERT(finalised_);
        return iplus_count;
      }

      af::versa<double, af::c_grid<2> > iminus_completeness_bins() {
        DIALS_ASSERT(finalised_);
        return iminus_count;
      }

      af::versa<double, af::c_grid<2> > ieither_completeness_bins() {
        DIALS_ASSERT(finalised_);
        return ieither_count;
      }

      af::versa<double, af::c_grid<2> > iboth_completeness_bins() {
        DIALS_ASSERT(finalised_);
        return iboth_count;
      }

      af::shared<double> iplus_completeness() {
        DIALS_ASSERT(finalised_);
        return iplus_comp_overall;
      }

      af::shared<double> iminus_completeness() {
        DIALS_ASSERT(finalised_);
        return iminus_comp_overall;
      }

      af::shared<double> ieither_completeness() {
        DIALS_ASSERT(finalised_);
        return ieither_comp_overall;
      }

      af::shared<double> iboth_completeness() {
        DIALS_ASSERT(finalised_);
        return iboth_comp_overall;
      }

    private:
      bool finalised_;

      af::shared<std::size_t> dose_;
      af::shared<double> d_star_sq_;
      cctbx::miller::binner const &binner_;
      std::size_t const n_steps_;

      af::versa<double, af::c_grid<2> > iplus_count, iminus_count, ieither_count,
        iboth_count;

      af::shared<double> iplus_comp_overall, iminus_comp_overall, ieither_comp_overall,
        iboth_comp_overall;
    };

    class RcpScpAccumulator {
    public:
      RcpScpAccumulator(af::const_ref<double> const &intensities,
                        af::const_ref<double> const &sigmas,
                        af::const_ref<std::size_t> const &dose,
                        af::const_ref<double> const &d_star_sq,
                        cctbx::miller::binner const &binner,
                        int n_steps)
          : finalised_(false),
            intensities_(intensities.begin(), intensities.end()),
            sigmas_(sigmas.begin(), sigmas.end()),
            dose_(dose.begin(), dose.end()),
            d_star_sq_(d_star_sq.begin(), d_star_sq.end()),
            binner_(binner),
            n_steps_(n_steps),
            A(af::c_grid<2>(binner.n_bins_used(), n_steps), 0.0),
            B(af::c_grid<2>(binner.n_bins_used(), n_steps), 0.0),
            isigma(af::c_grid<2>(binner.n_bins_used(), n_steps), 0.0),
            rcp_bins_(af::c_grid<2>(binner.n_bins_used(), n_steps), 0.0),
            scp_bins_(af::c_grid<2>(binner.n_bins_used(), n_steps), 0.0),
            count(af::c_grid<2>(binner.n_bins_used(), n_steps), 0),
            rcp_(n_steps, 0.0),
            scp_(n_steps, 0.0) {}

      void operator()(ObservationGroup group) {
        if (group.iplus().size()) {
          std::size_t i_bin = binner_.get_i_bin(d_star_sq_[group.iplus()[0]]);
          if (i_bin > binner_.n_bins_used()) {
          }
          DIALS_ASSERT(i_bin <= binner_.n_bins_used())(i_bin);
          if (i_bin == 0) {
            // outside "used" bins
            return;
          }
          accumulate(group.iplus(), i_bin - 1);
        }
        if (group.iminus().size()) {
          std::size_t i_bin = binner_.get_i_bin(d_star_sq_[group.iminus()[0]]);
          DIALS_ASSERT(i_bin <= binner_.n_bins_used())(i_bin);
          if (i_bin == 0) {
            // outside "used" bins
            return;
          }
          accumulate(group.iminus(), i_bin - 1);
        }
      }

    private:
      void accumulate(scitbx::af::shared<std::size_t> irefs, std::size_t i_bin) {
        DIALS_ASSERT(i_bin < binner_.n_bins_used());
        for (std::size_t i = 0; i < irefs.size(); i++) {
          std::size_t dose_i = dose_[irefs[i]];
          double I_i = intensities_[irefs[i]];
          double sigi_i = sigmas_[irefs[i]];
          for (std::size_t j = i + 1; j < irefs.size(); j++) {
            std::size_t dose_j = dose_[irefs[j]];
            double I_j = intensities_[irefs[j]];
            double sigi_j = sigmas_[irefs[j]];
            double A_part = std::fabs(I_i - I_j);
            double B_part = 0.5 * std::fabs(I_i + I_j);
            std::size_t dose_0 = std::max(dose_i, dose_j);
            DIALS_ASSERT(dose_0 < n_steps_);
            A(i_bin, dose_0) += A_part;
            B(i_bin, dose_0) += B_part;
            isigma(i_bin, dose_0) += ((I_i / sigi_i) + (I_j / sigi_j));
            count(i_bin, dose_0) += 2;
          }
        }
      }

    public:
      void finalise() {
        DIALS_ASSERT(!finalised_);

        // accumulate as a function of time

        for (std::size_t i_bin = 0; i_bin < binner_.n_bins_used(); i_bin++) {
          for (std::size_t i_step = 1; i_step < n_steps_; i_step++) {
            A(i_bin, i_step) += A(i_bin, i_step - 1);
            B(i_bin, i_step) += B(i_bin, i_step - 1);
            isigma(i_bin, i_step) += isigma(i_bin, i_step - 1);
            count(i_bin, i_step) += count(i_bin, i_step - 1);
          }
        }

        // accumulate as a function of dose and resolution

        for (std::size_t i_step = 0; i_step < n_steps_; i_step++) {
          double overall_top = 0;
          double overall_bottom = 0;
          double scp_overall = 0;

          for (std::size_t i_bin = 0; i_bin < binner_.n_bins_used(); i_bin++) {
            double top = A(i_bin, i_step);
            double bottom = B(i_bin, i_step);

            overall_top += top;
            overall_bottom += bottom;

            double rcp_tmp = 0;
            double scp_tmp = 0;

            if (bottom > 0) {
              rcp_tmp = top / bottom;
              if (count(i_bin, i_step) > 100) {
                double isig_tmp = isigma(i_bin, i_step) / count(i_bin, i_step);
                scp_tmp = rcp_tmp / (1.1284 / isig_tmp);
              }
            }
            rcp_bins_(i_bin, i_step) = rcp_tmp;
            scp_bins_(i_bin, i_step) = scp_tmp;
            scp_overall += scp_tmp / binner_.n_bins_used();
          }

          rcp_[i_step] = (overall_bottom > 0) ? (overall_top / overall_bottom) : 0.0;
          scp_[i_step] = scp_overall;
        }

        finalised_ = true;
      }

      af::versa<double, af::c_grid<2> > rcp_bins() {
        DIALS_ASSERT(finalised_);
        return rcp_bins_;
      }

      af::versa<double, af::c_grid<2> > scp_bins() {
        DIALS_ASSERT(finalised_);
        return scp_bins_;
      }

      af::shared<double> rcp() {
        DIALS_ASSERT(finalised_);
        return rcp_;
      }

      af::shared<double> scp() {
        DIALS_ASSERT(finalised_);
        return scp_;
      }

    private:
      bool finalised_;

      af::shared<double> intensities_;
      af::shared<double> sigmas_;
      af::shared<std::size_t> dose_;
      af::shared<double> d_star_sq_;
      cctbx::miller::binner const &binner_;
      std::size_t const n_steps_;

      af::versa<double, af::c_grid<2> > A, B, isigma, rcp_bins_, scp_bins_;

      af::versa<std::size_t, af::c_grid<2> > count;

      af::shared<double> rcp_, scp_;
    };

    class RdAccumulator {
    public:
      RdAccumulator(af::const_ref<double> const &intensities,
                    af::const_ref<std::size_t> const &dose,
                    int n_steps)
          : finalised_(false),
            intensities_(intensities.begin(), intensities.end()),
            dose_(dose.begin(), dose.end()),
            n_steps_(n_steps),
            rd_top(n_steps, 0.0),
            rd_bottom(n_steps, 0.0),
            rd_(n_steps, 0.0) {}

      void operator()(ObservationGroup group) {
        if (group.iplus().size()) {
          accumulate(group.iplus());
        }
        if (group.iminus().size()) {
          accumulate(group.iminus());
        }
      }

    private:
      void accumulate(scitbx::af::shared<std::size_t> irefs) {
        for (std::size_t i = 0; i < irefs.size(); i++) {
          int dose_i = dose_[irefs[i]];
          double I_i = intensities_[irefs[i]];
          for (std::size_t j = i + 1; j < irefs.size(); j++) {
            int dose_j = dose_[irefs[j]];
            double I_j = intensities_[irefs[j]];
            std::size_t d_dose = std::abs(dose_i - dose_j);
            DIALS_ASSERT(d_dose < n_steps_);
            rd_top[d_dose] += std::fabs(I_i - I_j);
            rd_bottom[d_dose] += 0.5 * (I_i + I_j);
          }
        }
      }

    public:
      void finalise() {
        DIALS_ASSERT(!finalised_);

        for (std::size_t i_step = 0; i_step < n_steps_; i_step++) {
          if (rd_bottom[i_step] > 0) {
            rd_[i_step] = rd_top[i_step] / rd_bottom[i_step];
          }
        }

        finalised_ = true;
      }

      af::shared<double> rd() {
        DIALS_ASSERT(finalised_);
        return rd_;
      }

    private:
      bool finalised_;

      af::shared<double> intensities_;
      af::shared<std::size_t> dose_;
      std::size_t const n_steps_;

      af::shared<double> rd_top, rd_bottom, rd_;
    };

  }  // namespace accumulator

  class ChefStatistics {
  public:
    ChefStatistics(scitbx::af::const_ref<cctbx::miller::index<> > const &miller_indices,
                   af::const_ref<double> const &intensities,
                   af::const_ref<double> const &sigmas,
                   af::const_ref<double> const &d_star_sq,
                   af::const_ref<std::size_t> const &dose,
                   af::const_ref<std::size_t> counts_complete,
                   cctbx::miller::binner const &binner,
                   sgtbx::space_group space_group,
                   bool anomalous_flag,
                   int n_steps)
        : observations(miller_indices, space_group, anomalous_flag),
          completeness_accumulator(dose, d_star_sq, binner, n_steps),
          rcp_scp_accumulator(intensities, sigmas, dose, d_star_sq, binner, n_steps),
          rd_accumulator(intensities, dose, n_steps) {
      typedef Observations::map_type map_t;
      map_t groups = observations.observation_groups();

      for (map_t::iterator it = groups.begin(); it != groups.end(); it++) {
        ObservationGroup group = it->second;
        completeness_accumulator(group);
        rcp_scp_accumulator(group);
        rd_accumulator(group);
      }

      completeness_accumulator.finalise(counts_complete);
      rcp_scp_accumulator.finalise();
      rd_accumulator.finalise();
    }

    af::versa<double, af::c_grid<2> > iplus_completeness_bins() {
      return completeness_accumulator.iplus_completeness_bins();
    }

    af::versa<double, af::c_grid<2> > iminus_completeness_bins() {
      return completeness_accumulator.iminus_completeness_bins();
    }

    af::versa<double, af::c_grid<2> > ieither_completeness_bins() {
      return completeness_accumulator.ieither_completeness_bins();
    }

    af::versa<double, af::c_grid<2> > iboth_completeness_bins() {
      return completeness_accumulator.iboth_completeness_bins();
    }

    af::shared<double> iplus_completeness() {
      return completeness_accumulator.iplus_completeness();
    }

    af::shared<double> iminus_completeness() {
      return completeness_accumulator.iminus_completeness();
    }

    af::shared<double> ieither_completeness() {
      return completeness_accumulator.ieither_completeness();
    }

    af::shared<double> iboth_completeness() {
      return completeness_accumulator.iboth_completeness();
    }

    af::versa<double, af::c_grid<2> > rcp_bins() {
      return rcp_scp_accumulator.rcp_bins();
    }

    af::versa<double, af::c_grid<2> > scp_bins() {
      return rcp_scp_accumulator.scp_bins();
    }

    af::shared<double> rcp() {
      return rcp_scp_accumulator.rcp();
    }

    af::shared<double> scp() {
      return rcp_scp_accumulator.scp();
    }

    af::shared<double> rd() {
      return rd_accumulator.rd();
    }

  private:
    Observations observations;

    accumulator::CompletenessAccumulator completeness_accumulator;
    accumulator::RcpScpAccumulator rcp_scp_accumulator;
    accumulator::RdAccumulator rd_accumulator;
  };

}}  // namespace dials::pychef

#endif
