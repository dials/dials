
#ifndef DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LEARNER_2D_H
#define DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LEARNER_2D_H

#include <scitbx/vec3.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/algorithms/integration/profile/grid_sampler_2d.h>
#include <dials/model/data/shoebox.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using scitbx::vec3;

  class ReferenceLearner2D {
  public:

    typedef af::versa< double, af::c_grid<2> > double2d;
    typedef af::versa< int, af::c_grid<2> > int2d;

    ReferenceLearner2D(const GridSampler2D &sampler,
                       const af::const_ref<std::size_t> &xsize,
                       const af::const_ref<std::size_t> &ysize,
                       std::size_t min_average)
      : sampler_(sampler),
        counts_(sampler.size(), 0),
        profiles_(sampler.size()),
        masks_(sampler.size()),
        min_average_(min_average) {
      DIALS_ASSERT(xsize.size() == ysize.size());
      DIALS_ASSERT(xsize.size() == sampler_.size());
      for (std::size_t i = 0; i < xsize.size(); ++i) {
        profiles_[i] = double2d(af::c_grid<2>(ysize[i], xsize[i]),0);
        masks_[i] = int2d(af::c_grid<2>(ysize[i], xsize[i]),0);
      }
    }

    void add(
        const af::const_ref< double, af::c_grid<2> > &data,
        const af::const_ref< double, af::c_grid<2> > &background,
        const af::const_ref< int, af::c_grid<2> > &mask,
        const vec3<double> &coord) {
      DIALS_ASSERT(data.accessor().all_eq(background.accessor()));
      DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
      add_internal(data, background, mask,
          sampler_.nearest(double2(coord[0], coord[1])));
    }

    void finish() {
      for (std::size_t i = 0; i < counts_.size(); ++i) {
        if (counts_[i] < min_average_) {
          average_neighbours(i);
          normalize(profiles_[i].ref(), masks_[i].const_ref());
        } else {
          normalize(profiles_[i].ref(), masks_[i].const_ref());
        }
      }
    }

    af::shared<std::size_t> counts() const {
      return af::shared<std::size_t>(counts_.begin(), counts_.end());
    }

    double2d profile(std::size_t index) const {
      DIALS_ASSERT(index < profiles_.size());
      return profiles_[index];
    }

    int2d mask(std::size_t index) const {
      DIALS_ASSERT(index < masks_.size());
      return masks_[index];
    }

    std::size_t size() const {
      return sampler_.size();
    }

    double2 coord(std::size_t index) const {
      return sampler_[index];
    }

  private:

    void add_internal(
        const af::const_ref< double, af::c_grid<2> > &data,
        const af::const_ref< double, af::c_grid<2> > &background,
        const af::const_ref< int, af::c_grid<2> > &mask,
        std::size_t index) {
      DIALS_ASSERT(index < counts_.size());
      double2d& p = profiles_[index];
      int2d& m = masks_[index];
      DIALS_ASSERT(p.accessor().all_eq(data.accessor()));
      for (std::size_t i = 0; i < p.size(); ++i) {
        p[i] += data[i] - background[i];
        m[i] |= mask[i];
      }
      counts_[index]++;
    }

    void normalize(
        af::ref< double, af::c_grid<2> > data,
        const af::const_ref< int, af::c_grid<2> > &mask) {
      double tot = 0;
      for (std::size_t i = 0; i < data.size(); ++i) {
        if (mask[i] | Foreground) {
          tot += data[i];
        }
      }
      DIALS_ASSERT(tot > 0);
      for (std::size_t i = 0; i < data.size(); ++i) {
        data[i] /= tot;
      }
    }

    void average_neighbours(std::size_t index) {
      af::shared<std::size_t> neighbours = sampler_.neighbours(index);
      double2d& p = profiles_[index];
      int2d& m = masks_[index];
      std::size_t ys1 = p.accessor()[0] / 2;
      std::size_t xs1 = p.accessor()[1] / 2;
      for (std::size_t i = 0; i < neighbours.size(); ++i) {
        std::size_t j = neighbours[i];
        DIALS_ASSERT(j < counts_.size());
        double2d data = profiles_[j];
        int2d mask = masks_[j];
        std::size_t ys2 = data.accessor()[0] / 2;
        std::size_t xs2 = data.accessor()[1] / 2;
        int y0 = ys1 - ys2;
        int x0 = xs1 - xs2;
        DIALS_ASSERT(p.accessor().all_eq(data.accessor()));
        for (int yy = 0; yy < ys2; ++yy) {
          for (int xx = 0; xx < xs2; ++xx) {
            int yyy = yy + y0;
            int xxx = xx + x0;
            if (yyy >= 0 && xxx >= 0 && yyy < ys1 && xxx < xs1) {
              p(yyy,xxx) += data(yy,xx);
              m(yyy,xxx) |= mask(yy,xx);
            }
          }
        }
        counts_[i] += counts_[j];
      }
    }

    GridSampler2D sampler_;
    af::shared<std::size_t> counts_;
    af::shared<double2d> profiles_;
    af::shared<int2d> masks_;
    std::size_t min_average_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_PROFILE_REFERENCE_LEARNER_2D_H
