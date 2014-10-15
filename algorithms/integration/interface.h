/*
 * interface.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H
#define DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <list>
#include <vector>
#include <dials/model/data/image.h>
#include <dials/model/data/shoebox.h>
#include <dials/array_family/reflection_table.h>
#include <dials/array_family/boost_python/flex_table_suite.h>


namespace dials { namespace algorithms {

  using model::Image;
  using model::Shoebox;

  /**
   * A class to compute integration jobs.
   */
  class IntegrationJobCalculator {
  public:

    /**
     * Compute the integration jobs
     * @param array_range The range of frames to process
     * @param block_size The number of frames in a job
     */
    IntegrationJobCalculator(
        vec2<int> array_range,
        double block_size) {
      int frame0 = array_range[0];
      int frame1 = array_range[1];
      DIALS_ASSERT(frame1 > frame0);
      int nframes = frame1 - frame0;
      DIALS_ASSERT(nframes > 0);
      if (block_size > nframes) {
        block_size = nframes;
      }
      DIALS_ASSERT(block_size > 0);
      if (block_size == 1) {
        for (int f = frame0; f < frame1; ++f) {
          jobs_.push_back(tiny<int,2>(f, f+1));
        }
      } else {
        int nblocks = (int)std::ceil(2.0 * nframes / block_size);
        DIALS_ASSERT(nblocks > 0 && nblocks <= nframes);
        int half_block_size = (int)std::ceil((double)nframes / (double)nblocks);
        af::shared<int> indices;
        indices.push_back(frame0);
        for (int i = 0; i < nblocks; ++i) {
          int frame = frame0 + (i + 1) * half_block_size;
          if (frame > frame1) {
            frame = frame1;
          }
          indices.push_back(frame);
          if (frame == frame1) {
            break;
          }
        }
        DIALS_ASSERT(indices.front() == frame0);
        DIALS_ASSERT(indices.back() == frame1);
        DIALS_ASSERT(indices.size() > 2);
        for (std::size_t i = 0; i < indices.size() - 2; ++i) {
          int i1 = indices[i];
          int i2 = indices[i+2];
          DIALS_ASSERT(i2 > i1);
          jobs_.push_back(tiny<int,2>(i1, i2));
        }
        DIALS_ASSERT(jobs_.size() > 0);
      }
    }

    /**
     * @returns The list of jobs.
     */
    af::shared< tiny<int,2> > jobs() const {
      return af::shared< tiny<int,2> >(&jobs_[0], &jobs_[0] + jobs_.size());
    }

  private:
    std::vector< tiny<int,2> > jobs_;
  };


  /**
   * A class to managing spliting and mergin data
   */
  class IntegrationManagerExecutor {
  public:

    IntegrationManagerExecutor(
        const IntegrationJobCalculator &jobcalculator,
        af::reflection_table reflections)
          : jobs_(jobcalculator.jobs()),
            data_(reflections) {
      DIALS_ASSERT(jobs_.size() > 0);

      // Set all the finished flags to false
      finished_.assign(jobs_.size(), false);

      // Generate indices of reflections to be integrated, used as reference
      // spots or passed just as data for each data block. If the reflection is
      // not to be integrated, it is added to each block which it overlaps. If
      // the reflection is a reference spot, it is added to each block in which
      // it is fully recorded. If the spot is to be integrated, it is added to
      // the block in which it is closest to the centre. If the reflection is
      // larger than block_size / 2, then it is not fully recorded in any block
      // and is unprocessed.
      typedef std::pair<std::size_t, bool> job_type;
      typedef std::vector<job_type> job_list_type;

      // Get some reflection data
      af::const_ref<int6> bbox = data_["bbox"];
      af::ref<std::size_t> flags = data_["flags"];

      // Check blocks are valid (can be overlapping)
      DIALS_ASSERT(jobs_[0][1] > jobs_[0][0]);
      for (std::size_t i = 1; i < jobs_.size(); ++i) {
        DIALS_ASSERT(jobs_[i][1] > jobs_[i][0]);
        DIALS_ASSERT(jobs_[i][0] > jobs_[i-1][0]);
        DIALS_ASSERT(jobs_[i][1] > jobs_[i-1][1]);
        DIALS_ASSERT(jobs_[i][0] <= jobs_[i-1][1]);
      }

      // Check all the reflections are in range
      int frame0 = jobs_.front()[0];
      int frame1 = jobs_.back()[1];
      DIALS_ASSERT(frame1 > frame0);
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        DIALS_ASSERT(bbox[i][1] > bbox[i][0]);
        DIALS_ASSERT(bbox[i][3] > bbox[i][2]);
        DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
        DIALS_ASSERT(bbox[i][4] >= frame0);
        DIALS_ASSERT(bbox[i][5] <= frame1);
      }

      // Create the lookups
      std::vector<std::size_t> lookup0(frame1 - frame0);
      std::vector<std::size_t> lookup1(frame1 - frame0);
      int frame = frame0;
      for (std::size_t i = 0; i < jobs_.size(); ++i) {
        tiny<int,2> b = jobs_[i];
        DIALS_ASSERT(frame >= b[0]);
        for (; frame < b[1]; ++frame) {
          lookup0[frame-frame0] = i;
        }
      }
      DIALS_ASSERT(frame == frame1);
      for (std::size_t i = 0; i < jobs_.size(); ++i) {
        std::size_t j = jobs_.size() - i - 1;
        tiny<int,2> b = jobs_[j];
        DIALS_ASSERT(frame <= b[1]);
        for (; frame > b[0]; --frame) {
          lookup1[frame-frame0-1] = j;
        }
      }
      DIALS_ASSERT(frame == frame0);

      // Check the lookups
      for (std::size_t i = 1; i < lookup0.size(); ++i) {
        DIALS_ASSERT(lookup0[i] >= lookup0[i-1]);
        DIALS_ASSERT(lookup1[i] >= lookup1[i-1]);
      }

      // Get which reflections to process in which job and task
      std::vector<job_list_type> indices(jobs_.size());
      for (std::size_t index = 0; index < bbox.size(); ++index) {
        int z0 = bbox[index][4];
        int z1 = bbox[index][5];
        std::size_t &f = flags[index];
        if (!(f & af::DontIntegrate)) {
          std::size_t j0 = lookup0[z0-frame0];
          std::size_t j1 = lookup1[z1-frame0-1];
          DIALS_ASSERT(j0 < jobs_.size());
          DIALS_ASSERT(j1 < jobs_.size());
          DIALS_ASSERT(j1 >= j0);
          DIALS_ASSERT(z0 >= jobs_[j0][0]);
          DIALS_ASSERT(z1 <= jobs_[j1][1]);
          std::size_t jmin = 0;
          double dmin = 0;
          bool inside = false;
          for (std::size_t j = j0; j <= j1; ++j) {
            int jz0 = jobs_[j][0];
            int jz1 = jobs_[j][1];
            if (z0 >= jz0 && z1 <= jz1) {
              if (f & af::ReferenceSpot) {
                indices[j].push_back(job_type(index, false));
              }
              double zc = (z1 + z0) / 2.0;
              double jc = (jz1 + jz0) / 2.0;
              double d = std::abs(zc - jc);
              if (!inside || d < dmin) {
                jmin = j;
                dmin = d;
                inside = true;
              }
            }
          }
          int jz0 = jobs_[jmin][0];
          int jz1 = jobs_[jmin][1];
          DIALS_ASSERT(inside == true);
          DIALS_ASSERT(z0 >= jz0 && z1 <= jz1);
          if (f & af::ReferenceSpot) {
            DIALS_ASSERT(indices[jmin].size() > 0);
            DIALS_ASSERT(indices[jmin].back().first == index);
            indices[jmin].back().second = true;
          } else {
            indices[jmin].push_back(job_type(index, true));
          }
        } else {
          ignored_.push_back(index);
        }
      }

      // Compute number of reflections in each task
      std::vector<std::size_t> num(jobs_.size(), 0);
      for (std::size_t i = 0; i < num.size(); ++i) {
        num[i] = indices[i].size();
      }

      // Compute offsests
      offset_.push_back(0);
      std::partial_sum(num.begin(), num.end(), std::back_inserter(offset_));

      // Compute indices
      indices_.resize(offset_.back());
      mask_.resize(offset_.back());
      std::size_t k = 0;
      for (std::size_t i = 0; i < indices.size(); ++i) {
        const job_list_type& ind = indices[i];
        for (std::size_t j = 0; j < ind.size(); ++j, ++k) {
          DIALS_ASSERT(k < indices_.size());
          indices_[k] = ind[j].first;
          mask_[k] = ind[j].second;
          DIALS_ASSERT(indices_[k] < bbox.size());
        }
      }
      DIALS_ASSERT(k == indices_.size());
    }

    /**
     * @returns The result data
     */
    af::reflection_table data() {
      DIALS_ASSERT(finished());
      return data_;
    }

    /**
     * @returns Is the process finished
     */
    bool finished() const {
      return finished_.all_eq(true);
    }

    /**
     * @returns The number of tasks
     */
    std::size_t size() const {
      return finished_.size();
    }

    /**
     * @returns The block indices
     */
    vec2<int> job(std::size_t index) const {
      DIALS_ASSERT(index < jobs_.size());
      return jobs_[index];
    }

    /**
     * @returns The list of reflections to not process.
     */
    af::shared<std::size_t> ignored() const {
      return af::shared<std::size_t>(&ignored_[0], &ignored_[0]+ignored_.size());
    }

    /**
     * @returns The reflections for a particular block.
     */
    af::reflection_table split(std::size_t index) {

      using namespace af::boost_python::flex_table_suite;

      // Check the input
      DIALS_ASSERT(index < finished_.size());
      af::const_ref<std::size_t> ind = indices(index);
      af::const_ref<bool> msk = mask(index);
      DIALS_ASSERT(ind.size() == msk.size());

      // Extract the reflection table
      af::reflection_table result = select_rows_index(data_, ind);

      // Extract the flags and set those reflections that are not to be
      // processed.
      af::ref<std::size_t> flags = result["flags"];
      for (std::size_t i = 0; i < flags.size(); ++i) {
        if (msk[i] == false) {
          flags[i] |= af::DontIntegrate;
        }
      }

      // Return the reflections
      return result;
    }

    /**
     * Accumulate the results.
     */
    void accumulate(std::size_t index, af::reflection_table result) {

      using namespace af::boost_python::flex_table_suite;

      // Check the input
      DIALS_ASSERT(index < finished_.size());
      DIALS_ASSERT(finished_[index] == false);
      af::const_ref<std::size_t> ind = indices(index);
      af::const_ref<bool> msk = mask(index);
      DIALS_ASSERT(ind.size() == msk.size());
      DIALS_ASSERT(ind.size() == result.size());

      // Set the result
      set_selected_rows_index_mask(data_, ind, msk, result);

      // Set finished flag
      finished_[index] = true;
    }

  private:

    /**
     * Get the indices for each job
     */
    af::const_ref<std::size_t> indices(std::size_t index) const {
      DIALS_ASSERT(index < offset_.size()-1);
      std::size_t i0 = offset_[index];
      std::size_t i1 = offset_[index+1];
      DIALS_ASSERT(i1 >= i0);
      std::size_t off = i0;
      std::size_t num = i1 - i0;
      DIALS_ASSERT(off + num <= indices_.size());
      return af::const_ref<std::size_t> (&indices_[off], num);
    }

    /**
     * Get the mask for each job
     */
    af::const_ref<bool> mask(std::size_t index) const {
      DIALS_ASSERT(index < offset_.size()-1);
      std::size_t i0 = offset_[index];
      std::size_t i1 = offset_[index+1];
      DIALS_ASSERT(i1 >= i0);
      std::size_t off = i0;
      std::size_t num = i1 - i0;
      DIALS_ASSERT(off + num <= mask_.size());
      return af::const_ref<bool> (&mask_[off], num);
    }

    af::shared< tiny<int,2> > jobs_;
    af::shared<bool> finished_;
    af::reflection_table data_;
    af::shared<std::size_t> offset_;
    af::shared<std::size_t> indices_;
    af::shared<bool> mask_;
   af::shared<std::size_t> ignored_;
  };


  /**
   * A class to produce a summary of integation stats
   */
  class Summary {
  public:

    /**
     * A class to produce a summary of per image stats
     */
    class ImageSummary {
    public:

      ImageSummary(af::reflection_table data,
                   int2 image_range) {

        // Check the input
        DIALS_ASSERT(image_range[1] > image_range[0]);
        DIALS_ASSERT(data.size() > 0);
        DIALS_ASSERT(data.is_consistent());
        DIALS_ASSERT(data.contains("bbox"));
        DIALS_ASSERT(data.contains("partiality"));
        DIALS_ASSERT(data.contains("intensity.sum.value"));
        DIALS_ASSERT(data.contains("intensity.sum.variance"));
        DIALS_ASSERT(data.contains("flags"));

        // Get some data
        af::const_ref<int6> bbox = data["bbox"];
        af::const_ref<double> partiality = data["partiality"];
        af::const_ref<double> isum = data["intensity.sum.value"];
        af::const_ref<double> vsum = data["intensity.sum.variance"];
        af::const_ref<std::size_t> flags = data["flags"];

        // Maybe get profile fitting results
        af::const_ref<double> iprf(0, 0);
        af::const_ref<double> vprf(0, 0);
        bool prf = false;
        if (data.contains("intensity.prf.value") &&
            data.contains("intensity.prf.variance")) {
          iprf = data["intensity.prf.value"];
          vprf = data["intensity.prf.variance"];
          prf = true;
        }

        // Allocate the arrays to the number of images
        std::size_t nimages = image_range[1] - image_range[0];
        full_.assign(nimages, 0);
        part_.assign(nimages, 0);
        sum_cnt_.assign(nimages, 0);
        prf_cnt_.assign(nimages, 0);
        prf_ios_.assign(nimages, 0);
        sum_ios_.assign(nimages, 0);

        // Fully recorded threshold
        const double EPS = 1e-7;
        const double full_value = 0.997300203937 - EPS;

        // Compute the statistics
        for (std::size_t i = 0; i < data.size(); ++i) {
          DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
          DIALS_ASSERT(bbox[i][4] >= image_range[0]);
          DIALS_ASSERT(bbox[i][5] <= image_range[1]);
          if (!(flags[i] & af::DontIntegrate)) {
            DIALS_ASSERT(partiality[i] >= 0 && partiality[i] <= 1);
            for (std::size_t j = bbox[i][4]; j < bbox[i][5]; ++j) {
              std::size_t k = j - image_range[0];
              if (partiality[i] >= full_value) {
                full_[k]++;
              } else {
                part_[k]++;
              }
            }
          }
          if (!(flags[i] & af::DontIntegrate) &&
               flags[i] & af::IntegratedSum &&
               vsum[i] > 0) {
            double ios = isum[i] / std::sqrt(vsum[i]);
            for (std::size_t j = bbox[i][4]; j < bbox[i][5]; ++j) {
              std::size_t k = j - image_range[0];
              sum_ios_[k] += ios;
              sum_cnt_[k]++;
            }
          }
          if (!(flags[i] & af::DontIntegrate) &&
               flags[i] & af::IntegratedPrf &&
               prf && vprf[i] > 0) {
            double ios = iprf[i] / std::sqrt(vprf[i]);
            for (std::size_t j = bbox[i][4]; j < bbox[i][5]; ++j) {
              std::size_t k = j - image_range[0];
              prf_ios_[k] += ios;
              prf_cnt_[k]++;
            }
          }
        }

        // Compute mean of ios for sum and profile fitting
        for (std::size_t i = 0; i < nimages; ++i) {
          if (sum_cnt_[i] > 0) sum_ios_[i] /= sum_cnt_[i];
          if (prf_cnt_[i] > 0) prf_ios_[i] /= prf_cnt_[i];
        }
      }

      af::shared<std::size_t> full() const {
        return full_;
      }

      af::shared<std::size_t> part() const {
        return part_;
      }

      af::shared<std::size_t> sum_count() const {
        return sum_cnt_;
      }

      af::shared<std::size_t> prf_count() const {
        return prf_cnt_;
      }

      af::shared<double> prf_ios() const {
        return prf_ios_;
      }

      af::shared<double> sum_ios() const {
        return sum_ios_;
      }

    private:
      af::shared<std::size_t> full_;
      af::shared<std::size_t> part_;
      af::shared<std::size_t> sum_cnt_;
      af::shared<std::size_t> prf_cnt_;
      af::shared<double> prf_ios_;
      af::shared<double> sum_ios_;
    };

    /**
     * Class to produce statistics at resolution bins
     */
    class ResolutionSummary {
    public:

      ResolutionSummary(af::reflection_table data,
                        std::size_t nbins) {

        // Check the input
        DIALS_ASSERT(data.size() > 0);
        DIALS_ASSERT(data.is_consistent());
        DIALS_ASSERT(data.contains("d"));
        DIALS_ASSERT(data.contains("partiality"));
        DIALS_ASSERT(data.contains("intensity.sum.value"));
        DIALS_ASSERT(data.contains("intensity.sum.variance"));
        DIALS_ASSERT(data.contains("flags"));

        // Get some data
        af::const_ref<double> d = data["d"];
        af::const_ref<double> partiality = data["partiality"];
        af::const_ref<double> isum = data["intensity.sum.value"];
        af::const_ref<double> vsum = data["intensity.sum.variance"];
        af::const_ref<std::size_t> flags = data["flags"];

        // Maybe get profile fitting results
        af::const_ref<double> iprf(0, 0);
        af::const_ref<double> vprf(0, 0);
        bool prf = false;
        if (data.contains("intensity.prf.value") &&
            data.contains("intensity.prf.variance")) {
          iprf = data["intensity.prf.value"];
          vprf = data["intensity.prf.variance"];
          prf = true;
        }

        // Compute the resolution bins
        double d_min = af::min(d);
        double d_max = af::max(d);
        DIALS_ASSERT(d_max > d_min);
        bin_range_ = (d_max - d_min) / nbins;
        for (std::size_t i = 0; i <= nbins; ++i) {
          bins_.push_back(d_min + bin_range_ * i);
        }

        // Allocate other arrays
        sum_ios_.assign(nbins, 0);
        sum_ios_full_.assign(nbins, 0);
        sum_ios_part_.assign(nbins, 0);
        prf_ios_.assign(nbins, 0);
        prf_ios_full_.assign(nbins, 0);
        prf_ios_part_.assign(nbins, 0);
        sum_cnt_.assign(nbins, 0);
        sum_cnt_full_.assign(nbins, 0);
        sum_cnt_part_.assign(nbins, 0);
        prf_cnt_.assign(nbins, 0);
        prf_cnt_full_.assign(nbins, 0);
        prf_cnt_part_.assign(nbins, 0);

        // Fully recorded threshold
        const double EPS = 1e-7;
        const double full_value = 0.997300203937 - EPS;

        // Compute the stats
        for (std::size_t i = 0; i < data.size(); ++i) {
          DIALS_ASSERT(d[i] >= d_min && d[i] <= d_max);
          int j = std::floor((d[i] - d_min) / bin_range_);
          DIALS_ASSERT(j >= 0 && j <= nbins);
          if (j == nbins) j = nbins-1;
          if (!(flags[i] & af::DontIntegrate)) {
            if (flags[i] & af::IntegratedSum &&
                vsum[i] > 0) {
              double ios = isum[i] / std::sqrt(vsum[i]);
              if (partiality[i] >= full_value) {
                sum_ios_full_[j] += ios;
                sum_cnt_full_[j]++;
              } else {
                sum_ios_part_[j] += ios;
                sum_cnt_part_[j]++;
              }
              sum_ios_[j] += ios;
              sum_cnt_[j]++;
            }
            if (flags[i] & af::IntegratedPrf && prf && vprf[i] > 0) {
              double ios = iprf[i] / std::sqrt(vprf[i]);
              if (partiality[i] >= full_value) {
                prf_ios_full_[j] += ios;
                prf_cnt_full_[j]++;
              } else {
                prf_ios_part_[j] += ios;
                prf_cnt_part_[j]++;
              }
              prf_ios_[j] += ios;
              prf_cnt_[j]++;
            }
          }
        }

        // Take the averages
        for (std::size_t i = 0; i < nbins; ++i) {
          if (sum_cnt_[i] > 0) sum_ios_[i] /= sum_cnt_[i];
          if (sum_cnt_full_[i] > 0) sum_ios_full_[i] /= sum_cnt_full_[i];
          if (sum_cnt_part_[i] > 0) sum_ios_part_[i] /= sum_cnt_part_[i];
          if (prf_cnt_[i] > 0) prf_ios_[i] /= prf_cnt_[i];
          if (prf_cnt_full_[i] > 0) prf_ios_full_[i] /= prf_cnt_full_[i];
          if (prf_cnt_part_[i] > 0) prf_ios_part_[i] /= prf_cnt_part_[i];
        }
      }

      std::size_t nbins() const {
        return bins_.size();
      }

      double bin_range() const {
        return bin_range_;
      }

      af::shared<double> bins() const {
        return bins_;
      }

      af::shared<double> sum_ios() const {
        return sum_ios_;
      }

      af::shared<double> sum_ios_full() const {
        return sum_ios_full_;
      }

      af::shared<double> sum_ios_part() const {
        return sum_ios_part_;
      }

      af::shared<double> prf_ios() const {
        return prf_ios_;
      }

      af::shared<double> prf_ios_full() const {
        return prf_ios_full_;
      }

      af::shared<double> prf_ios_part() const {
        return prf_ios_part_;
      }

    private:

      double bin_range_;
      af::shared<double> bins_;
      af::shared<double> sum_ios_;
      af::shared<double> sum_ios_full_;
      af::shared<double> sum_ios_part_;
      af::shared<double> prf_ios_;
      af::shared<double> prf_ios_full_;
      af::shared<double> prf_ios_part_;
      af::shared<std::size_t> sum_cnt_;
      af::shared<std::size_t> sum_cnt_full_;
      af::shared<std::size_t> sum_cnt_part_;
      af::shared<std::size_t> prf_cnt_;
      af::shared<std::size_t> prf_cnt_full_;
      af::shared<std::size_t> prf_cnt_part_;
    };

    /**
     * A class to produce summary stats for whole dataset
     */
    class WholeSummary {
    public:

      WholeSummary(af::reflection_table data)
        : sum_ios_(0),
          sum_ios_full_(0),
          sum_ios_part_(0),
          prf_ios_(0),
          prf_ios_full_(0),
          prf_ios_part_(0),
          sum_cnt_(0),
          sum_cnt_full_(0),
          sum_cnt_part_(0),
          prf_cnt_(0),
          prf_cnt_full_(0),
          prf_cnt_part_(0) {

        // Check the input
        DIALS_ASSERT(data.size() > 0);
        DIALS_ASSERT(data.is_consistent());
        DIALS_ASSERT(data.contains("partiality"));
        DIALS_ASSERT(data.contains("intensity.sum.value"));
        DIALS_ASSERT(data.contains("intensity.sum.variance"));
        DIALS_ASSERT(data.contains("flags"));

        // Get some data
        af::const_ref<double> partiality = data["partiality"];
        af::const_ref<double> isum = data["intensity.sum.value"];
        af::const_ref<double> vsum = data["intensity.sum.variance"];
        af::const_ref<std::size_t> flags = data["flags"];

        // Maybe get profile fitting results
        af::const_ref<double> iprf(0, 0);
        af::const_ref<double> vprf(0, 0);
        bool prf = false;
        if (data.contains("intensity.prf.value") &&
            data.contains("intensity.prf.variance")) {
          iprf = data["intensity.prf.value"];
          vprf = data["intensity.prf.variance"];
          prf = true;
        }

        // Fully recorded threshold
        const double EPS = 1e-7;
        const double full_value = 0.997300203937 - EPS;

        // Compute the stats
        for (std::size_t i = 0; i < data.size(); ++i) {
          if (!(flags[i] & af::DontIntegrate)) {
            if (flags[i] & af::IntegratedSum && vsum[i] > 0) {
              double ios = isum[i] / std::sqrt(vsum[i]);
              if (partiality[i] >= full_value) {
                sum_ios_full_ += ios;
                sum_cnt_full_++;
              } else {
                sum_ios_part_ += ios;
                sum_cnt_part_++;
              }
              sum_ios_ += ios;
              sum_cnt_++;
            }
            if (flags[i] & af::IntegratedPrf && prf && vprf[i] > 0) {
              double ios = iprf[i] / std::sqrt(vprf[i]);
              if (partiality[i] >= full_value) {
                prf_ios_full_ += ios;
                prf_cnt_full_++;
              } else {
                prf_ios_part_ += ios;
                prf_cnt_part_++;
              }
              prf_ios_ += ios;
              prf_cnt_++;
            }
          }
        }

        // Take the averages
        if (sum_cnt_ > 0) sum_ios_ /= sum_cnt_;
        if (sum_cnt_full_ > 0) sum_ios_full_ /= sum_cnt_full_;
        if (sum_cnt_part_ > 0) sum_ios_part_ /= sum_cnt_part_;
        if (prf_cnt_ > 0) prf_ios_ /= prf_cnt_;
        if (prf_cnt_full_ > 0) prf_ios_full_ /= prf_cnt_full_;
        if (prf_cnt_part_ > 0) prf_ios_part_ /= prf_cnt_part_;
      }

      double sum_ios() const {
        return sum_ios_;
      }

      double sum_ios_full() const {
        return sum_ios_full_;
      }

      double sum_ios_part() const {
        return sum_ios_part_;
      }

      double prf_ios() const {
        return prf_ios_;
      }

      double prf_ios_full() const {
        return prf_ios_full_;
      }

      double prf_ios_part() const {
        return prf_ios_part_;
      }

    private:
      double sum_ios_;
      double sum_ios_full_;
      double sum_ios_part_;
      double prf_ios_;
      double prf_ios_full_;
      double prf_ios_part_;
      std::size_t sum_cnt_;
      std::size_t sum_cnt_full_;
      std::size_t sum_cnt_part_;
      std::size_t prf_cnt_;
      std::size_t prf_cnt_full_;
      std::size_t prf_cnt_part_;
    };

    Summary(af::reflection_table data,
            int2 image_range,
            std::size_t nresbins)
      : img_summary_(data, image_range),
        res_summary_(data, nresbins),
        who_summary_(data) {
    }

    ImageSummary image_summary() {
      return img_summary_;
    }

    ResolutionSummary resolution_summary() {
      return res_summary_;
    }

    WholeSummary whole_summary() {
      return who_summary_;
    }

  private:

    ImageSummary img_summary_;
    ResolutionSummary res_summary_;
    WholeSummary who_summary_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H
