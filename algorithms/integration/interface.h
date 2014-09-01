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
#include <list>
#include <vector>
#include <boost/function.hpp>
#include <cctbx/miller/index_generator.h>
#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/space_group.h>
#include <dials/model/data/image.h>
#include <dials/model/data/shoebox.h>
#include <dials/array_family/reflection_table.h>
#include <dials/array_family/boost_python/flex_table_suite.h>
#include <dials/algorithms/reflection_basis/coordinate_system.h>


namespace dials { namespace algorithms {

  using model::Image;
  using model::Shoebox;

  class IntegrationTask3DExecutor {
  public:

    typedef af::reflection_table rtable;
    typedef boost::function<rtable (rtable)> callback_type;

    IntegrationTask3DExecutor(
            af::reflection_table data,
            const af::const_ref< tiny<int,2> > &jobs,
            std::size_t npanels,
            callback_type callback) 
        : data_(data),
          jobs_(jobs.begin(), jobs.end()),
          npanels_(npanels),
          process_(callback) {
    
      // Check the jobs are valid. Jobs must be ordered in increasing frame
      // number and can overlap but 1 jobs must not be fully contained in
      // another. There must also be no gaps.  
      DIALS_ASSERT(jobs.size() > 0);
      for (std::size_t i = 0; i < jobs.size()-1; ++i) {
        DIALS_ASSERT(jobs[i][0] <  jobs[i][1]);
        DIALS_ASSERT(jobs[i][0] <  jobs[i+1][0]);
        DIALS_ASSERT(jobs[i][1] <  jobs[i+1][1]);
        DIALS_ASSERT(jobs[i][1] >= jobs[i+1][0]);
      }

      // Get the range of frames
      frame0_ = jobs_.front()[0];
      frame1_ = jobs_.back()[1];
      frame_ = frame0_;
      DIALS_ASSERT(frame0_ < frame1_);
      nframes_ = frame1_ - frame0_;

      // Set the current active job
      active_job_ = 0;
      
      // Check the reflection table contains a few required fields
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.contains("job_id"));
      DIALS_ASSERT(data.contains("panel"));
      DIALS_ASSERT(data.contains("bbox"));

      // Get some data from the reflection table
      af::const_ref< std::size_t > job_id = data["job_id"];
      af::const_ref< std::size_t > panel = data["panel"];
      af::const_ref< int6 > bbox = data["bbox"];

      // Check the jobs ids are valid
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        std::size_t j = job_id[i];
        int z0 = bbox[i][4];
        int z1 = bbox[i][5];
        DIALS_ASSERT(z1 > z0);
        DIALS_ASSERT(j < jobs.size());
        int j0 = jobs[j][0];
        int j1 = jobs[j][1];
        DIALS_ASSERT(j1 > j0);
        DIALS_ASSERT(z0 >= j0 && z1 <= j1);
      }

      // Initialise the shoeboxes
      shoebox_ = data_["shoebox"];
      for (std::size_t i = 0; i < shoebox_.size(); ++i) {
        shoebox_[i] = Shoebox<>(panel[i], bbox[i]);
      }

      // Initialise the offsets and indices for each frame/panel
      initialise_indices();
      initialise_job_indices();
      malloc(0);
    }

    int frame0() const {
      return frame0_;
    }

    int frame1() const {
      return frame1_;
    }

    int frame() const {
      return frame_;
    }

    std::size_t nframes() const {
      return nframes_;
    }
    
    std::size_t job() const {
      return active_job_;
    }

    void next(const Image &image) {
      DIALS_ASSERT(!finished());
      next_image(image);
      if (ready()) {
        merge(process_(split()));
        next_job();
      }
    }

    af::reflection_table data() {
      DIALS_ASSERT(finished());
      data_.erase("shoebox");
      return data_;
    }

    bool finished() const {
      return frame_ == frame1_; 
    }

  private:

    void initialise_indices() {
      typedef std::list<std::size_t> list_type;
      typedef list_type::iterator iterator;

      // Temporary index arrays
      std::size_t num = nframes_ * npanels_;
      std::size_t none = shoebox_.size();
      std::vector<iterator>  temp_off(num + 1);
      std::list<std::size_t> temp_ind;
      temp_off[0] = temp_ind.begin();
      for (std::size_t i = 1; i < temp_off.size(); ++i) {
        temp_off[i] = temp_ind.insert(temp_ind.end(), none);
      }

      // Insert all the partial indices
      for (std::size_t i = 0; i < shoebox_.size(); ++i) {
        std::size_t p = shoebox_[i].panel;
        int6 &b = shoebox_[i].bbox;
        for (int z = b[4]; z < b[5]; ++z) {
          int j = p + (z - frame0_)*npanels_+1;
          temp_ind.insert(temp_off[j], i);
        }
      }

      // Copy the temporary indices
      offset_.resize(temp_off.size());
      indices_.resize(temp_ind.size());
      offset_[0] = 0;
      for (std::size_t i = 0; i < temp_off.size()-1; ++i) {
        std::size_t num = std::distance(temp_off[i], temp_off[i+1]);
        offset_[i+1] = offset_[i] + num;
      }
      DIALS_ASSERT(offset_.back() == indices_.size());
      std::copy(temp_ind.begin(), temp_ind.end(), indices_.begin());
    }

    void initialise_job_indices() {
      
      af::const_ref<std::size_t> job_id = data_["job_id"];
     
      std::vector<std::size_t> num(jobs_.size(), 0);
      std::vector<std::size_t> count(jobs_.size(), 0);
      for (std::size_t i = 0; i < job_id.size(); ++i) {
        num[job_id[i]]++;
      }
      job_offset_.resize(jobs_.size() + 1);
      job_indices_.resize(job_id.size());
      job_offset_[0] = 0;
      for (std::size_t i = 1; i < job_offset_.size(); ++i) {
        job_offset_[i] = job_offset_[i-1] + num[i-1];
      }
      for (std::size_t i = 0; i < job_id.size(); ++i) {
        std::size_t j = job_id[i];
        std::size_t k = job_offset_[j] + count[j];
        job_indices_[k] = i;
        count[j]++;
      }
      for (std::size_t i = 0; i < num.size(); ++i) {
        DIALS_ASSERT(count[i] == num[i]);
      }
    }

    void next_image(const Image &image) {
      DIALS_ASSERT(image.npanels() == npanels_);
      for (std::size_t p = 0; p < image.npanels(); ++p) {
        af::const_ref<std::size_t> ind = indices(frame_, p);
        af::const_ref< int, af::c_grid<2> > data = image.data(p);
        af::const_ref< bool, af::c_grid<2> > mask = image.mask(p);
        DIALS_ASSERT(data.accessor().all_eq(mask.accessor()));
        for (std::size_t i = 0; i < ind.size(); ++i) {
          DIALS_ASSERT(ind[i] < shoebox_.size());
          Shoebox<>& sbox = shoebox_[ind[i]];
          int6 b = sbox.bbox;
          DIALS_ASSERT(b[1] > b[0]);
          DIALS_ASSERT(b[3] > b[2]);
          DIALS_ASSERT(b[5] > b[4]);
          DIALS_ASSERT(frame_ >= b[4] && frame_ < b[5]);
          int x0 = b[0];
          int x1 = b[1];
          int y0 = b[2];
          int y1 = b[3];
          int z0 = b[4];
          std::size_t xs = x1 - x0;
          std::size_t ys = y1 - y0;
          std::size_t z = frame_ - z0;
          DIALS_ASSERT(x0 >= 0 && y0 >= 0);
          DIALS_ASSERT(y1 <= data.accessor()[0]);
          DIALS_ASSERT(x1 <= data.accessor()[1]);
          DIALS_ASSERT(sbox.is_consistent());
          for (std::size_t y = 0; y < ys; ++y) {
            for (std::size_t x = 0; x < xs; ++x) {
              sbox.data(z, y, x) = data(y+y0,x+x0);
              sbox.mask(z, y, x) = mask(y+y0,x+x0) ? Valid : 0;
            }
          }
        }
      }
      frame_++;
    }

    af::const_ref<std::size_t> indices(int frame, std::size_t panel) const {
      std::size_t j0 = panel+(frame-frame0_)*npanels_;
      DIALS_ASSERT(offset_.size() > 0);
      DIALS_ASSERT(j0 < offset_.size()-1);
      std::size_t i0 = offset_[j0];
      std::size_t i1 = offset_[j0+1];
      DIALS_ASSERT(i1 > i0);
      std::size_t off = i0;
      std::size_t num = i1 - off - 1;
      DIALS_ASSERT(off + num <= indices_.size());
      return af::const_ref<std::size_t>(&indices_[off], num);
    }

    bool ready() const {
      DIALS_ASSERT(active_job_ < jobs_.size());
      int j0 = jobs_[active_job_][0];
      int j1 = jobs_[active_job_][1];
      DIALS_ASSERT(frame_ >= j0 && frame_ <= j1);
      return frame_ == j1;
    }

    void next_job() {
      free(active_job_);
      active_job_++;
      malloc(active_job_);
    }

    void free(std::size_t job) {
      af::const_ref<std::size_t> ind = job_indices(job);
      for (std::size_t i = 0; i < ind.size(); ++i) {
        shoebox_[ind[i]].deallocate();
      }
    }

    void malloc(std::size_t job) {
      int z0 = jobs_[job][0];
      int z1 = jobs_[job][1];
      for (std::size_t j = job; j < jobs_.size(); ++j) {
        int z2 = jobs_[j][0];
        if (z2 < z1) {
          af::const_ref<std::size_t> ind = job_indices(j);
          for (std::size_t i = 0; i < ind.size(); ++i) {
            shoebox_[ind[i]].allocate();
          }
        }
      }
    }

    af::reflection_table split() const {
      using namespace dials::af::boost_python::flex_table_suite;
      return select_rows_index(data_, job_indices(active_job_));
    }

    void merge(af::reflection_table result) {
      using namespace dials::af::boost_python::flex_table_suite;
      set_selected_rows_index(data_, job_indices(active_job_), result);
    }

    af::const_ref<std::size_t> job_indices(std::size_t index) const {
      DIALS_ASSERT(index < job_offset_.size()-1); 
      std::size_t i0 = job_offset_[index];
      std::size_t i1 = job_offset_[index+1];
      DIALS_ASSERT(i1 > i0);
      std::size_t off = i0;
      std::size_t num = i1 - i0;
      DIALS_ASSERT(off + num <= indices_.size());
      return af::const_ref<std::size_t> (&job_indices_[off], num);
    }

    af::reflection_table data_;
    std::vector< tiny<int,2> > jobs_;
    int frame0_;
    int frame1_;
    int frame_;
    std::size_t nframes_;
    std::size_t npanels_;
    std::size_t active_job_;
    callback_type process_;
    af::shared< Shoebox<> > shoebox_;
    std::vector<std::size_t> offset_;
    std::vector<std::size_t> indices_;
    std::vector<std::size_t> job_offset_;
    std::vector<std::size_t> job_indices_;
  };

  /* /** */
  /*  * A class to managing spliting and mergin data */
  /*  *1/ */
  /* class IntegrationManagerData3D { */
  /* public: */

  /*   IntegrationManagerData3D( */
  /*       af::reflection_table reflections, */
  /*       vec2<double> oscillation, */
  /*       vec2<int> array_range, */
  /*       double block_size, */
  /*       std::size_t max_procs */
  /*       std::size_t max_tasks) */
  /*         : data_(reflections) { */

  /*     // Compute the blocks */
  /*     compute_blocks(oscillation, array_range, block_size); */
  /*     finished_.assign(blocks_.size(), false); */

  /*     // Get the bounding boxes and flags */
  /*     af::const_ref<int6> bbox = reflections["bbox"]; */
  /*     af::ref<std::size_t> flags = reflections["flags"]; */

  /*     // Generate indices of reflections to be integrated, used as reference */
  /*     // spots or passed just as data for each data block. If the reflection is */
  /*     // not to be integrated, it is added to each block which it overlaps. If */
  /*     // the reflection is a reference spot, it is added to each block in which */
  /*     // it is fully recorded. If the spot is to be integrated, it is added to */
  /*     // the block in which it is closest to the centre. If the reflection is */
  /*     // larger than block_size / 2, then it is not fully recorded in any block */
  /*     // and is unprocessed. */
  /*     to_process_.resize(blocks_.size()); */
  /*     to_include_.resize(blocks_.size()); */
  /*     for (std::size_t i = 0; i < bbox.size(); ++i) { */
  /*       int z0 = bbox[i][4]; */
  /*       int z1 = bbox[i][5]; */
  /*       std::size_t &f = flags[i]; */
  /*       if (!(f & af::DontIntegrate)) { */
  /*         std::size_t jmin = 0; */
  /*         double dmin = 0; */
  /*         for (std::size_t j = 0; j < blocks_.size(); ++j) { */
  /*           int bz0 = blocks_[j][0]; */
  /*           int bz1 = blocks_[j][1]; */
  /*           if (f & af::ReferenceSpot) { */
  /*             if (z0 >= bz0 && z1 <= bz1) { */
  /*               to_include_[j].push_back(i); */
  /*             } */
  /*           } */
  /*           double zc = (z1 + z0) / 2.0; */
  /*           double bc = (bz1 + bz0) / 2.0; */
  /*           double d = (zc - bc)*(zc - bc); */
  /*           if (j == 0 || d < dmin) { */
  /*             jmin = j; */
  /*             dmin = d; */
  /*           } */
  /*         } */
  /*         int bz0 = blocks_[jmin][0]; */
  /*         int bz1 = blocks_[jmin][1]; */
  /*         if (z0 >= bz0 && z1 <= bz1) { */
  /*           to_process_[jmin].push_back(i); */
  /*         } else { */
  /*           to_not_process_.push_back(i); */
  /*           f |= af::DontIntegrate; */
  /*         } */
  /*       } */
  /*       /1* if (f & af::DontIntegrate) { *1/ */
  /*       /1*   for (std::size_t j = 0; j < blocks_.size(); ++j) { *1/ */
  /*       /1*     int bz0 = blocks_[j][0]; *1/ */
  /*       /1*     int bz1 = blocks_[j][1]; *1/ */
  /*       /1*     if (!(z1 <= bz0 || z0 >= bz1)) { *1/ */
  /*       /1*       to_include_[j].push_back(i); *1/ */
  /*       /1*     } *1/ */
  /*       /1*   } *1/ */
  /*       /1* } *1/ */
  /*     } */
  /*   } */

  /*   /** */
  /*    * @returns The result data */
  /*    *1/ */
  /*   af::reflection_table data() { */
  /*     DIALS_ASSERT(finished()); */
  /*     return data_; */
  /*   } */

  /*   /** */
  /*    * @returns Is the process finished */
  /*    *1/ */
  /*   bool finished() const { */
  /*     return finished_.all_eq(true); */
  /*   } */

  /*   /** */
  /*    * @returns The number of tasks */
  /*    *1/ */
  /*   std::size_t size() const { */
  /*     return finished_.size(); */
  /*   } */

  /*   /** */
  /*    * @returns The block indices */
  /*    *1/ */
  /*   vec2<int> block(std::size_t index) const { */
  /*     DIALS_ASSERT(index < blocks_.size()); */
  /*     return blocks_[index]; */
  /*   } */

  /*   /** */
  /*    * @returns The number of jobs within the block. */
  /*    *1/ */
  /*   std::size_t njobs(std::size_t index) const { */
  /*     DIALS_ASSERT(index < njobs_.size()); */
  /*     return njobs_[index]; */
  /*   } */

  /*   /** */
  /*    * @returns The list of reflections to not process. */
  /*    *1/ */
  /*   af::shared<std::size_t> to_not_process() const { */
  /*     return af::shared<std::size_t>( */
  /*         &to_not_process_[0], */
  /*         &to_not_process_[0] + to_not_process_.size()); */
  /*   } */

  /*   /** */
  /*    * @returns The list of reflections to include. */
  /*    *1/ */
  /*   af::shared<std::size_t> to_include(std::size_t index) const { */
  /*     DIALS_ASSERT(index < blocks_.size()); */
  /*     return af::shared<std::size_t>( */
  /*         &to_include_[index][0], */
  /*         &to_include_[index][0] + to_include_[index].size()); */
  /*   } */

  /*   /** */
  /*    * @returns The list of reflections to process. */
  /*    *1/ */
  /*   af::shared<std::size_t> to_process(std::size_t index) const { */
  /*     DIALS_ASSERT(index < blocks_.size()); */
  /*     return af::shared<std::size_t>( */
  /*         &to_process_[index][0], */
  /*         &to_process_[index][0] + to_process_[index].size()); */
  /*   } */

  /*   /** */
  /*    * @returns The reflections for a particular block. */
  /*    *1/ */
  /*   af::reflection_table split(std::size_t index) { */

  /*     using namespace af::boost_python::flex_table_suite; */

  /*     // Check the input */
  /*     DIALS_ASSERT(index < finished_.size()); */

  /*     // Get the indices of reflections to select */
  /*     std::vector<std::size_t> &process = to_process_[index]; */
  /*     std::vector<std::size_t> &include = to_include_[index]; */
  /*     af::shared<std::size_t> indices(process.size() + include.size()); */
  /*     std::copy(process.begin(), process.end(), indices.begin()); */
  /*     std::copy(include.begin(), include.end(), indices.begin() + process.size()); */

  /*     // Extract the reflection table */
  /*     af::reflection_table result = select_rows_index( */
  /*         data_, indices.const_ref()); */

  /*     // Extract the flags and set those reflections that are not to be */
  /*     // processed. */
  /*     af::ref<std::size_t> bk_id = result["bk_id"]; */
  /*     af::ref<std::size_t> flags = result["flags"]; */
  /*     for (std::size_t i = 0; i < bk_id.size(); ++i) { */
  /*       bk_id[i] = i; */
  /*     } */
  /*     for (std::size_t i = process.size(); i < flags.size(); ++i) { */
  /*       flags[i] |= af::DontIntegrate; */
  /*     } */

  /*     // Return the reflections */
  /*     return result; */
  /*   } */

  /*   /** */
  /*    * Accumulate the results. */
  /*    *1/ */
  /*   void accumulate(std::size_t index, af::reflection_table result) { */

  /*     using namespace af::boost_python::flex_table_suite; */

  /*     // Check the input */
  /*     DIALS_ASSERT(index < finished_.size()); */
  /*     DIALS_ASSERT(finished_[index] == false); */

  /*     // Get the indices of reflections to select */
  /*     std::vector<std::size_t> &process = to_process_[index]; */
  /*     std::vector<std::size_t> &include = to_include_[index]; */
  /*     af::const_ref<std::size_t> indices(&process[0], process.size()); */

  /*     // Resize the input to only select those which should have been processed. */
  /*     // Check that the book-keeping indices are as expected */
  /*     DIALS_ASSERT(process.size() + include.size() == result.size()); */
  /*     result.resize(process.size()); */
  /*     af::const_ref<std::size_t> bk_id = result["bk_id"]; */
  /*     for (std::size_t i = 0; i < bk_id.size(); ++i) { */
  /*       DIALS_ASSERT(bk_id[i] == i); */
  /*     } */
  /*     result.erase("bk_id"); */

  /*     // Set the result */
  /*     set_selected_rows_index(data_, indices, result); */

  /*     // Set finished flag */
  /*     finished_[index] = true; */
  /*   } */

  /* private: */

  /*   /** */
  /*    * Compute the blocks */
  /*    *1/ */
  /*   void compute_blocks( */
  /*       vec2<double> oscillation, */
  /*       vec2<int> array_range, */
  /*       double block_size, */
  /*       std::size_t max_procs, */
  /*       std::size_t max_tasks) { */
  /*     DIALS_ASSERT(block_size > 0); */
  /*     DIALS_ASSERT(max_procs > 0); */
  /*     DIALS_ASSERT(max_tasks > 0); */
  /*     double dphi = oscillation[1]; */
  /*     int frame0 = array_range[0]; */
  /*     int frame1 = array_range[1]; */
  /*     DIALS_ASSERT(frame1 > frame0); */
  /*     int nframes = frame1 - frame0; */
  /*     double half_block_size = block_size / 2.0; */
  /*     DIALS_ASSERT(half_block_size >= std::abs(dphi)); */
  /*     DIALS_ASSERT(half_block_size <= std::abs(nframes * dphi)); */
  /*     double half_block_length_f = half_block_size / dphi; */
  /*     int nblocks = (int)std::ceil(nframes / half_block_length_f); */
  /*     DIALS_ASSERT(nblocks > 0 && nblocks <= nframes); */
  /*     int half_block_length = (int)std::ceil((double)nframes / (double)nblocks); */
  /*     af::shared<int> indices; */
  /*     indices.push_back(frame0); */
  /*     for (int i = 0; i < nblocks; ++i) { */
  /*       int frame = frame0 + (i + 1) * half_block_length; */
  /*       if (frame > frame1) { */
  /*         frame = frame1; */
  /*       } */
  /*       indices.push_back(frame); */
  /*       if (frame == frame1) { */
  /*         break; */
  /*       } */
  /*     } */
  /*     DIALS_ASSERT(indices.front() == frame0); */
  /*     DIALS_ASSERT(indices.back() == frame1); */
  /*     DIALS_ASSERT(indices.size() > 2); */
  /*     std::vector< vec2<int> > jobs; */
  /*     for (std::size_t i = 0; i < indices.size() - 2; ++i) { */
  /*       jobs.push_back(vec2<int>(indices[i], indices[i+2])); */
  /*     } */
  /*     DIALS_ASSERT(jobs.size() > 0); */
  /*   } */

  /*   af::reflection_table data_; */
  /*   std::vector< vec2<int> > blocks_; */
  /*   std::vector< std::size_t> njobs_; */
  /*   std::vector<bool> finished_; */
  /*   std::vector< std::vector<std::size_t> > to_process_; */
  /*   std::vector< std::vector<std::size_t> > to_include_; */
  /*   std::vector<std::size_t> to_not_process_; */
  /* }; */

  /**
   * A class to managing spliting and mergin data
   */
  class IntegrationManagerData3D {
  public:

    IntegrationManagerData3D(
        af::reflection_table reflections,
        vec2<double> oscillation,
        vec2<int> array_range,
        double block_size)
          : data_(reflections) {

      // Compute the blocks
      compute_blocks(oscillation, array_range, block_size);
      finished_.assign(blocks_.size(), false);

      // Get the bounding boxes and flags
      af::const_ref<int6> bbox = reflections["bbox"];
      af::ref<std::size_t> flags = reflections["flags"];

      // Generate indices of reflections to be integrated, used as reference
      // spots or passed just as data for each data block. If the reflection is
      // not to be integrated, it is added to each block which it overlaps. If
      // the reflection is a reference spot, it is added to each block in which
      // it is fully recorded. If the spot is to be integrated, it is added to
      // the block in which it is closest to the centre. If the reflection is
      // larger than block_size / 2, then it is not fully recorded in any block
      // and is unprocessed.
      to_process_.resize(blocks_.size());
      to_include_.resize(blocks_.size());
      for (std::size_t i = 0; i < bbox.size(); ++i) {
        int z0 = bbox[i][4];
        int z1 = bbox[i][5];
        std::size_t &f = flags[i];
        if (!(f & af::DontIntegrate)) {
          std::size_t jmin = 0;
          double dmin = 0;
          for (std::size_t j = 0; j < blocks_.size(); ++j) {
            int bz0 = blocks_[j][0];
            int bz1 = blocks_[j][1];
            if (f & af::ReferenceSpot) {
              if (z0 >= bz0 && z1 <= bz1) {
                to_include_[j].push_back(i);
              }
            }
            double zc = (z1 + z0) / 2.0;
            double bc = (bz1 + bz0) / 2.0;
            double d = (zc - bc)*(zc - bc);
            if (j == 0 || d < dmin) {
              jmin = j;
              dmin = d;
            }
          }
          int bz0 = blocks_[jmin][0];
          int bz1 = blocks_[jmin][1];
          if (z0 >= bz0 && z1 <= bz1) {
            to_process_[jmin].push_back(i);
          } else {
            to_not_process_.push_back(i);
            f |= af::DontIntegrate;
          }
        }
        /* if (f & af::DontIntegrate) { */
        /*   for (std::size_t j = 0; j < blocks_.size(); ++j) { */
        /*     int bz0 = blocks_[j][0]; */
        /*     int bz1 = blocks_[j][1]; */
        /*     if (!(z1 <= bz0 || z0 >= bz1)) { */
        /*       to_include_[j].push_back(i); */
        /*     } */
        /*   } */
        /* } */
      }
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
    vec2<int> block(std::size_t index) const {
      DIALS_ASSERT(index < blocks_.size());
      return blocks_[index];
    }

    /**
     * @returns The list of reflections to not process.
     */
    af::shared<std::size_t> to_not_process() const {
      return af::shared<std::size_t>(
          &to_not_process_[0],
          &to_not_process_[0] + to_not_process_.size());
    }

    /**
     * @returns The list of reflections to include.
     */
    af::shared<std::size_t> to_include(std::size_t index) const {
      DIALS_ASSERT(index < blocks_.size());
      return af::shared<std::size_t>(
          &to_include_[index][0],
          &to_include_[index][0] + to_include_[index].size());
    }

    /**
     * @returns The list of reflections to process.
     */
    af::shared<std::size_t> to_process(std::size_t index) const {
      DIALS_ASSERT(index < blocks_.size());
      return af::shared<std::size_t>(
          &to_process_[index][0],
          &to_process_[index][0] + to_process_[index].size());
    }

    /**
     * @returns The reflections for a particular block.
     */
    af::reflection_table split(std::size_t index) {

      using namespace af::boost_python::flex_table_suite;

      // Check the input
      DIALS_ASSERT(index < finished_.size());

      // Get the indices of reflections to select
      std::vector<std::size_t> &process = to_process_[index];
      std::vector<std::size_t> &include = to_include_[index];
      af::shared<std::size_t> indices(process.size() + include.size());
      std::copy(process.begin(), process.end(), indices.begin());
      std::copy(include.begin(), include.end(), indices.begin() + process.size());

      // Extract the reflection table
      af::reflection_table result = select_rows_index(
          data_, indices.const_ref());

      // Extract the flags and set those reflections that are not to be
      // processed.
      af::ref<std::size_t> bk_id = result["bk_id"];
      af::ref<std::size_t> flags = result["flags"];
      for (std::size_t i = 0; i < bk_id.size(); ++i) {
        bk_id[i] = i;
      }
      for (std::size_t i = process.size(); i < flags.size(); ++i) {
        flags[i] |= af::DontIntegrate;
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

      // Get the indices of reflections to select
      std::vector<std::size_t> &process = to_process_[index];
      std::vector<std::size_t> &include = to_include_[index];
      af::const_ref<std::size_t> indices(&process[0], process.size());

      // Resize the input to only select those which should have been processed.
      // Check that the book-keeping indices are as expected
      DIALS_ASSERT(process.size() + include.size() == result.size());
      result.resize(process.size());
      af::const_ref<std::size_t> bk_id = result["bk_id"];
      for (std::size_t i = 0; i < bk_id.size(); ++i) {
        DIALS_ASSERT(bk_id[i] == i);
      }
      result.erase("bk_id");

      // Set the result
      set_selected_rows_index(data_, indices, result);

      // Set finished flag
      finished_[index] = true;
    }

  private:

    /**
     * Compute the blocks
     */
    void compute_blocks(
        vec2<double> oscillation,
        vec2<int> array_range,
        double block_size) {
      double dphi = oscillation[1];
      int frame0 = array_range[0];
      int frame1 = array_range[1];
      DIALS_ASSERT(frame1 > frame0);
      int nframes = frame1 - frame0;
      double half_block_size = block_size / 2.0;
      DIALS_ASSERT(half_block_size >= std::abs(dphi));
      DIALS_ASSERT(half_block_size <= std::abs(nframes * dphi));
      double half_block_length_f = half_block_size / dphi;
      int nblocks = (int)std::ceil(nframes / half_block_length_f);
      DIALS_ASSERT(nblocks > 0 && nblocks <= nframes);
      int half_block_length = (int)std::ceil((double)nframes / (double)nblocks);
      af::shared<int> indices;
      indices.push_back(frame0);
      for (int i = 0; i < nblocks; ++i) {
        int frame = frame0 + (i + 1) * half_block_length;
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
        blocks_.push_back(vec2<int>(indices[i], indices[i+2]));
      }
      DIALS_ASSERT(blocks_.size() > 0);
    }

    af::shared< vec2<int> > blocks_;
    af::shared<bool> finished_;
    af::reflection_table data_;
    std::vector< std::vector<std::size_t> > to_process_;
    std::vector< std::vector<std::size_t> > to_include_;
    std::vector<std::size_t> to_not_process_;
  };

}}

#endif // DIALS_ALGORITHMS_INTEGRATION_INTERFACE_H
