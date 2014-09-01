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


  class IntegrationTask3DSpec {
  public:

    IntegrationTask3DSpec(
            af::reflection_table data,
            std::size_t npanels,
            const af::const_ref< tiny<int,2> > &jobs,
            const af::const_ref< std::size_t > &off,
            const af::const_ref< std::size_t > &ind)
        : data_(data),
          npanels_(npanels),
          jobs_(jobs.begin(), jobs.end()),
          offset_(off.begin(), off.end()),
          indices_(ind.begin(), ind.end()) {

      // Check some input
      DIALS_ASSERT(data.size() > 0);
      DIALS_ASSERT(jobs.size() > 0);
      DIALS_ASSERT(off.size() > 0);
      DIALS_ASSERT(ind.size() > 0);

      // Check the jobs are valid. Jobs must be ordered in increasing frame
      // number and can overlap but 1 jobs must not be fully contained in
      // another. There must also be no gaps.
      for (std::size_t i = 0; i < jobs.size()-1; ++i) {
        DIALS_ASSERT(jobs[i][0] <  jobs[i][1]);
        DIALS_ASSERT(jobs[i][0] <  jobs[i+1][0]);
        DIALS_ASSERT(jobs[i][1] <  jobs[i+1][1]);
        DIALS_ASSERT(jobs[i][1] >= jobs[i+1][0]);
      }

      // Get the range of frames
      frame0_ = jobs_.front()[0];
      frame1_ = jobs_.back()[1];
      DIALS_ASSERT(frame0_ < frame1_);
      nframes_ = frame1_ - frame0_;

      // Check the reflection table contains a few required fields
      DIALS_ASSERT(data.is_consistent());
      DIALS_ASSERT(data.contains("panel"));
      DIALS_ASSERT(data.contains("bbox"));

      // Get some data from the reflection table
      af::const_ref< std::size_t > panel = data["panel"];
      af::const_ref< int6 > bbox = data["bbox"];

      // Check the jobs ids are valid
      DIALS_ASSERT(offset_.size() == njobs()+1);
      DIALS_ASSERT(offset_.front() == 0);
      DIALS_ASSERT(offset_.back() == indices_.size());
      for (std::size_t i = 0; i < njobs(); ++i) {
        af::const_ref<std::size_t> ind = indices(i);
        for (std::size_t j = 0; j < ind.size(); ++j) {
          std::size_t k = ind[j];
          std::size_t p = panel[k];
          int z0 = bbox[k][4];
          int z1 = bbox[k][5];
          DIALS_ASSERT(p < npanels_);
          DIALS_ASSERT(z1 > z0);
          int j0 = job(i)[0];
          int j1 = job(i)[1];
          DIALS_ASSERT(j1 > j0);
          DIALS_ASSERT(z0 >= j0 && z1 <= j1);
        }
      }
    }

    af::reflection_table data() const {
      return data_;
    }

    tiny<int,2> job(std::size_t index) const {
      return jobs_[index];
    }

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

    int frame0() const {
      return frame0_;
    }

    int frame1() const {
      return frame1_;
    }

    std::size_t nframes() const {
      return nframes_;
    }

    std::size_t npanels() const {
      return npanels_;
    }

    std::size_t njobs() const {
      return jobs_.size();
    }

  private:

    af::reflection_table data_;
    std::size_t npanels_;
    std::vector< tiny<int,2> > jobs_;
    std::vector< std::size_t > offset_;
    std::vector< std::size_t > indices_;
    int frame0_;
    int frame1_;
    std::size_t nframes_;
  };


  class IntegrationTask3DAllocator {
  public:

    IntegrationTask3DAllocator(const IntegrationTask3DSpec &spec) {

      typedef std::pair<std::size_t,std::size_t> active_type;

      std::vector<active_type> active(spec.data().size());
      std::vector<bool> first(active.size(), true);

      for (std::size_t i = 0; i < spec.njobs(); ++i) {
        af::const_ref<std::size_t> ind = spec.indices(i);
        for (std::size_t j = 0; j < ind.size(); ++j) {
          std::size_t k = ind[j];
          active_type &a = active[k];
          if (first[k]) {
            a.first  = i;
            a.second = i;
            first[k] = false;
          } else {
            if (i < a.first)  a.first  = i;
            if (i > a.second) a.second = i;
          }
        }
      }

      std::vector<std::size_t> count1(spec.njobs(), 0);
      std::vector<std::size_t> count2(spec.njobs(), 0);
      std::vector<std::size_t> num1(spec.njobs(), 0);
      std::vector<std::size_t> num2(spec.njobs(), 0);
      for (std::size_t i = 0; i < active.size(); ++i) {
        num1[active[i].first]++;
        num2[active[i].second]++;
      }
      malloc_offset_.resize(spec.njobs() + 1);
      malloc_indices_.resize(spec.data().size());
      free_offset_.resize(spec.njobs() + 1);
      free_indices_.resize(spec.data().size());
      malloc_offset_[0] = 0;
      free_offset_[0] = 0;
      for (std::size_t i = 1; i < malloc_offset_.size(); ++i) {
        malloc_offset_[i] = malloc_offset_[i-1] + num1[i-1];
        free_offset_[i] = free_offset_[i-1] + num2[i-1];
      }
      DIALS_ASSERT(malloc_offset_.back() == malloc_indices_.size());
      DIALS_ASSERT(free_offset_.back() == free_indices_.size());
      for (std::size_t i = 0; i < active.size(); ++i) {
        std::size_t i1 = active[i].first;
        std::size_t i2 = active[i].second;
        std::size_t k1 = malloc_offset_[i1] + count1[i1];
        std::size_t k2 = malloc_offset_[i2] + count2[i2];
        malloc_indices_[k1] = i;
        free_indices_[k2] = i;
        count1[i1]++;
        count2[i2]++;
      }
      for (std::size_t i =0 ;i < count1.size(); ++i) {
        DIALS_ASSERT(count1[i] == num1[i]);
        DIALS_ASSERT(count2[i] == num2[i]);
      }
    }

    af::const_ref<std::size_t> malloc(std::size_t index) const {
      DIALS_ASSERT(index < malloc_offset_.size()-1);
      std::size_t i0 = malloc_offset_[index];
      std::size_t i1 = malloc_offset_[index+1];
      DIALS_ASSERT(i1 >= i0);
      std::size_t off = i0;
      std::size_t num = i1 - i0;
      DIALS_ASSERT(off + num <= malloc_indices_.size());
      return af::const_ref<std::size_t> (&malloc_indices_[off], num);
    }

    af::const_ref<std::size_t> free(std::size_t index) const {
      DIALS_ASSERT(index < free_offset_.size()-1);
      std::size_t i0 = free_offset_[index];
      std::size_t i1 = free_offset_[index+1];
      DIALS_ASSERT(i1 >= i0);
      std::size_t off = i0;
      std::size_t num = i1 - i0;
      DIALS_ASSERT(off + num <= free_indices_.size());
      return af::const_ref<std::size_t> (&free_indices_[off], num);
    }

  private:

    std::vector<std::size_t> malloc_offset_;
    std::vector<std::size_t> malloc_indices_;
    std::vector<std::size_t> free_offset_;
    std::vector<std::size_t> free_indices_;
  };

  class IntegrationTask3DExecutor {
  public:

    typedef af::reflection_table rtable;
    typedef boost::function<rtable (rtable)> callback_type;

    IntegrationTask3DExecutor(
            const IntegrationTask3DSpec &spec,
            callback_type callback)
        : spec_(spec),
          allocator_(spec),
          process_(callback) {

      frame_ = spec.frame0();

      // Set the current active job
      begin_active_ = 0;
      end_active_ = 0;

      // Initialise the shoeboxes
      af::const_ref< std::size_t > panel = spec_.data()["panel"];
      af::const_ref< int6 > bbox = spec_.data()["bbox"];
      shoebox_ = spec_.data()["shoebox"];
      for (std::size_t i = 0; i < shoebox_.size(); ++i) {
        shoebox_[i] = Shoebox<>(panel[i], bbox[i]);
      }

      // Initialise the offsets and indices for each frame/panel
      initialise_indices();
      /* malloc(0); */
    }

    int frame0() const {
      return spec_.frame0();
    }

    int frame1() const {
      return spec_.frame1();
    }

    int frame() const {
      return frame_;
    }

    std::size_t nframes() const {
      return spec_.nframes();
    }

    std::size_t npanels() const {
      return spec_.npanels();
    }

    void next(const Image &image) {
      DIALS_ASSERT(!finished());
      if (first_image_in_job()) {
        begin_job();
      }
      next_image(image);
      if (last_image_in_job()) {
        end_job();
      }
    }

    af::reflection_table data() {
      DIALS_ASSERT(finished());
      spec_.data().erase("shoebox");
      return spec_.data();
    }

    bool finished() const {
      return frame_ == frame1();
    }

  private:

    void initialise_indices() {
      typedef std::list<std::size_t> list_type;
      typedef list_type::iterator iterator;

      // Temporary index arrays
      std::size_t num = nframes() * npanels();
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
          int j = p + (z - frame0())*npanels()+1;
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

    void next_image(const Image &image) {
      DIALS_ASSERT(image.npanels() == npanels());
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
      std::size_t j0 = panel+(frame-frame0())*npanels();
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

    bool first_image_in_job() const {
      if (begin_active_ < spec_.njobs()) {
        int j0 = spec_.job(begin_active_)[0];
        DIALS_ASSERT(frame_ <= j0);
        return frame_ == j0;
      }
      return false;
    }

    bool last_image_in_job() const {
      DIALS_ASSERT(end_active_ < spec_.njobs());
      int j0 = spec_.job(end_active_)[0];
      int j1 = spec_.job(end_active_)[1];
      DIALS_ASSERT(frame_ >= j0 && frame_ <= j1);
      return frame_ == j1;
    }

    void begin_job() {
      malloc(begin_active_);
      begin_active_++;
    }

    void end_job() {
      free(end_active_);
      merge(process_(split(end_active_)), end_active_);
      end_active_++;
    }

    void free(std::size_t job) {
      af::const_ref<std::size_t> ind = allocator_.free(job);
      for (std::size_t i = 0; i < ind.size(); ++i) {
        shoebox_[ind[i]].deallocate();
      }
    }

    void malloc(std::size_t job) {
      af::const_ref<std::size_t> ind = allocator_.malloc(job);
      for (std::size_t i = 0; i < ind.size(); ++i) {
        shoebox_[ind[i]].allocate();
      }
      /* int z0 = spec_.job(job)[0]; */
      /* int z1 = spec_.job(job)[1]; */
      /* for (std::size_t j = job; j < njobs(); ++j) { */
      /*   int z2 = spec_.job(j)[0]; */
      /*   if (z2 < z1) { */
      /*     af::const_ref<std::size_t> ind = free(job_indices(j); */
      /*     for (std::size_t i = 0; i < ind.size(); ++i) { */
      /*       shoebox_[ind[i]].allocate(); */
      /*     } */
      /*   } */
      /* } */
    }

    af::reflection_table split(std::size_t job) const {
      using namespace dials::af::boost_python::flex_table_suite;
      return select_rows_index(spec_.data(), spec_.indices(job));
    }

    void merge(af::reflection_table result, std::size_t job) {
      using namespace dials::af::boost_python::flex_table_suite;
      af::reflection_table data = spec_.data();
      set_selected_rows_index(data, spec_.indices(job), result);
    }

    IntegrationTask3DSpec spec_;
    IntegrationTask3DAllocator allocator_;
    int frame_;
    std::size_t begin_active_;
    std::size_t end_active_;
    callback_type process_;
    af::shared< Shoebox<> > shoebox_;
    std::vector<std::size_t> offset_;
    std::vector<std::size_t> indices_;
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
