/*
 * interface.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/integration/processor.h>
#include <dials/algorithms/integration/integrator.h>
#include <dials/algorithms/integration/manager.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  /**
   * Split the reflection table where the blocks are given.
   */
  inline void job_list_split(const JobList &self, af::reflection_table data) {
    // Check the input
    DIALS_ASSERT(data.is_consistent());
    DIALS_ASSERT(data.contains("bbox"));
    DIALS_ASSERT(data.contains("id"));
    DIALS_ASSERT(self.size() > 0);
    DIALS_ASSERT(data.size() > 0);

    // Get the bounding boxes
    af::const_ref<int6> bbox = data["bbox"];
    af::const_ref<int> id = data["id"];

    // Check all the reflections are in range
    for (std::size_t i = 0; i < bbox.size(); ++i) {
      DIALS_ASSERT(id[i] >= 0);
      DIALS_ASSERT(bbox[i][1] > bbox[i][0]);
      DIALS_ASSERT(bbox[i][3] > bbox[i][2]);
      DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
    }

    // Create the lookup
    JobRangeLookup lookup(self);

    // Split the reflections
    af::shared<int6> bbox_new;
    af::shared<std::size_t> indices;
    for (std::size_t i = 0; i < bbox.size(); ++i) {
      int z0 = bbox[i][4];
      int z1 = bbox[i][5];
      std::size_t eid = id[i];
      std::size_t j0 = lookup.first(eid, z0);
      std::size_t j1 = lookup.last(eid, z1 - 1);
      DIALS_ASSERT(j0 < self.size());
      DIALS_ASSERT(j1 < self.size());
      DIALS_ASSERT(j1 >= j0);
      DIALS_ASSERT(z0 >= self[j0].frames()[0]);
      DIALS_ASSERT(z1 <= self[j1].frames()[1]);
      bool inside = false;
      for (std::size_t j = j0; j <= j1; ++j) {
        int jz0 = self[j].frames()[0];
        int jz1 = self[j].frames()[1];
        if (z0 >= jz0 && z1 <= jz1) {
          inside = true;
          break;
        }
      }
      if (inside) {
        bbox_new.push_back(bbox[i]);
        indices.push_back(i);
      } else {
        int6 b = bbox[i];
        std::vector<int> divisions;
        for (std::size_t j = j0; j <= j1; ++j) {
          divisions.push_back(self[j].frames()[0]);
          divisions.push_back(self[j].frames()[1]);
        }
        std::size_t k = 1;
        for (std::size_t j = 1; j < divisions.size(); ++j) {
          if (divisions[j] > divisions[j - 1]) {
            divisions[k] = divisions[j];
            k++;
          } else if (divisions[j] == divisions[j - 1]) {
            continue;
          } else {
            int a = divisions[j];
            int b = divisions[j - 1];
            int c = (a + b) / 2;
            DIALS_ASSERT(c >= a);
            DIALS_ASSERT(c < b);
            divisions[k] = c;
          }
        }
        divisions.resize(k);
        divisions[0] = b[4];
        k = 1;
        for (std::size_t j = 1; j < divisions.size(); ++j) {
          if (divisions[j] >= b[5]) {
            break;
          } else if (divisions[j] > divisions[j - 1]) {
            k++;
          } else {
            continue;
          }
        }
        divisions[k++] = b[5];
        divisions.resize(k);
        for (std::size_t j = 1; j < divisions.size(); ++j) {
          DIALS_ASSERT(divisions[j] > divisions[j - 1]);
        }
        for (std::size_t j = 1; j < divisions.size(); ++j) {
          b[5] = divisions[j];
          DIALS_ASSERT(b[5] > b[4]);
          bbox_new.push_back(b);
          indices.push_back(i);
          b[4] = b[5];
        }
      }
    }

    // Resize the reflection table
    DIALS_ASSERT(bbox_new.size() == indices.size());
    data.resize(bbox_new.size());

    // Reorder the reflections
    af::boost_python::flex_table_suite::reorder(data, indices.const_ref());

    // Set the new bounding boxes
    af::boost_python::flex_table_suite::setitem_column(
      data, "bbox", bbox_new.const_ref());
    af::boost_python::flex_table_suite::setitem_column(
      data, "partial_id", indices.const_ref());
  }

  /**
   * Compute the memory for each job
   */
  af::shared<std::size_t> job_list_shoebox_memory(const JobList &self,
                                                  af::reflection_table data,
                                                  bool flatten) {
    // Check the input
    DIALS_ASSERT(data.is_consistent());
    DIALS_ASSERT(data.contains("bbox"));
    DIALS_ASSERT(data.contains("id"));
    DIALS_ASSERT(data.contains("flags"));
    DIALS_ASSERT(self.size() > 0);
    DIALS_ASSERT(data.size() > 0);

    // Get the bounding boxes
    af::const_ref<int6> bbox = data["bbox"];
    af::const_ref<int> id = data["id"];
    af::const_ref<std::size_t> flags = data["flags"];

    // Check all the reflections are in range
    for (std::size_t i = 0; i < bbox.size(); ++i) {
      DIALS_ASSERT(id[i] >= 0);
      DIALS_ASSERT(bbox[i][1] > bbox[i][0]);
      DIALS_ASSERT(bbox[i][3] > bbox[i][2]);
      DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
    }

    // Create the lookup
    ReflectionLookup lookup(id, flags, bbox, self);

    // Compute the memory for each job
    af::shared<std::size_t> result(lookup.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      int frame0 = lookup.job(i).frames()[0];
      int frame1 = lookup.job(i).frames()[1];
      DIALS_ASSERT(frame1 > frame0);
      af::const_ref<std::size_t> indices = lookup.indices(i);
      af::shared<std::size_t> memory_to_alloc(frame1 - frame0, 0);
      af::shared<std::size_t> memory_to_free(frame1 - frame0, 0);
      for (std::size_t j = 0; j < indices.size(); ++j) {
        DIALS_ASSERT(indices[j] < bbox.size());
        int6 b = bbox[indices[j]];
        DIALS_ASSERT(b[4] >= frame0);
        DIALS_ASSERT(b[5] <= frame1);
        DIALS_ASSERT(b[5] > b[4]);
        DIALS_ASSERT(b[3] > b[2]);
        DIALS_ASSERT(b[1] > b[0]);
        std::size_t xsize = b[1] - b[0];
        std::size_t ysize = b[3] - b[2];
        std::size_t zsize = b[5] - b[4];
        std::size_t size = xsize * ysize;
        if (!flatten) {
          size *= zsize;
        }
        std::size_t nbytes = size
                             * (sizeof(Shoebox<>::float_type)
                                + sizeof(Shoebox<>::float_type) + sizeof(int));
        memory_to_alloc[b[4] - frame0] += nbytes;
        memory_to_free[b[5] - frame0 - 1] += nbytes;
      }
      std::size_t max_memory_usage = 0;
      std::size_t cur_memory_usage = 0;
      for (int frame = frame0; frame < frame1; ++frame) {
        std::size_t j = frame - frame0;
        cur_memory_usage += memory_to_alloc[j];
        DIALS_ASSERT(memory_to_free[j] <= cur_memory_usage);
        max_memory_usage = std::max(max_memory_usage, cur_memory_usage);
        cur_memory_usage -= memory_to_free[j];
      }
      DIALS_ASSERT(cur_memory_usage == 0);
      result[i] = max_memory_usage;
    }

    // Return the result
    return result;
  }

  /**
   * Wrapper class to allow python function to inherit
   */
  struct ExecutorWrapper : Executor, wrapper<Executor> {
    void process(int frame, af::reflection_table data) {
      this->get_override("process")(frame, data);
    }
  };

  BOOST_PYTHON_MODULE(dials_algorithms_integration_integrator_ext) {
    class_<GroupList::Group>("Group", no_init)
      .def("index", &GroupList::Group::index)
      .def("nindex", &GroupList::Group::nindex)
      .def("expr", &GroupList::Group::expr)
      .def("nexpr", &GroupList::Group::nexpr)
      .def("frames", &GroupList::Group::frames)
      .def("nframes", &GroupList::Group::nframes);

    class_<GroupList>("GroupList")
      .def("__len__", &GroupList::size)
      .def("__getitem__", &GroupList::operator[], return_internal_reference<>());

    class_<JobList::Job>("Job", no_init)
      .def("index", &JobList::Job::index)
      .def("expr", &JobList::Job::expr)
      .def("nexpr", &JobList::Job::nexpr)
      .def("frames", &JobList::Job::frames)
      .def("nframes", &JobList::Job::nframes);

    class_<JobList>("JobList")
      .def(init<tiny<int, 2>, const af::const_ref<tiny<int, 2> > &>())
      .def("add", &JobList::add)
      .def("__len__", &JobList::size)
      .def("__getitem__", &JobList::operator[], return_internal_reference<>())
      .def("split", &job_list_split)
      .def("shoebox_memory", &job_list_shoebox_memory);

    class_<ReflectionManager>("ReflectionManager", no_init)
      .def(init<const JobList &, af::reflection_table>((arg("jobs"), arg("data"))))
      .def("__len__", &ReflectionManager::size)
      .def("finished", &ReflectionManager::finished)
      .def("accumulate", &ReflectionManager::accumulate)
      .def("split", &ReflectionManager::split)
      .def("job", &ReflectionManager::job, return_internal_reference<>())
      .def("data", &ReflectionManager::data)
      .def("num_reflections", &ReflectionManager::num_reflections);

    class_<ReflectionManagerPerImage>("ReflectionManagerPerImage", no_init)
      .def(init<int2, af::reflection_table>((arg("frames"), arg("data"))))
      .def("__len__", &ReflectionManagerPerImage::size)
      .def("finished", &ReflectionManagerPerImage::finished)
      .def("accumulate", &ReflectionManagerPerImage::accumulate)
      .def("split", &ReflectionManagerPerImage::split)
      .def("frames", &ReflectionManagerPerImage::frames)
      .def("data", &ReflectionManagerPerImage::data)
      .def("num_reflections", &ReflectionManagerPerImage::num_reflections);

    class_<ExecutorWrapper, boost::noncopyable>("Executor")
      .def("process", pure_virtual(&Executor::process))
      .enable_pickling();

    class_<ShoeboxProcessor>("ShoeboxProcessor", no_init)
      .def(init<af::reflection_table, std::size_t, int, int, bool>())
      .def("compute_max_memory_usage", &ShoeboxProcessor::compute_max_memory_usage)
      .def("next", &ShoeboxProcessor::next<double>)
      .def("next", &ShoeboxProcessor::next<int>)
      .def("frame0", &ShoeboxProcessor::frame0)
      .def("frame1", &ShoeboxProcessor::frame1)
      .def("frame", &ShoeboxProcessor::frame)
      .def("nframes", &ShoeboxProcessor::nframes)
      .def("npanels", &ShoeboxProcessor::npanels)
      .def("finished", &ShoeboxProcessor::finished)
      .def("extract_time", &ShoeboxProcessor::extract_time)
      .def("process_time", &ShoeboxProcessor::process_time);
  }

}}}  // namespace dials::algorithms::boost_python
