/*
 * helpers.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPOT_FINDING_HELPERS_H
#define DIALS_ALGORITHMS_SPOT_FINDING_HELPERS_H

#include <dials/array_family/reflection_table.h>
#include <dials/array_family/boost_python/flex_table_suite.h>
#include <dials/algorithms/image/connected_components/connected_components.h>

namespace dials { namespace algorithms {

  using dials::model::Shoebox;

  /**
   * Relabel some shoeboxes
   */
  class Labeller {
  public:
    Labeller() {}

    void add(const Shoebox<> &shoebox) {
      // Add to the list of shoeboxes
      shoeboxes_.push_back(shoebox);
    }

    af::shared<Shoebox<> > shoeboxes() {
      typedef Shoebox<>::float_type float_type;

      // Adjacency list type
      typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>
        AdjacencyList;

      // Create the edges between shoeboxes
      AdjacencyList graph(shoeboxes_.size());
      for (std::size_t s1 = 0; s1 < shoeboxes_.size() - 1; ++s1) {
        std::size_t panel1 = shoeboxes_[s1].panel;
        int6 bbox1 = shoeboxes_[s1].bbox;
        for (std::size_t s2 = s1 + 1; s2 < shoeboxes_.size(); ++s2) {
          std::size_t panel2 = shoeboxes_[s2].panel;
          int6 bbox2 = shoeboxes_[s2].bbox;
          if (panel1 == panel2 && bbox1[5] == bbox2[4]
              && (bbox1[0] < bbox2[1] && bbox1[1] > bbox2[0] && bbox1[2] < bbox2[3]
                  && bbox1[3] > bbox2[2])) {
            af::const_ref<int, af::c_grid<3> > mask1 = shoeboxes_[s1].mask.const_ref();
            af::const_ref<int, af::c_grid<3> > mask2 = shoeboxes_[s2].mask.const_ref();
            int x0 = std::max(bbox1[0], bbox2[0]);
            int x1 = std::min(bbox1[1], bbox2[1]);
            int y0 = std::max(bbox1[2], bbox2[2]);
            int y1 = std::min(bbox1[3], bbox2[3]);
            int k1 = mask1.accessor()[0] - 1;
            int k2 = 0;
            for (int y = y0; y < y1; ++y) {
              for (int x = x0; x < x1; ++x) {
                int i1 = x - bbox1[0];
                int i2 = x - bbox2[0];
                int j1 = y - bbox1[2];
                int j2 = y - bbox2[2];
                DIALS_ASSERT(i1 >= 0 && i1 < mask1.accessor()[2]);
                DIALS_ASSERT(i2 >= 0 && i2 < mask2.accessor()[2]);
                DIALS_ASSERT(j1 >= 0 && j1 < mask1.accessor()[1]);
                DIALS_ASSERT(j2 >= 0 && j2 < mask2.accessor()[1]);
                if ((mask1(k1, j1, i1) & Foreground)
                    && (mask2(k2, j2, i2) & Foreground)) {
                  boost::add_edge(s1, s2, graph);
                }
              }
            }
          }
        }
      }

      // Do the connected components
      af::shared<int> labels(num_vertices(graph), af::init_functor_null<int>());
      DIALS_ASSERT(labels.size() == shoeboxes_.size());
      int num = boost::connected_components(graph, &labels[0]);
      DIALS_ASSERT(num <= labels.size());

      // Get the number of labels and allocate the array
      std::size_t max_label = af::max(labels.const_ref());
      af::shared<Shoebox<> > result(max_label + 1, Shoebox<>());

      // Initialise the bboxes
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i].panel = -1;
      }

      // Set the shoeboxes
      for (std::size_t i = 0; i < labels.size(); ++i) {
        int l = labels[i];
        Shoebox<> s = shoeboxes_[i];
        if (result[l].panel == -1) {
          result[l].panel = s.panel;
          result[l].bbox = s.bbox;
        } else {
          if (s.bbox[0] < result[l].bbox[0]) result[l].bbox[0] = s.bbox[0];
          if (s.bbox[1] > result[l].bbox[1]) result[l].bbox[1] = s.bbox[1];
          if (s.bbox[2] < result[l].bbox[2]) result[l].bbox[2] = s.bbox[2];
          if (s.bbox[3] > result[l].bbox[3]) result[l].bbox[3] = s.bbox[3];
          if (s.bbox[4] < result[l].bbox[4]) result[l].bbox[4] = s.bbox[4];
          if (s.bbox[5] > result[l].bbox[5]) result[l].bbox[5] = s.bbox[5];
        }
      }

      // Allocate all the arrays
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i].allocate();
      }

      // Set all the mask and data points
      for (std::size_t i = 0; i < labels.size(); ++i) {
        int l = labels[i];
        int6 bbox1 = shoeboxes_[i].bbox;
        int6 bbox2 = result[l].bbox;
        DIALS_ASSERT(result[l].is_consistent());
        DIALS_ASSERT(shoeboxes_[i].is_consistent());
        af::const_ref<int, af::c_grid<3> > mask1 = shoeboxes_[i].mask.const_ref();
        af::const_ref<float_type, af::c_grid<3> > data1 =
          shoeboxes_[i].data.const_ref();
        af::ref<int, af::c_grid<3> > mask2 = result[l].mask.ref();
        af::ref<float_type, af::c_grid<3> > data2 = result[l].data.ref();
        for (int z = bbox1[4]; z < bbox1[5]; ++z) {
          for (int y = bbox1[2]; y < bbox1[3]; ++y) {
            for (int x = bbox1[0]; x < bbox1[1]; ++x) {
              int i1 = x - bbox1[0];
              int i2 = x - bbox2[0];
              int j1 = y - bbox1[2];
              int j2 = y - bbox2[2];
              int k1 = z - bbox1[4];
              int k2 = z - bbox2[4];
              DIALS_ASSERT(i1 >= 0 && i1 < mask1.accessor()[2]);
              DIALS_ASSERT(i2 >= 0 && i2 < mask2.accessor()[2]);
              DIALS_ASSERT(j1 >= 0 && j1 < mask1.accessor()[1]);
              DIALS_ASSERT(j2 >= 0 && j2 < mask2.accessor()[1]);
              DIALS_ASSERT(k1 >= 0 && k1 < mask1.accessor()[0]);
              DIALS_ASSERT(k2 >= 0 && k2 < mask2.accessor()[0]);
              if (mask1(k1, j1, i1) & Foreground) {
                mask2(k2, j2, i2) = mask1(k1, j1, i1);
                data2(k2, j2, i2) = data1(k1, j1, i1);
              }
            }
          }
        }
      }

      return result;
    }

    void clear() {
      shoeboxes_ = af::shared<Shoebox<> >();
    }

  private:
    af::shared<Shoebox<> > shoeboxes_;
  };

  /**
   * A class to combine strong spot lists (i.e. perform the connected component
   * labelling at the boundary of two spot lists from the same sweep.
   */
  class StrongSpotCombiner {
  public:
    /**
     * Initialise everything
     */
    StrongSpotCombiner() : all_minz_(0), all_maxz_(0) {}

    /**
     * Add a reflection table to the combiner
     * @param rlist The reflection table
     */
    void add(af::const_ref<Shoebox<> > shoebox) {
      // Find the min and max frame
      int minz = shoebox[0].bbox[4];
      int maxz = shoebox[0].bbox[5];
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        if (shoebox[i].bbox[4] < minz) {
          minz = shoebox[i].bbox[4];
        }
        if (shoebox[i].bbox[5] > maxz) {
          maxz = shoebox[i].bbox[5];
        }
      }

      // If this is the first reflection table, then copy, otherwise, check that
      // the frame frames match.
      if (all_maxz_ == all_minz_) {
        all_minz_ = minz;
        all_maxz_ = maxz;
      } else {
        DIALS_ASSERT(minz == all_maxz_);
        all_maxz_ = maxz;
      }

      // Loop through all the shoeboxes. If the shoebox is not at the edge of
      // the image, then add it to the finished list. Otherwise, add it to the
      // labeller.
      for (std::size_t i = 0; i < shoebox.size(); ++i) {
        if (minz < shoebox[i].bbox[4] && shoebox[i].bbox[5] < maxz) {
          finished_.push_back(shoebox[i]);
        } else {
          labeller_.add(shoebox[i]);
        }
      }
    }

    /**
     * @returns The final reflection table
     */
    af::shared<Shoebox<> > shoeboxes() {
      af::shared<Shoebox<> > labelled = labeller_.shoeboxes();
      finished_.insert(finished_.end(), labelled.begin(), labelled.end());
      labeller_.clear();
      return finished_;
    }

  private:
    Labeller labeller_;
    af::shared<Shoebox<> > finished_;
    int all_minz_;
    int all_maxz_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPOT_FINDING_HELPERS_H
