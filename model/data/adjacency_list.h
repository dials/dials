/*
 * adjacency_list.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_ADJACENCY_LIST_H
#define DIALS_MODEL_ADJACENCY_LIST_H

#include <iostream>
#include <deque>
#include <vector>
#include <dials/error.h>

namespace dials { namespace model {

  // Create the adjacency list type
  // typedef boost::adjacency_list<
  //  boost::listS,
  //  boost::vecS,
  //  boost::undirectedS> AdjacencyList;

  // typedef boost_adaptbx::graph_type::adjacency_list_undirected_vecS_setS_type
  // AdjacencyList;

  class AdjacencyList {
  public:
    typedef std::pair<std::size_t, std::size_t> edge_descriptor;
    typedef std::deque<edge_descriptor> edge_list;
    typedef edge_list::const_iterator edge_iterator;
    typedef std::pair<edge_iterator, edge_iterator> edge_iterator_range;

    AdjacencyList(std::size_t num_vertices)
        : offset_(num_vertices + 1), num_vertices_(num_vertices), consistent_(false) {}

    std::size_t source(edge_descriptor edge) const {
      DIALS_ASSERT(consistent_);
      return edge.first;
    }

    std::size_t target(edge_descriptor edge) const {
      DIALS_ASSERT(consistent_);
      return edge.second;
    }

    edge_iterator_range edges() const {
      DIALS_ASSERT(consistent_);
      return edge_iterator_range(edges_.begin(), edges_.end());
    }

    edge_iterator_range edges(std::size_t i) const {
      DIALS_ASSERT(consistent_);
      DIALS_ASSERT(i < num_vertices());
      DIALS_ASSERT(i < offset_.size() - 1);
      std::size_t o1 = offset_[i];
      std::size_t o2 = offset_[i + 1];
      DIALS_ASSERT(o2 >= o1);
      DIALS_ASSERT(o2 <= edges_.size());
      return edge_iterator_range(edges_.begin() + o1, edges_.begin() + o2);
    }

    void add_edge(std::size_t a, std::size_t b) {
      consistent_ = false;
      DIALS_ASSERT(num_vertices());
      DIALS_ASSERT(num_vertices());
      edges_.push_back(edge_descriptor(a, b));
      edges_.push_back(edge_descriptor(b, a));
    }

    void finish() {
      std::sort(edges_.begin(), edges_.end());
      std::vector<std::size_t> num(num_vertices());
      for (std::size_t i = 0; i < edges_.size(); ++i) {
        num[edges_[i].first]++;
      }
      offset_[0] = 0;
      for (std::size_t i = 0; i < num.size(); ++i) {
        offset_[i + 1] = offset_[i] + num[i];
      }
      DIALS_ASSERT(offset_.back() == edges_.size());
      consistent_ = true;
    }

    std::size_t num_vertices() const {
      return num_vertices_;
    }

    std::size_t num_edges() const {
      DIALS_ASSERT((edges_.size() & 1) == 0);
      return edges_.size() / 2;
    }

    std::size_t vertex_num_edges(std::size_t i) const {
      DIALS_ASSERT(i < offset_.size() - 1);
      std::size_t o1 = offset_[i];
      std::size_t o2 = offset_[i + 1];
      DIALS_ASSERT(o2 >= o1);
      return o2 - o1;
    }

  private:
    edge_list edges_;
    std::vector<std::size_t> offset_;
    std::size_t num_vertices_;
    bool consistent_;
  };

}}  // namespace dials::model

#endif /* DIALS_MODEL_DATA_ADJACENCY_LIST_H */
