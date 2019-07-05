/*
 * index.h
 *
 *  Copyright (C) 2014 Diamond Light Source
 *
 *  Author: Richard Gildea
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_INDEXING_H
#define DIALS_ALGORITHMS_INDEXING_H
#include <vector>
#include <map>
#include <algorithm>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/math/utils.h>
#include <cctbx/miller.h>
#include <scitbx/constants.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/connected_components.hpp>
#include <annlib_adaptbx/ann_adaptor.h>

namespace dials { namespace algorithms {

  class AssignIndices {
  public:
    AssignIndices(af::const_ref<scitbx::vec3<double> > const& reciprocal_space_points,
                  af::const_ref<double> const& phi,
                  af::const_ref<scitbx::mat3<double> > const& UB_matrices,
                  double tolerance = 0.3)
        : miller_indices_(reciprocal_space_points.size(),
                          cctbx::miller::index<>(0, 0, 0)),
          crystal_ids_(reciprocal_space_points.size(), -1) {
      DIALS_ASSERT(reciprocal_space_points.size() == phi.size());

      typedef std::pair<const cctbx::miller::index<>, const std::size_t> pair_t;
      typedef std::multimap<cctbx::miller::index<>, std::size_t> map_t;

      map_t hkl_to_rlp_map;

      std::vector<af::shared<cctbx::miller::index<> > > hkl_ints;
      std::vector<af::shared<double> > lengths_sq;

      const double pi_4 = scitbx::constants::pi / 4;

      // loop over crystals and assign one hkl per crystal per reflection
      for (int i_lattice = 0; i_lattice < UB_matrices.size(); i_lattice++) {
        scitbx::mat3<double> A = UB_matrices[i_lattice];
        scitbx::mat3<double> A_inv = A.inverse();
        af::shared<cctbx::miller::index<> > hkl_ints_(
          af::reserve(reciprocal_space_points.size()));
        af::shared<double> lengths_sq_(af::reserve(reciprocal_space_points.size()));
        for (int i_ref = 0; i_ref < reciprocal_space_points.size(); i_ref++) {
          scitbx::vec3<double> rlp = reciprocal_space_points[i_ref];
          scitbx::vec3<double> hkl_f = A_inv * rlp;
          cctbx::miller::index<> hkl_i;
          for (std::size_t j = 0; j < 3; j++) {
            hkl_i[j] = scitbx::math::iround(hkl_f[j]);
          }
          scitbx::vec3<double> diff = hkl_f - scitbx::vec3<double>(hkl_i);
          hkl_ints_.push_back(hkl_i);
          lengths_sq_.push_back(diff.length_sq());
        }
        hkl_ints.push_back(hkl_ints_);
        lengths_sq.push_back(lengths_sq_);
      }

      // loop over all reflections and choose the best hkl (and consequently
      // crystal) for each reflection
      double tolerance_sq = tolerance * tolerance;
      for (int i_ref = 0; i_ref < reciprocal_space_points.size(); i_ref++) {
        af::shared<double> n;
        af::shared<cctbx::miller::index<> > potential_hkls;
        for (int i_lattice = 0; i_lattice < UB_matrices.size(); i_lattice++) {
          n.push_back(lengths_sq[i_lattice][i_ref]);
          potential_hkls.push_back(hkl_ints[i_lattice][i_ref]);
        }
        int i_best_lattice = af::min_index(n.const_ref());
        if (n[i_best_lattice] > tolerance_sq) {
          continue;
        }
        cctbx::miller::index<> hkl = potential_hkls[i_best_lattice];
        if (hkl[0] == 0 && hkl[1] == 0 && hkl[2] == 0) {
          continue;
        }
        miller_indices_[i_ref] = hkl;
        hkl_to_rlp_map.insert(pair_t(hkl, i_ref));
        crystal_ids_[i_ref] = i_best_lattice;
      }

      cctbx::miller::index<> curr_hkl(0, 0, 0);
      std::vector<std::size_t> i_same_hkl;
      // if more than one spot can be assigned the same miller index then
      // choose the closest one
      for (map_t::iterator it = hkl_to_rlp_map.begin(); it != hkl_to_rlp_map.end();
           it++) {
        if (it->first == cctbx::miller::index<>(0, 0, 0)) {
          continue;
        }
        if (it->first != curr_hkl) {
          if (i_same_hkl.size() > 1) {
            for (int i = 0; i < i_same_hkl.size(); i++) {
              const std::size_t i_ref = i_same_hkl[i];
              for (int j = i + 1; j < i_same_hkl.size(); j++) {
                const std::size_t j_ref = i_same_hkl[j];
                int crystal_i = crystal_ids_[i_ref];
                int crystal_j = crystal_ids_[j_ref];
                if (crystal_i != crystal_j) {
                  continue;
                } else if (crystal_i == -1) {
                  continue;
                }
                double phi_i = phi[i_ref];
                double phi_j = phi[j_ref];
                if (std::abs(phi_i - phi_j) > pi_4) {
                  continue;
                }
                if (lengths_sq[crystal_j][j_ref] < lengths_sq[crystal_i][i_ref]) {
                  miller_indices_[i_ref] = cctbx::miller::index<>(0, 0, 0);
                  crystal_ids_[i_ref] = -1;
                } else {
                  miller_indices_[j_ref] = cctbx::miller::index<>(0, 0, 0);
                  crystal_ids_[j_ref] = -1;
                }
              }
            }
          }
          curr_hkl = it->first;
          i_same_hkl.clear();
        }
        i_same_hkl.push_back(it->second);
      }
    }

    af::shared<cctbx::miller::index<> > miller_indices() {
      return miller_indices_;
    }

    af::shared<int> crystal_ids() {
      return crystal_ids_;
    }

  private:
    af::shared<cctbx::miller::index<> > miller_indices_;
    af::shared<int> crystal_ids_;
  };

  typedef struct edge_ {
    std::size_t i;
    std::size_t j;
    double l_ij;
    cctbx::miller::index<> h_ij;
  } MyEdge;

  template <typename Edge>
  struct record_dfs_order : public boost::default_dfs_visitor {
    record_dfs_order(std::vector<Edge>& e) : edges(e) {}

    template <typename Graph>
    void tree_edge(Edge e, const Graph& G) const {
      edges.push_back(e);
    }
    std::vector<Edge>& edges;
  };

  class AssignIndicesLocal {
  public:
    AssignIndicesLocal(
      af::const_ref<scitbx::vec3<double> > const& reciprocal_space_points,
      af::const_ref<double> const& phi,
      af::const_ref<scitbx::mat3<double> > const& UB_matrices,
      const double epsilon = 0.05,
      const double delta = 5,
      const double l_min = 0.8,
      const int nearest_neighbours = 20)
        : miller_indices_(reciprocal_space_points.size(),
                          cctbx::miller::index<>(0, 0, 0)),
          crystal_ids_(reciprocal_space_points.size(), -1) {
      DIALS_ASSERT(reciprocal_space_points.size() == phi.size());

      using annlib_adaptbx::AnnAdaptor;
      using namespace boost;

      typedef adjacency_list<vecS, vecS, undirectedS, no_property, MyEdge> Graph;

      typedef Graph::vertex_descriptor Vertex;
      typedef graph_traits<Graph>::edge_descriptor Edge;

      // convert into a single array for input to AnnAdaptor
      // based on flex.vec_3.as_double()
      // scitbx/array_family/boost_python/flex_vec3_double.cpp
      af::shared<double> rlps_double(reciprocal_space_points.size() * 3,
                                     af::init_functor_null<double>());
      double* r = rlps_double.begin();
      for (std::size_t i = 0; i < reciprocal_space_points.size(); i++) {
        for (std::size_t j = 0; j < 3; j++) {
          *r++ = reciprocal_space_points[i][j];
        }
      }

      AnnAdaptor ann = AnnAdaptor(rlps_double, 3, nearest_neighbours);
      ann.query(rlps_double);

      typedef std::pair<const int, const int> pair_t;

      const double one_over_epsilon = 1.0 / epsilon;

      // loop over crystals and assign one hkl per crystal per reflection
      for (int i_lattice = 0; i_lattice < UB_matrices.size(); i_lattice++) {
        scitbx::mat3<double> const& A = UB_matrices[i_lattice];
        scitbx::mat3<double> const& A_inv = A.inverse();

        Graph G(reciprocal_space_points.size());

        for (std::size_t i = 0; i < reciprocal_space_points.size(); i++) {
          std::size_t i_k = i * nearest_neighbours;
          for (std::size_t i_ann = 0; i_ann < nearest_neighbours; i_ann++) {
            std::size_t i_k_plus_i_ann = i_k + i_ann;
            std::size_t j = ann.nn[i_k_plus_i_ann];
            if (boost::edge(i, j, G).second) {
              continue;
            }
            scitbx::vec3<double> d_r =
              reciprocal_space_points[i] - reciprocal_space_points[j];
            scitbx::vec3<double> h_f = A_inv * d_r;
            scitbx::vec3<double> h_ij;
            for (std::size_t ii = 0; ii < 3; ii++) {
              h_ij[ii] = scitbx::math::nearest_integer(h_f[ii]);
            }
            scitbx::vec3<double> d_h = h_f - h_ij;

            // calculate l_ij
            double exponent = 0;
            for (std::size_t ii = 0; ii < 3; ii++) {
              exponent += std::pow(
                std::max(std::abs(d_h[ii]) - epsilon, 0.) * one_over_epsilon, 2);
              exponent += std::pow(std::max(std::abs(h_ij[ii]) - delta, 0.), 2);
            }
            exponent *= -2;
            double l_ij = 1 - std::exp(exponent);

            // add the edge
            Edge e1;
            bool b1;
            boost::tie(e1, b1) = add_edge(i, j, G);

            // set the edge properties
            G[e1].h_ij = h_ij;
            G[e1].l_ij = l_ij;
            G[e1].i = i;
            G[e1].j = j;
          }
        }

        // create a graph for the MST
        Graph MST(reciprocal_space_points.size());

        // compute the minimum spanning tree
        std::vector<Edge> spanning_tree;
        kruskal_minimum_spanning_tree(G,
                                      std::back_inserter(spanning_tree),
                                      boost::weight_map(boost::get(&MyEdge::l_ij, G)));

        for (size_t i = 0; i < spanning_tree.size(); ++i) {
          // get the edge
          Edge e1 = spanning_tree[i];

          // add the edge to MST
          Edge e2;
          bool b2;
          boost::tie(e2, b2) = add_edge(source(e1, G), target(e1, G), MST);

          // copy over edge properties
          MST[e2].h_ij = G[e1].h_ij;
          MST[e2].l_ij = G[e1].l_ij;
          MST[e2].i = G[e1].i;
          MST[e2].j = G[e1].j;
        }

        // Get the connected components in case there is more than one tree
        std::vector<int> component(num_vertices(MST));
        connected_components(MST, &component[0]);

        // Record the order of edges for a depth first search so we can walk
        // the tree later
        std::vector<Edge> ordered_edges;
        record_dfs_order<Edge> vis(ordered_edges);
        depth_first_search(MST, visitor(vis));

        // Walk the tree(s), incrementing the subtree_id for each node if the
        // edge weight l_ij >= l_min, or if we start a new component
        std::size_t next_subtree = 0;
        int last_component = -1;
        af::shared<std::size_t> subtree_ids_(reciprocal_space_points.size(), 0);
        af::shared<cctbx::miller::index<> > hkl_ints_(reciprocal_space_points.size());

        for (std::vector<Edge>::iterator it = ordered_edges.begin();
             it != ordered_edges.end();
             ++it) {
          Edge e = *it;
          std::size_t i = source(e, MST);
          std::size_t j = target(e, MST);
          if (component[i] != last_component) {
            subtree_ids_[i] = next_subtree;
            next_subtree += 1;
          }
          last_component = component[i];
          cctbx::miller::index<> h_ij = MST[e].h_ij;
          if (i == MST[e].j) {
            // edge is other direction, negate h_ij
            h_ij = -h_ij;
          }
          double l_ij = MST[e].l_ij;
          hkl_ints_[j] = hkl_ints_[i] - h_ij;
          if (l_ij < l_min) {
            subtree_ids_[j] = subtree_ids_[i];
          } else {
            subtree_ids_[j] = subtree_ids_[i] + 1;
            next_subtree += 1;
          }
        }

        std::set<std::size_t> unique_subtree_ids(subtree_ids_.begin(),
                                                 subtree_ids_.end());
        std::size_t largest_subtree_id = 0;
        std::size_t largest_subtree_size = 0;
        for (std::set<std::size_t>::iterator it = unique_subtree_ids.begin();
             it != unique_subtree_ids.end();
             ++it) {
          std::size_t size = std::count(subtree_ids_.begin(), subtree_ids_.end(), *it);
          if (size > largest_subtree_size) {
            largest_subtree_size = size;
            largest_subtree_id = *it;
          }
        }

        for (std::size_t i = 0; i < reciprocal_space_points.size(); i++) {
          if (subtree_ids_[i] != largest_subtree_id) {
            continue;
          } else if (crystal_ids_[i] == -2) {
            continue;
          } else if (crystal_ids_[i] == -1) {
            miller_indices_[i] = hkl_ints_[i];
            crystal_ids_[i] = i_lattice;
          } else {
            crystal_ids_[i] = -2;
            miller_indices_[i] = cctbx::miller::index<>(0, 0, 0);
          }
        }
      }
    }

    af::shared<cctbx::miller::index<> > miller_indices() {
      return miller_indices_;
    }

    af::shared<int> crystal_ids() {
      return crystal_ids_;
    }

  private:
    af::shared<cctbx::miller::index<> > miller_indices_;
    af::shared<std::size_t> subtree_ids_;
    af::shared<int> crystal_ids_;
  };

}}  // namespace dials::algorithms

#endif
