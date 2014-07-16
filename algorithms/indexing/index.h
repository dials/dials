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
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <annlib_adaptbx/ann_adaptor.h>


namespace dials { namespace algorithms {

  class AssignIndices {

  public:
    AssignIndices(
      af::const_ref<scitbx::vec3<double> > const & reciprocal_space_points,
      af::const_ref<scitbx::mat3<double> > const & UB_matrices,
      double tolerance=0.3
      )
      : miller_indices_(
          reciprocal_space_points.size(), cctbx::miller::index<>(0,0,0)),
        crystal_ids_(reciprocal_space_points.size(), -1),
        n_rejects_(0) {

      typedef std::pair<const cctbx::miller::index<>, const std::size_t> pair_t;
      typedef std::multimap<cctbx::miller::index<>, std::size_t> map_t;

      map_t hkl_to_rlp_map;

      std::vector<af::shared<cctbx::miller::index<> > > hkl_ints;
      std::vector<af::shared<double> > lengths_sq;

      // loop over crystals and assign one hkl per crystal per reflection
      for (int i_lattice=0; i_lattice<UB_matrices.size(); i_lattice++) {
        scitbx::mat3<double> A = UB_matrices[i_lattice];
        scitbx::mat3<double> A_inv = A.inverse();
        af::shared<cctbx::miller::index<> > hkl_ints_(
          af::reserve(reciprocal_space_points.size()));
        af::shared<double> lengths_sq_(
          af::reserve(reciprocal_space_points.size()));
        for (int i_ref=0; i_ref<reciprocal_space_points.size(); i_ref++) {
          scitbx::vec3<double> rlp = reciprocal_space_points[i_ref];
          scitbx::vec3<double> hkl_f = A_inv * rlp;
          cctbx::miller::index<> hkl_i;
          for (std::size_t j=0; j<3; j++) {
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
      for (int i_ref=0; i_ref<reciprocal_space_points.size(); i_ref++) {
        af::shared<double> n;
        af::shared<cctbx::miller::index<> > potential_hkls;
        for (int i_lattice=0; i_lattice<UB_matrices.size(); i_lattice++) {
          n.push_back(lengths_sq[i_lattice][i_ref]);
          potential_hkls.push_back(hkl_ints[i_lattice][i_ref]);
        }
        int i_best_lattice = af::min_index(n.const_ref());
        if (n[i_best_lattice] > tolerance_sq) {
          n_rejects_++;
          continue;
        }
        cctbx::miller::index<> hkl = potential_hkls[i_best_lattice];
        miller_indices_[i_ref] = hkl;
        hkl_to_rlp_map.insert(pair_t(hkl, i_ref));
        crystal_ids_[i_ref] = i_best_lattice;
      }

      cctbx::miller::index<> curr_hkl(0,0,0);
      std::vector<std::size_t> i_same_hkl;
      // if more than one spot can be assigned the same miller index then
      // choose the closest one
      for (map_t::iterator it=hkl_to_rlp_map.begin();
           it != hkl_to_rlp_map.end(); it++) {
        if (it->first == cctbx::miller::index<>(0,0,0)) { continue; }
        if (it->first != curr_hkl) {
          if (i_same_hkl.size() > 1) {
            for (int i=0; i<i_same_hkl.size(); i++) {
              const std::size_t i_ref = i_same_hkl[i];
              for (int j=i+1; j<i_same_hkl.size(); j++) {
                const std::size_t j_ref = i_same_hkl[j];
                int crystal_i = crystal_ids_[i_ref];
                int crystal_j = crystal_ids_[j_ref];
                //DIALS_ASSERT(hkl_ints[crystal_i][i_ref] == hkl_ints[crystal_j][j_ref]);
                if (crystal_i != crystal_j) { continue; }
                else if (crystal_i == -1) { continue; }
                if (lengths_sq[crystal_j][j_ref] < lengths_sq[crystal_i][i_ref]) {
                  miller_indices_[i_ref] = cctbx::miller::index<>(0,0,0);
                  crystal_ids_[i_ref] = -1;
                }
                else {
                  miller_indices_[j_ref] = cctbx::miller::index<>(0,0,0);
                  crystal_ids_[j_ref] = -1;
                }
              }
            }
          }
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

    std::size_t n_rejects() {
      return n_rejects_;
    }

  private:
    af::shared<cctbx::miller::index<> > miller_indices_;
    af::shared<int> crystal_ids_;
    std::size_t n_rejects_;
  };


  template <typename Edge>
  struct record_dfs_order : public boost::default_dfs_visitor
  {
    record_dfs_order(std::vector<Edge> &e)
      : edges(e) { }

    template <typename Graph>
    void tree_edge(Edge e, const Graph& G) const {
      edges.push_back(e);
    }
    std::vector<Edge>& edges;
  };

  class AssignIndicesLocal {

  public:
    AssignIndicesLocal(
      af::const_ref<scitbx::vec3<double> > const & reciprocal_space_points,
      af::const_ref<scitbx::mat3<double> > const & UB_matrices,
      const double epsilon=0.05,
      const double delta=5,
      const double l_min=0.8
      )
      : miller_indices_(
          reciprocal_space_points.size(), cctbx::miller::index<>(0,0,0)),
        subtree_ids_(reciprocal_space_points.size(), 0),
        crystal_ids_(reciprocal_space_points.size(), -1),
        n_rejects_(0) {

      // Currently only support 1 lattice
      //DIALS_ASSERT(UB_matrices.size() == 1);

      using annlib_adaptbx::AnnAdaptor;
      using namespace boost;
      typedef property<edge_weight_t, double> EdgeWeightProperty;
      typedef adjacency_list<vecS, vecS, undirectedS, no_property,
        EdgeWeightProperty > Graph;

      typedef graph_traits<Graph>::vertex_descriptor Vertex;
      typedef graph_traits<Graph>::edge_descriptor Edge;

      typedef boost::property_map< Graph, boost::vertex_index_t>::type VertexIndexMap;
      typedef boost::property_map< Graph, boost::edge_weight_t>::type WeightMap;

      // convert into a single array for input to AnnAdaptor
      // based on flex.vec_3.as_double()
      // scitbx/array_family/boost_python/flex_vec3_double.cpp
      af::shared<double> rlps_double(
        reciprocal_space_points.size()*3, af::init_functor_null<double>());
      double* r = rlps_double.begin();
      for(std::size_t i=0;i<reciprocal_space_points.size();i++) {
        for(std::size_t j=0;j<3;j++) {
          *r++ = reciprocal_space_points[i][j];
        }
      }

      int nn_window = 10;
      AnnAdaptor ann = AnnAdaptor(rlps_double, 3, nn_window);
      ann.query(rlps_double);

      scitbx::mat3<double> const &A = UB_matrices[0];
      scitbx::mat3<double> const &A_inv = A.inverse();

      Graph G(reciprocal_space_points.size());

      typedef std::pair<const int, const int> pair_t;

      std::map<pair_t, cctbx::miller::index<> > edge_to_h_ij;
      std::map<pair_t, double> edge_to_l_ij;

      WeightMap weight = get(edge_weight, G);

      double sum_l_ij = 0;

      for (std::size_t i=0; i < reciprocal_space_points.size(); i++) {
        std::size_t i_k = i * nn_window;
        for (std::size_t i_ann=0; i_ann < nn_window; i_ann++) {
          std::size_t i_k_plus_i_ann = i_k + i_ann;
          std::size_t j = ann.nn[i_k_plus_i_ann];
          if (boost::edge(i, j, G).second) {
            continue;
          }
          scitbx::vec3<double>
          d_r = reciprocal_space_points[i] - reciprocal_space_points[j];
          scitbx::vec3<double> h_f = A_inv * d_r;
          scitbx::vec3<double> h_ij;
          for (std::size_t ii=0; ii<3; ii++) {
            h_ij[ii] = scitbx::math::nearest_integer(h_f[ii]);
          }
          scitbx::vec3<double> d_h = h_f - h_ij;

          double exponent = 0;
          for (std::size_t ii=0; ii<3; ii++) {
            exponent += std::pow(std::max(std::abs(d_h[ii]) - epsilon, 0.)/epsilon, 2);
            exponent += std::pow(std::max(std::abs(h_ij[ii]) - delta, 0.),2);
          }
          exponent *= -2;
          double l_ij = 1 - std::exp(exponent);
          sum_l_ij += l_ij;
          std::pair<Edge, bool> myPair = add_edge(i, j, G);
          weight[myPair.first] = l_ij;

          edge_to_h_ij.insert(
            std::pair<pair_t, cctbx::miller::index<> >(pair_t(i, j), h_ij));
          edge_to_h_ij.insert(
            std::pair<pair_t, cctbx::miller::index<> >(pair_t(j, i), -h_ij));
          edge_to_l_ij.insert(std::pair<pair_t, double>(pair_t(i, j), l_ij));
          edge_to_l_ij.insert(std::pair<pair_t, double>(pair_t(j, i), l_ij));
        }
      }

      std::vector< Edge > spanning_tree;
      kruskal_minimum_spanning_tree(G, std::back_inserter(spanning_tree));

      //create a graph for the MST
      Graph MST(reciprocal_space_points.size());

      WeightMap weightMST = get(edge_weight, MST);

      for(size_t i = 0; i < spanning_tree.size(); ++i){
        //get the edge
        Edge e = spanning_tree[i];
        //get the vertices
        add_edge(source(e, G), target(e, G), MST);
        weightMST[e] = weight[e];
      }

      std::vector<Edge> ordered_edges;

      record_dfs_order<Edge> vis(ordered_edges);
      depth_first_search(MST, visitor(vis));

      for (std::vector<Edge>::iterator it = ordered_edges.begin() ;
           it != ordered_edges.end(); ++it) {
        Edge e = *it;
        std::size_t i = source(e, MST);
        std::size_t j = target(e, MST);
        cctbx::miller::index<> h_ij = edge_to_h_ij[pair_t(i, j)];
        double l_ij = edge_to_l_ij[pair_t(i, j)];
        miller_indices_[j] = miller_indices_[i] - h_ij;
        if (l_ij < l_min) {
          subtree_ids_[j] = subtree_ids_[i];
        }
        else {
          subtree_ids_[j] = subtree_ids_[i] + 1;
        }
      }

    }

    af::shared<cctbx::miller::index<> > miller_indices() {
      return miller_indices_;
    }

    af::shared<std::size_t> subtree_ids() {
      return subtree_ids_;
    }

    af::shared<int> crystal_ids() {
      return crystal_ids_;
    }

    std::size_t n_rejects() {
      return n_rejects_;
    }

  private:
    af::shared<cctbx::miller::index<> > miller_indices_;
    af::shared<std::size_t> subtree_ids_;
    af::shared<int> crystal_ids_;
    std::size_t n_rejects_;
  };


}}

#endif
