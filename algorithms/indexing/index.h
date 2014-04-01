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
#include <scitbx/array_family/flex_types.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/math/utils.h>
#include <cctbx/miller.h>
#include <dials/array_family/scitbx_shared_and_versa.h>


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

}}

#endif
