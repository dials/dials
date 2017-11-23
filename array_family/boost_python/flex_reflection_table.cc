/*
 * flex_reflection_table.cc
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
#include <numeric>
#include <dials/array_family/boost_python/flex_table_suite.h>
#include <dials/array_family/reflection_table.h>
#include <dials/array_family/reflection.h>
#include <dials/model/data/shoebox.h>
#include <dials/model/data/observation.h>
#include <dials/algorithms/profile_model/gaussian_rs/coordinate_system.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <cctbx/miller.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;
  using flex_table_suite::flex_table_wrapper;
  using dials::model::Shoebox;
  using dials::model::Observation;
  using dials::algorithms::profile_model::gaussian_rs::CoordinateSystem;

  /**
   * Construct a reflection table from a list of observations and shoeboxes
   * @param o The observation
   * @param s The shoeboxes
   * @returns A reflection table
   */
  template <typename T>
  T* make_from_observation_and_shoebox(
      const af::const_ref<Observation> &o,
      const af::const_ref< Shoebox<> > &s) {
    // FIXME Should remove as we don't know whether intensity is summed or
    // profile fitted. In any case observation model should be deprecated.
    DIALS_ASSERT(o.size() == s.size());

    // The reflection table
    T result(o.size());
    af::shared<std::size_t>    panel  = result["panel"];
    af::shared< vec3<double> > xyzval = result["xyzobs.px.value"];
    af::shared< vec3<double> > xyzvar = result["xyzobs.px.variance"];
    af::shared<double>         iraw   = result["intensity.sum.value"];
    af::shared<double>         irawv  = result["intensity.sum.variance"];
    af::shared<int6>           bbox   = result["bbox"];
    af::shared< Shoebox<> >    sbox   = result["shoebox"];

    // Copy all the values
    for (std::size_t i = 0; i < result.nrows(); ++i) {

      // Check panel numbers
      DIALS_ASSERT(o[i].panel == s[i].panel);
      panel[i] = o[i].panel;

      // Copy observation info
      xyzval[i] = o[i].centroid.px.position;
      xyzvar[i] = o[i].centroid.px.std_err_sq;
      iraw[i]   = o[i].intensity.observed.value;
      irawv[i]  = o[i].intensity.observed.variance;

      // Copy shoebox info
      bbox[i] = s[i].bbox;
      sbox[i].panel = s[i].panel;
      sbox[i].bbox = s[i].bbox;
      sbox[i].data = s[i].data;
      sbox[i].mask = s[i].mask;
      sbox[i].background = s[i].background;
    }

    // Return the new reflection table
    return new T(result);
  }

  /**
   * A function to print some help about keys
   */
  template <typename T>
  std::string help_keys(const T &self) {
    std::string result =
      "Standard column names:\n"
      "======================\n"
      "\n"
      " Columns in the reflection table can have any name and type;\n"
      " however, it is helpful to have a set of standard data columns\n"
      " which can be used by different algorithms. These are shown below.\n"
      "\n"
      " General properties\n"
      " ------------------\n"
      "\n"
      "  flags:                  bit mask status flags\n"
      "  id:                     experiment id\n"
      "  panel:                  the detector panel index\n"
      "\n"
      " Predicted properties\n"
      " --------------------\n"
      "\n"
      "  miller_index:           miller indices\n"
      "  entering:               reflection entering/exiting\n"
      "  s1:                     the diffracted beam vector\n"
      "  xyzcal.mm:              the predicted location (mm, mm, rad)\n"
      "  xyzcal.px:              the predicted location (px, px, frame)\n"
      "  ub_matrix:              predicted crystal setting\n"
      "\n"
      " Observed properties\n"
      " -------------------\n"
      "\n"
      "  xyzobs.px.value:        centroid pixel position\n"
      "  xyzobs.px.variance:     centroid pixel variance\n"
      "  xyzobs.mm.value:        centroid millimetre position\n"
      "  xyzobs.mm.variance:     centroid millimetre variance\n"
      "  rlp:                    reciprocal lattice point\n"
      "  intensity.sum.value:    raw intensity value\n"
      "  intensity.sum.variance: raw intensity variance\n"
      "  intensity.prf.value:    profile fitted intensity value\n"
      "  intensity.prf.variance: profile fitted intensity variance\n"
      "  lp:                     LP correction (multiplicative)\n"
      "  profile.correlation:    correlation in profile fitting\n"
      "\n"
      " Shoebox properties\n"
      " ------------------\n"
      "\n"
      "  bbox:                   bounding box\n"
      "  shoebox:                shoebox data/mask/background struct\n"
      "\n"
      ;
    return result;
  }

  /**
   * Do ray intersections for all items
   */
  template <typename T>
  af::shared< vec2<double> > compute_ray_intersections(
      const T &self, const model::Detector &detector) {
    af::shared< vec2<double> > result(self.nrows());
    af::const_ref< vec3<double> > s1 =
      self.template get< vec3<double> >("s1").const_ref();
    af::const_ref< std::size_t > panel =
      self.template get< std::size_t >("panel").const_ref();
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = detector[panel[i]].get_ray_intersection(s1[i]);
    }
    return result;
  }

  /**
   * Get where the flag value is set
   */
  template <typename T>
  af::shared<bool> get_flags(const T &self, std::size_t value, bool all) {
    af::shared<bool> result(self.nrows());
    af::shared<std::size_t> flags = self.template get<std::size_t>("flags");
    DIALS_ASSERT(flags.size() == result.size());
    if (all) {
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = (flags[i] & value) == value;
      }
    } else {
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = (flags[i] & value) != 0;
      }
    }
    return result;
  }

  /**
   * Set the flags.
   */
  template <typename T>
  void set_flags_by_mask(T self, const af::const_ref<bool> mask,
      std::size_t value) {
    DIALS_ASSERT(mask.size() == self.nrows());
    af::shared<std::size_t> flags = self.template get<std::size_t>("flags");
    for (std::size_t i = 0; i < mask.size(); ++i) {
      if (mask[i]) {
        flags[i] |= value;
      }
    }
  }

  /**
   * Set the flags.
   */
  template <typename T>
  void set_flags_by_index(T self, const af::const_ref<std::size_t> index,
      std::size_t value) {
    af::shared<std::size_t> flags = self.template get<std::size_t>("flags");
    for (std::size_t i = 0; i < index.size(); ++i) {
      DIALS_ASSERT(index[i] < flags.size());
      flags[index[i]] |= value;
    }
  }

  /**
   * Unset the flags.
   */
  template <typename T>
  void unset_flags_by_mask(T self, const af::const_ref<bool> mask,
      std::size_t value) {
    DIALS_ASSERT(mask.size() == self.nrows());
    af::shared<std::size_t> flags = self.template get<std::size_t>("flags");
    for (std::size_t i = 0; i < mask.size(); ++i) {
      if (mask[i]) {
        flags[i] &= ~value;
      }
    }
  }

  /**
   * Unset the flags.
   */
  template <typename T>
  void unset_flags_by_index(T self, const af::const_ref<std::size_t> index,
      std::size_t value) {
    af::shared<std::size_t> flags = self.template get<std::size_t>("flags");
    for (std::size_t i = 0; i < index.size(); ++i) {
      DIALS_ASSERT(index[i] < flags.size());
      flags[index[i]] &= ~value;
    }
  }

  /**
   * Split the reflection table in partials.
   */
  template <typename T>
  void split_partials(T self) {

    // Check the input
    DIALS_ASSERT(self.is_consistent());
    DIALS_ASSERT(self.contains("bbox"));
    DIALS_ASSERT(self.contains("shoebox") == false);

    // Compute the number of partials
    af::const_ref<int6> bbox = self["bbox"];
    std::size_t num_full = bbox.size();
    std::size_t num_partial = 0;
    for (std::size_t i = 0; i < bbox.size(); ++i) {
      DIALS_ASSERT(bbox[i][1] > bbox[i][0]);
      DIALS_ASSERT(bbox[i][3] > bbox[i][2]);
      DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
      num_partial += (bbox[i][5] - bbox[i][4]);
    }
    DIALS_ASSERT(num_partial >= num_full);

    // If num partial is the same as num full exit early
    if (num_partial == num_full) {
      return;
    }

    // Create the new bounding boxes and indices
    af::shared<int6> bbox_new(num_partial);
    af::shared<std::size_t> indices(num_partial);
    std::size_t j = 0;
    for (std::size_t i = 0; i < num_full; ++i) {
      int6 b = bbox[i];
      for (int z = bbox[i][4]; z < bbox[i][5]; ++z) {
        DIALS_ASSERT(j < num_partial);
        bbox_new[j] = b;
        bbox_new[j][4] = z;
        bbox_new[j][5] = z + 1;
        indices[j] = i;
        j++;
      }
    }
    DIALS_ASSERT(j == num_partial);

    // Resize the reflection table
    self.resize(num_partial);

    // Reorder the reflections
    flex_table_suite::reorder(self, indices.const_ref());

    // Set the new bounding boxes
    flex_table_suite::setitem_column(self, "bbox", bbox_new.const_ref());
    flex_table_suite::setitem_column(self, "partial_id", indices.const_ref());
  }

  /**
   * Split the reflection table in partials.
   */
  template <typename T>
  void split_partials_with_shoebox(T self) {

    // Check the input
    DIALS_ASSERT(self.is_consistent());
    DIALS_ASSERT(self.contains("bbox"));
    DIALS_ASSERT(self.contains("shoebox"));
    DIALS_ASSERT(self.contains("panel"));

    // Compute the number of partials
    af::const_ref< Shoebox<> > sbox = self["shoebox"];
    af::const_ref<int6> bbox = self["bbox"];
    af::const_ref<std::size_t> panel = self["panel"];
    std::size_t num_full = bbox.size();
    std::size_t num_partial = 0;
    for (std::size_t i = 0; i < bbox.size(); ++i) {
      DIALS_ASSERT(bbox[i][1] > bbox[i][0]);
      DIALS_ASSERT(bbox[i][3] > bbox[i][2]);
      DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
      num_partial += (bbox[i][5] - bbox[i][4]);
    }
    DIALS_ASSERT(num_partial >= num_full);

    // If num partial is the same as num full exit early
    if (num_partial == num_full) {
      return;
    }

    // Create the new bounding boxes and indices
    af::shared< Shoebox<> > sbox_new(num_partial);
    af::shared<int6> bbox_new(num_partial);
    af::shared<std::size_t> indices(num_partial);
    std::size_t j = 0;
    for (std::size_t i = 0; i < num_full; ++i) {
      const Shoebox<> &s1 = sbox[i];
      int6 b = bbox[i];
      DIALS_ASSERT(s1.is_consistent());
      DIALS_ASSERT(s1.bbox[0] == b[0] && s1.bbox[1] == b[1]);
      DIALS_ASSERT(s1.bbox[2] == b[2] && s1.bbox[3] == b[3]);
      DIALS_ASSERT(s1.bbox[4] == b[4] && s1.bbox[5] == b[5]);
      std::size_t first = 0;
      for (int z = bbox[i][4]; z < bbox[i][5]; ++z) {
        DIALS_ASSERT(j < num_partial);
        bbox_new[j] = b;
        bbox_new[j][4] = z;
        bbox_new[j][5] = z + 1;
        indices[j] = i;
        Shoebox<> s2(panel[i], bbox_new[j]);
        s2.allocate();
        DIALS_ASSERT(s2.is_consistent());
        std::size_t last = first + s2.data.size();
        std::copy(
            s1.data.begin() + first,
            s1.data.begin() + last,
            s2.data.begin());
        std::copy(
            s1.mask.begin() + first,
            s1.mask.begin() + last,
            s2.mask.begin());
        std::copy(
            s1.background.begin() + first,
            s1.background.begin() + last,
            s2.background.begin());
        sbox_new[j] = s2;
        j++;
        first = last;
      }
    }
    DIALS_ASSERT(j == num_partial);

    // Resize the reflection table
    self.resize(num_partial);

    // Reorder the reflections
    flex_table_suite::reorder(self, indices.const_ref());

    // Set the new bounding boxes
    flex_table_suite::setitem_column(self, "bbox", bbox_new.const_ref());
    flex_table_suite::setitem_column(self, "shoebox", sbox_new.const_ref());
    flex_table_suite::setitem_column(self, "partial_id", indices.const_ref());
  }

  /**
   * Split the reflection table in partials and return indices.
   */
  template <typename T>
  af::shared<std::size_t> split_partial_indices(T self) {

    // Check the input
    DIALS_ASSERT(self.is_consistent());
    DIALS_ASSERT(self.contains("bbox"));

    // Compute the number of partials
    af::const_ref<int6> bbox = self["bbox"];
    std::size_t num_full = bbox.size();
    std::size_t num_partial = 0;
    for (std::size_t i = 0; i < bbox.size(); ++i) {
      DIALS_ASSERT(bbox[i][1] > bbox[i][0]);
      DIALS_ASSERT(bbox[i][3] > bbox[i][2]);
      DIALS_ASSERT(bbox[i][5] > bbox[i][4]);
      num_partial += (bbox[i][5] - bbox[i][4]);
    }
    DIALS_ASSERT(num_partial >= num_full);

    // If num partial is the same as num full exit early
    af::shared<std::size_t> indices(num_partial);
    if (num_partial == num_full) {
      for (std::size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
      }
      return indices;
    }

    // Create the new bounding boxes and indices
    std::size_t j = 0;
    for (std::size_t i = 0; i < num_full; ++i) {
      for (int z = bbox[i][4]; z < bbox[i][5]; ++z) {
        DIALS_ASSERT(j < num_partial);
        indices[j] = i;
        j++;
      }
    }
    DIALS_ASSERT(j == num_partial);
    return indices;
  }

  /**
   * Split reflection table by experiment id
   */
  template <typename T>
  boost::python::list split_by_experiment_id(T self) {
    DIALS_ASSERT(self.contains("id"));

    // Get the id array
    af::const_ref<int> id = self["id"];

    // Get the number of experiments
    std::size_t num_expr = 0;
    for (std::size_t i = 0; i < id.size(); ++i) {
      DIALS_ASSERT(id[i] >= 0);
      if (id[i] >= num_expr) num_expr = id[i] + 1;
    }

    // Get the number of each
    std::vector<std::size_t> num(num_expr, 0);
    for (std::size_t i = 0; i < id.size(); ++i) {
      num[id[i]]++;
    }

    // Compute the indices
    std::vector<std::size_t> indices(id.size());
    std::vector<std::size_t> offset(1, 0);
    std::partial_sum(num.begin(), num.end(), std::back_inserter(offset));
    num.assign(num.size(), 0);
    for (std::size_t i = 0; i < indices.size(); ++i) {
      std::size_t j = id[i];
      DIALS_ASSERT(j < offset.size() - 1);
      std::size_t off1 = offset[j];
      std::size_t off2 = offset[j+1];
      DIALS_ASSERT(off2 > off1);
      DIALS_ASSERT(off2 <= indices.size());
      std::size_t k = off1 + num[j];
      DIALS_ASSERT(j < off2);
      num[j]++;
      indices[k] = i;
    }

    // For each experiment if select the reflections in the list
    boost::python::list result;
    for (std::size_t i = 0; i < offset.size()-1; ++i) {
      std::size_t off1 = offset[i];
      std::size_t off2 = offset[i+1];
      DIALS_ASSERT(off2 >= off1);
      std::size_t off = off1;
      std::size_t num = off2 - off1;
      if (num > 0) {
        DIALS_ASSERT(off + num <= indices.size());
        result.append(flex_table_suite::select_rows_index(
              self, const_ref<std::size_t>(&indices[off], num)));
      }
    }

    // Return the result
    return result;
  }

  /**
   * Split reflection table by experiment id
   */
  template <typename T>
  boost::python::list split_indices_by_experiment_id(
      T self, std::size_t num_expr) {
    DIALS_ASSERT(self.size() > 0);
    DIALS_ASSERT(num_expr > 0);
    DIALS_ASSERT(self.contains("id"));

    // Get the id array
    af::const_ref<int> id = self["id"];

    // Get the number of each
    std::vector<std::size_t> num(num_expr, 0);
    for (std::size_t i = 0; i < id.size(); ++i) {
      DIALS_ASSERT(id[i] >= 0);
      DIALS_ASSERT(id[i] < num_expr);
      num[id[i]]++;
    }

    // Compute the indices
    std::vector<std::size_t> indices(id.size());
    std::vector<std::size_t> offset(1, 0);
    std::partial_sum(num.begin(), num.end(), std::back_inserter(offset));
    num.assign(num.size(), 0);
    for (std::size_t i = 0; i < indices.size(); ++i) {
      std::size_t exp_id = id[i];
      DIALS_ASSERT(exp_id < offset.size() - 1);
      std::size_t off1 = offset[exp_id];
      std::size_t off2 = offset[exp_id+1];
      DIALS_ASSERT(off2 > off1);
      DIALS_ASSERT(off2 <= indices.size());
      std::size_t k = off1 + num[exp_id];
      DIALS_ASSERT(k < off2);
      num[exp_id]++;
      indices[k] = i;
    }

    // For each experiment if select the reflections in the list
    boost::python::list result;
    for (std::size_t i = 0; i < offset.size()-1; ++i) {
      std::size_t off1 = offset[i];
      std::size_t off2 = offset[i+1];
      DIALS_ASSERT(off2 >= off1);
      std::size_t off = off1;
      std::size_t num = off2 - off1;
      if (num > 0) {
        DIALS_ASSERT(off + num <= indices.size());
        result.append(af::shared<std::size_t>(&indices[off], &indices[off]+num));
      } else {
        result.append(af::shared<std::size_t>());
      }
    }

    // Return the result
    return result;
  }

  /**
   * Compute phi range of reflection
   */
  template <typename T>
  af::shared< vec2<double> > compute_phi_range(
      T self,
      vec3<double> m2,
      vec3<double> s0,
      double sigma_m,
      double n_sigma) {

    // Check contents
    DIALS_ASSERT(sigma_m > 0);
    DIALS_ASSERT(n_sigma > 0);
    DIALS_ASSERT(self.contains("s1"));
    DIALS_ASSERT(self.contains("xyzcal.mm"));

    af::const_ref< vec3<double> > s1 = self["s1"];
    af::const_ref< vec3<double> > xyz = self["xyzcal.mm"];
    af::shared< vec2<double> > result(self.size());
    double delta_m = sigma_m * n_sigma;

    for (std::size_t i = 0; i < self.size(); ++i) {

      // Create the coordinate system for the reflection
      CoordinateSystem xcs(m2, s0, s1[i], xyz[i][2]);

      /// Calculate the rotation angles at the following XDS
      // e3 coordinates: -delta_m, +delta_m
      double phi1 = xcs.to_rotation_angle_fast(-delta_m);
      double phi2 = xcs.to_rotation_angle_fast(+delta_m);
      if (phi2 < phi1) {
        std::swap(phi1,phi2);
      }
      result[i] = vec2<double>(phi1,phi2);
    }
    return result;
  }

  /**
   * A visitor to convert an item to an object
   */
  struct item_to_object_visitor : public boost::static_visitor<object> {
    template <typename T>
    object operator () (T &data) {
      return object(data);
    }
  };

  /**
   * Get an item from the reflection
   * @param self The reflection
   * @param name The item name
   */
  boost::python::object Reflection_get(const Reflection &self, std::string name) {
    Reflection::mapped_type item = self[name];
    item_to_object_visitor visitor;
    return item.apply_visitor(visitor);
  }

  /**
   * Set a bool item in the reflection
   * @param self The reflection
   * @param name The name of the item
   * @param item The item
   */
  void Reflection_set_bool(Reflection &self, std::string name, bool item) {
    self[name] = Reflection::data_type(item);
  }

  /**
   * Set a int item in the reflection
   * @param self The reflection
   * @param name The name of the item
   * @param item The item
   */
  void Reflection_set_int(Reflection &self, std::string name, int item) {
    self[name] = Reflection::data_type(item);
  }

  /**
   * Set a std::size_t item in the reflection
   * @param self The reflection
   * @param name The name of the item
   * @param item The item
   */
  void Reflection_set_size_t(Reflection &self, std::string name, std::size_t item) {
    self[name] = Reflection::data_type(item);
  }

  /**
   * Set a double item in the reflection
   * @param self The reflection
   * @param name The name of the item
   * @param item The item
   */
  void Reflection_set_double(Reflection &self, std::string name, double item) {
    self[name] = Reflection::data_type(item);
  }

  /**
   * Set a string item in the reflection
   * @param self The reflection
   * @param name The name of the item
   * @param item The item
   */
  void Reflection_set_string(Reflection &self, std::string name, std::string item) {
    self[name] = Reflection::data_type(item);
  }

  /**
   * Set a vec2<double> item in the reflection
   * @param self The reflection
   * @param name The name of the item
   * @param item The item
   */
  void Reflection_set_vec2_double(Reflection &self, std::string name, vec2<double> item) {
    self[name] = Reflection::data_type(item);
  }

  /**
   * Set a vec3<double> item in the reflection
   * @param self The reflection
   * @param name The name of the item
   * @param item The item
   */
  void Reflection_set_vec3_double(Reflection &self, std::string name, vec3<double> item) {
    self[name] = Reflection::data_type(item);
  }

  /**
   * Set a mat3<double> item in the reflection
   * @param self The reflection
   * @param name The name of the item
   * @param item The item
   */
  void Reflection_set_mat3_double(Reflection &self, std::string name, mat3<double> item) {
    self[name] = Reflection::data_type(item);
  }

  /**
   * Set a int6 item in the reflection
   * @param self The reflection
   * @param name The name of the item
   * @param item The item
   */
  void Reflection_set_int6(Reflection &self, std::string name, int6 item) {
    self[name] = Reflection::data_type(item);
  }

  /**
   * Set a miller_index item in the reflection
   * @param self The reflection
   * @param name The name of the item
   * @param item The item
   */
  void Reflection_set_miller_index(Reflection &self, std::string name, cctbx::miller::index<> item) {
    self[name] = Reflection::data_type(item);
  }

  /**
   * Set a Shoebox<> item in the reflection
   * @param self The reflection
   * @param name The name of the item
   * @param item The item
   */
  void Reflection_set_shoebox(Reflection &self, std::string name, Shoebox<> item) {
    self[name] = Reflection::data_type(item);
  }


  /**
   * Convert reflection table to list of reflections
   * @param self The reflection table
   * @returns The list of reflections
   */
  boost::python::list reflection_table_to_list_of_reflections(
      reflection_table self) {
    af::shared<Reflection> array = reflection_table_to_array(self);
    boost::python::list result;
    for (std::size_t i = 0; i < array.size(); ++i) {
      result.append(array[i]);
    }
    return result;
  }


  /**
   * Struct to facilitate wrapping reflection table type
   */
  template <typename T>
  struct flex_reflection_table_wrapper : public flex_table_wrapper<T> {

    typedef flex_table_wrapper<T> base_type;
    typedef typename base_type::flex_types flex_types;
    typedef typename base_type::flex_table_type flex_table_type;
    typedef typename base_type::class_type class_type;

    /**
     * Wrap the reflection table class
     */
    static
    class_type wrap(const char *name) {

      // Wrap with flex table bindings
      class_type result = base_type::wrap(name);

      // Add functions
      result
        .def("__init__", make_constructor(
          &make_from_observation_and_shoebox<flex_table_type>))
        .def("help_keys",
          &help_keys<flex_table_type>)
        .def("compute_ray_intersections",
          &compute_ray_intersections<flex_table_type>)
        .def("get_flags",
          &get_flags<flex_table_type>, (
            boost::python::arg("value"),
            boost::python::arg("all") = true))
        .def("set_flags",
          &set_flags_by_mask<flex_table_type>)
        .def("set_flags",
          &set_flags_by_index<flex_table_type>)
        .def("unset_flags",
          &unset_flags_by_mask<flex_table_type>)
        .def("unset_flags",
          &unset_flags_by_index<flex_table_type>)
        .def("split_partials",
          &split_partials<flex_table_type>)
        .def("split_partials_with_shoebox",
          &split_partials_with_shoebox<flex_table_type>)
        .def("split_partial_indices",
          &split_partial_indices<flex_table_type>)
        .def("split_by_experiment_id",
          &split_by_experiment_id<flex_table_type>)
        .def("split_indices_by_experiment_id",
          &split_indices_by_experiment_id<flex_table_type>)
        .def("compute_phi_range",
          &compute_phi_range<flex_table_type>)
        ;

      // Create the flags enum in the reflection table scope
      scope in_table = result;
      enum_<Flags>("flags")
        .value("predicted", Predicted)
        .value("observed", Observed)
        .value("indexed", Indexed)
        .value("used_in_refinement", UsedInRefinement)
        .value("strong", Strong)
        .value("reference_spot", ReferenceSpot)
        .value("dont_integrate", DontIntegrate)
        .value("integrated_sum", IntegratedSum)
        .value("integrated_prf", IntegratedPrf)
        .value("integrated", Integrated)
        .value("overloaded", Overloaded)
        .value("overlapped_bg", OverlappedBg)
        .value("overlapped_fg", OverlappedFg)
        .value("in_powder_ring", InPowderRing)
        .value("foreground_includes_bad_pixels", ForegroundIncludesBadPixels)
        .value("background_includes_bad_pixels", BackgroundIncludesBadPixels)
        .value("includes_bad_pixels", IncludesBadPixels)
        .value("bad_shoebox", BadShoebox)
        .value("bad_spot", BadSpot)
        .value("used_in_modelling", UsedInModelling)
        .value("centroid_outlier", CentroidOutlier)
        .value("failed_during_background_modelling", FailedDuringBackgroundModelling)
        .value("failed_during_summation", FailedDuringSummation)
        .value("failed_during_profile_fitting", FailedDuringProfileFitting)
        .value("bad_reference", BadReference)
        ;

      // return the wrapped class
      return result;
    }
  };

  void export_flex_reflection_table() {

    // Set the do
    docstring_options local_docstring_options;
    local_docstring_options.enable_user_defined();
    local_docstring_options.enable_py_signatures();
    local_docstring_options.disable_cpp_signatures();

    // Define all the types we want to support in the table
    typedef reflection_table::mapped_type flex_types;

    // Export the reflection table
    flex_reflection_table_wrapper<flex_types>::wrap("reflection_table");

    // Export the reflection object
    class_<Reflection>("Reflection")
      .def("get",
          &Reflection_get)
      .def("set_bool",
          &Reflection_set_bool)
      .def("set_int",
          &Reflection_set_int)
      .def("set_size_t",
          &Reflection_set_size_t)
      .def("set_double",
          &Reflection_set_double)
      .def("set_string",
          &Reflection_set_string)
      .def("set_vec2_double",
          &Reflection_set_vec2_double)
      .def("set_vec3_double",
          &Reflection_set_vec3_double)
      .def("set_mat3_double",
          &Reflection_set_mat3_double)
      .def("set_int6",
          &Reflection_set_int6)
      .def("set_miller_index",
          &Reflection_set_miller_index)
      .def("set_shoebox",
          &Reflection_set_shoebox)
      ;

    // Helper function
    def("reflection_table_to_list_of_reflections", &reflection_table_to_list_of_reflections);
  }

}}} // namespace dials::af::boost_python
