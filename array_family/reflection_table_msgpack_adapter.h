/*
 * reflection_table_msgpack_adapter.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ARRAY_FAMILY_REFLECTION_TABLE_MSGPACK_ADAPTER_H
#define DIALS_ARRAY_FAMILY_REFLECTION_TABLE_MSGPACK_ADAPTER_H

#include <scitbx/array_family/shared.h>
#include <dials/array_family/reflection_table.h>
#include <msgpack.hpp>

namespace msgpack {
MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS) {
namespace adaptor {

  /**
   * Helper struct to give a column type a name
   */
  template <typename T>
  struct column_type {
  };

  /**
   * Helper struct for bool column type
   */
  template <>
  struct column_type <bool> {
    static
    std::string name() {
      return "bool";
    }
  };

  /**
   * Helper struct for int column type
   */
  template <>
  struct column_type <int> {
    static
    std::string name() {
      return "int";
    }
  };

  /**
   * Helper struct for std::size_t column type
   */
  template <>
  struct column_type <std::size_t> {
    static
    std::string name() {
      return "std::size_t";
    }
  };

  /**
   * Helper struct for double column type
   */
  template <>
  struct column_type <double> {
    static
    std::string name() {
      return "double";
    }
  };

  /**
   * Helper struct for std::string column type
   */
  template <>
  struct column_type <std::string> {
    static
    std::string name() {
      return "std::string";
    }
  };

  /**
   * Helper struct for vec2<double> column type
   */
  template <>
  struct column_type < scitbx::vec2<double> > {
    static
    std::string name() {
      return "vec2<double>";
    }
  };

  /**
   * Helper struct for vec3<double> column type
   */
  template <>
  struct column_type < scitbx::vec3<double> > {
    static
    std::string name() {
      return "vec3<double>";
    }
  };

  /**
   * Helper struct for mat3<double> column type
   */
  template <>
  struct column_type < scitbx::mat3<double> > {
    static
    std::string name() {
      return "mat3<double>";
    }
  };

  /**
   * Helper struct for int6 column type
   */
  template <>
  struct column_type < scitbx::af::tiny<int,6> > {
    static
    std::string name() {
      return "int6";
    }
  };

  /**
   * Helper struct for miller index column type
   */
  template <>
  struct column_type < cctbx::miller::index<> > {
    static
    std::string name() {
      return "cctbx::miller::index<>";
    }
  };

  /**
   * Helper struct for Shoebox column type
   */
  template <>
  struct column_type < dials::af::Shoebox<> > {
    static
    std::string name() {
      return "Shoebox<>";
    }
  };

  /**
   * Pack a const_ref into a msgpack array
   */
  template <typename T>
  struct pack< scitbx::af::const_ref<T> > {
    template <typename Stream>
    msgpack::packer<Stream>& operator()(
        msgpack::packer<Stream>& o,
        const scitbx::af::const_ref<T>& v) const {
      typedef typename scitbx::af::const_ref<T>::const_iterator iterator;
      o.pack_array(v.size());
      for (iterator it = v.begin(); it != v.end(); ++it) {
        o.pack(*it);
      }
      return o;
    }
  };

  /**
   * Pack a shared into a msgpack array
   */
  template <typename T>
  struct pack< scitbx::af::shared<T> > {
    template <typename Stream>
    msgpack::packer<Stream>& operator()(
        msgpack::packer<Stream>& o,
        const scitbx::af::shared<T>& v) const {
      return o.pack(v.const_ref());
    }
  };

  /**
   * Pack a versa into a msgpack array and preserve accessor size
   */
  template <typename T, typename Accessor>
  struct pack< scitbx::af::versa<T, Accessor> > {
    template <typename Stream>
    msgpack::packer<Stream>& operator()(
        msgpack::packer<Stream>& o,
        const scitbx::af::versa<T, Accessor>& v) const {
      o.pack_array(2);
      o.pack(v.accessor().const_ref());
      o.pack(scitbx::af::const_ref<T>(&v[0], v.size()));
      return o;
    }
  };

  /**
   * Pack an af::tiny into a msgpack array
   */
  template <typename T, std::size_t N>
  struct pack< scitbx::af::tiny<T,N> > {
    template <typename Stream>
    msgpack::packer<Stream>& operator()(
        msgpack::packer<Stream>& o,
        const scitbx::af::tiny<T,N>& v) const {
      return o.pack(v.const_ref());
    }
  };

  /**
   * Pack a vec2 into a msgpack array
   */
  template <typename T>
  struct pack< scitbx::vec2<T> > {
    template <typename Stream>
    msgpack::packer<Stream>& operator()(
        msgpack::packer<Stream>& o,
        const scitbx::vec2<T>& v) const {
      return o.pack(v.const_ref());
    }
  };

  /**
   * Pack a vec3 into a msgpack array
   */
  template <typename T>
  struct pack< scitbx::vec3<T> > {
    template <typename Stream>
    msgpack::packer<Stream>& operator()(
        msgpack::packer<Stream>& o,
        const scitbx::vec3<T>& v) const {
      return o.pack(v.const_ref());
    }
  };

  /**
   * Pack a mat3 into a msgpack array
   */
  template <typename T>
  struct pack< scitbx::mat3<T> > {
    template <typename Stream>
    msgpack::packer<Stream>& operator()(
        msgpack::packer<Stream>& o,
        const scitbx::mat3<T>& v) const {
      return o.pack(v.const_ref());
    }
  };

  /**
   * Pack a cctbx::miller::index<> into a msgpack array
   */
  template <typename T>
  struct pack< cctbx::miller::index<T> > {
    template <typename Stream>
    msgpack::packer<Stream>& operator()(
        msgpack::packer<Stream>& o,
        const cctbx::miller::index<T>& v) const {
      return o.pack(v.const_ref());
    }
  };

  /**
   * Pack a Shoebox<> structure into an array like:
   * [ panel, bbox, data, mask, background ]
   */
  template <typename T>
  struct pack< dials::af::Shoebox<T> > {
    template <typename Stream>
    msgpack::packer<Stream>& operator()(
        msgpack::packer<Stream>& o,
        const dials::af::Shoebox<T>& v) const {
      o.pack_array(5);
      o.pack(v.panel);
      o.pack(v.bbox);
      o.pack(v.data);
      o.pack(v.mask);
      o.pack(v.background);
      return o;
    }
  };

  /**
   * A visitor to help with packing a column variant type.
   * Pack the column into an array like: [ name, data ]
   */
  template <typename Stream>
  struct packer_visitor : boost::static_visitor<void> {
    packer<Stream>& o_;
    packer_visitor(packer<Stream>& o)
      : o_(o) {}
    template <typename T>
    void operator()(T const& value) const {
      o_.pack_array(2);
      o_.pack(column_type<typename T::value_type>::name());
      o_.pack(value);
    }
  };

  /**
   * Pack the reflection table into an array with a map like:
   * [
   *  "dials::af::reflection_table",
   *  VERSION,
   *  N_ROWS,
   *  DATA
   *  ]
   *
   * The first entry identifies the data as a reflection table.
   * The second entry gives the version number in case this changes
   * The third entry gives the expected number of rows in the table
   * The fourth entry is a map with key value pairs corresponding to the names
   * and arrays of the column data.
   *
   */
  template <>
  struct pack< dials::af::reflection_table > {
    template <typename Stream>
    msgpack::packer<Stream>& operator()(
        msgpack::packer<Stream>& o,
        const dials::af::reflection_table& v) const {
      typedef typename dials::af::reflection_table::const_iterator iterator;
      o.pack_array(4);
      o.pack("dials::af::reflection_table");
      o.pack(1);
      o.pack(v.nrows());
      o.pack_map(v.ncols());
      for (iterator it = v.begin(); it != v.end(); ++it) {
        o.pack(it->first);
        boost::apply_visitor(packer_visitor<Stream>(o), it->second);
      }
      return o;
    }
  };

  /**
   * Convert a msgpack array into a fixed size scitbx::af::ref
   */
  template<typename T>
  struct convert< scitbx::af::ref<T> > {
    msgpack::object const& operator()(
        msgpack::object const& o,
        scitbx::af::ref<T>& v) const {
      typedef typename scitbx::af::ref<T>::iterator iterator;

      // Ensure the type is an array
      if (o.type != msgpack::type::ARRAY) {
        throw DIALS_ERROR("msgpack type is not an array");
      }

      // Ensure it is of the correct size
      if (o.via.array.size != v.size()) {
        throw DIALS_ERROR("msgpack array does not have correct dimensions");
      }

      // Convert the values in the array
      if (o.via.array.size > 0) {
        msgpack::object* first = o.via.array.ptr;
        msgpack::object* last = first + o.via.array.size;
        iterator out = v.begin();
        for (msgpack::object *it = first; it != last; ++it) {
          it->convert(*out++);
        }
      }
      return o;
    }
  };

  /**
   * Convert a msgpack array into a variable size scitbx::af::shared
   */
  template<typename T>
  struct convert< scitbx::af::shared<T> > {
    msgpack::object const& operator()(
        msgpack::object const& o,
        scitbx::af::shared<T>& v) const {
      typedef typename scitbx::af::shared<T>::iterator iterator;

      // Ensure the type is an array
      if (o.type != msgpack::type::ARRAY) {
        throw DIALS_ERROR("msgpack type is not an array");
      }

      // Resize the array
      v.resize(o.via.array.size);

      // Convert the values in the array
      if (o.via.array.size > 0) {
        msgpack::object* first = o.via.array.ptr;
        msgpack::object* last = first + o.via.array.size;
        iterator out = v.begin();
        for (msgpack::object *it = first; it != last; ++it) {
          it->convert(*out++);
        }
      }
      return o;
    }
  };

  /**
   * Convert a msgpack array into a scitbx::af::versa with accessor
   */
  template<typename T, typename Accessor>
  struct convert< scitbx::af::versa<T, Accessor> > {
    msgpack::object const& operator()(
        msgpack::object const& o,
        scitbx::af::versa<T, Accessor>& v) const {

      // Ensure type is an array
      if (o.type != msgpack::type::ARRAY) {
        throw DIALS_ERROR("msgpack type is not an array");
      }

      // Ensure that we have an accessor and data element
      if (o.via.array.size != 2) {
        throw DIALS_ERROR("msgpack array does not have correct dimensions");
      }

      // Read the accessor element
      Accessor grid;
      scitbx::af::ref<std::size_t> grid_ref = grid.ref();
      o.via.array.ptr[0].convert(grid_ref);

      // Resize the versa
      v.resize(grid);

      // Read the data
      scitbx::af::ref<T> data_ref(&v[0], v.size());
      o.via.array.ptr[1].convert(data_ref);
      return o;
    }
  };

  /**
   * Convert a msgpack array to an int6
   */
  template <typename T, std::size_t N>
  struct convert< scitbx::af::tiny<T,N> > {
    msgpack::object const& operator()(
        msgpack::object const& o,
        scitbx::af::tiny<T,N>& v) const {
      scitbx::af::ref<T> r = v.ref();
      o.convert(r);
      return o;
    }
  };

  /**
   * Convert a msgpack array to a vec2
   */
  template <typename T>
  struct convert< scitbx::vec2<T> > {
    msgpack::object const& operator()(
        msgpack::object const& o,
        scitbx::vec2<T>& v) const {
      scitbx::af::ref<T> r = v.ref();
      o.convert(r);
      return o;
    }
  };

  /**
   * Convert a msgpack array to a vec3
   */
  template <typename T>
  struct convert< scitbx::vec3<T> > {
    msgpack::object const& operator()(
        msgpack::object const & o,
        scitbx::vec3<T>& v) const {
      scitbx::af::ref<T> r = v.ref();
      o.convert(r);
      return o;
    }
  };

  /**
   * Convert a msgpack array to a mat3
   */
  template <typename T>
  struct convert< scitbx::mat3<T> > {
    msgpack::object const& operator()(
        msgpack::object const& o,
        scitbx::mat3<T>& v) const {
      scitbx::af::ref<T> r = v.ref();
      o.convert(r);
      return o;
    }
  };

  /**
   * Convert a msgpack array to a cctbx::miller::index
   */
  template <typename T>
  struct convert< cctbx::miller::index<T> > {
    msgpack::object const& operator()(
        msgpack::object const& o,
        cctbx::miller::index<T>& v) const {
      scitbx::af::ref<T> r = v.ref();
      o.convert(r);
      return o;
    }
  };

  /**
   * Convert a msgpack array to a Shoebox structure.
   * The msgpack array will have a structure like:
   * [ panel, bbox, data, mask, background ]
   */
  template<typename T>
  struct convert< dials::af::Shoebox<T> > {
    msgpack::object const& operator()(
        msgpack::object const& o,
        dials::af::Shoebox<>& v) const {

      // Check the type is an array
      if (o.type != msgpack::type::ARRAY) {
        throw DIALS_ERROR("msgpack type is not an array");
      }

      // Check the size is 5
      if (o.via.array.size != 5) {
        throw DIALS_ERROR("msgpack array does not have correct dimensions");
      }

      // Read the shoebox structure.
      o.via.array.ptr[0].convert(v.panel);
      o.via.array.ptr[1].convert(v.bbox);
      o.via.array.ptr[2].convert(v.data);
      o.via.array.ptr[3].convert(v.mask);
      o.via.array.ptr[4].convert(v.background);
      return o;
    }
  };

  /**
   * Convert a msgpack array into a column. The array will have a structure like
   * [ type, data ]
   */
  template<>
  struct convert< dials::af::reflection_table::mapped_type > {
    msgpack::object const& operator()(
        msgpack::object const& o,
        dials::af::reflection_table::mapped_type& v) const {

      // Check the type is an array
      if (o.type != msgpack::type::ARRAY) {
        throw DIALS_ERROR("msgpack type is not an array");
      }

      // Check there are 2 elements
      if (o.via.array.size != 2) {
        throw DIALS_ERROR("msgpack array does not have correct dimensions");
      }

      // Read the type name from the first element
      std::string name;
      o.via.array.ptr[0].convert(name);

      // Read an af::shared<T> from the second element
      if (name == "bool") {
        v = extract<bool>(o.via.array.ptr[1]);
      } else if (name == "int") {
        v = extract<int>(o.via.array.ptr[1]);
      } else if (name == "std::size_t") {
        v = extract<std::size_t>(o.via.array.ptr[1]);
      } else if (name == "double") {
        v = extract<double>(o.via.array.ptr[1]);
      } else if (name == "std::string") {
        v = extract<std::string>(o.via.array.ptr[1]);
      } else if (name == "vec2<double>") {
        v = extract< scitbx::vec2<double> >(o.via.array.ptr[1]);
      } else if (name == "vec3<double>") {
        v = extract< scitbx::vec3<double> >(o.via.array.ptr[1]);
      } else if (name == "mat3<double>") {
        v = extract< scitbx::mat3<double> >(o.via.array.ptr[1]);
      } else if (name == "int6") {
        v = extract<scitbx::af::int6>(o.via.array.ptr[1]);
      } else if (name == "cctbx::miller::index<>") {
        v = extract< cctbx::miller::index<> >(o.via.array.ptr[1]);
      } else if (name == "Shoebox<>") {
        v = extract< dials::af::Shoebox<> >(o.via.array.ptr[1]);
      } else {
        throw DIALS_ERROR("Unexpected column type");
      }
      return o;
    }

    template <typename T>
    scitbx::af::shared<T> extract(msgpack::object const& o) const {
      scitbx::af::shared<T> data;
      o.convert(data);
      return data;
    }
  };

  /**
   * Convert a msgpack structure into a reflection table. The msgpack structure
   * will be like:
   * [
   *  "dials::af::reflection_table",
   *  VERSION,
   *  N_ROWS,
   *  DATA
   *  ]
   *
   * The first entry identifies the data as a reflection table.
   * The second entry gives the version number in case this changes
   * The third entry gives the expected number of rows in the table
   * The fourth entry is a map with key value pairs corresponding to the names
   * and arrays of the column data.
   */
  template <>
  struct convert<dials::af::reflection_table> {
    msgpack::object const& operator()(
        msgpack::object const& o,
        dials::af::reflection_table& v) const {
      typedef typename dials::af::reflection_table::key_type key_type;
      typedef typename dials::af::reflection_table::mapped_type mapped_type;

      // Check the type is an array
      if (o.type != msgpack::type::ARRAY) {
        throw DIALS_ERROR("msgpack type is not an array");
      }

      // Check there are 4 elements
      if (o.via.array.size != 4) {
        throw DIALS_ERROR("msgpack array does not have correct dimensions");
      }

      // Read the metadata
      std::string type;
      std::size_t version;
      std::size_t nrows;
      o.via.array.ptr[0].convert(type);
      o.via.array.ptr[1].convert(version);
      o.via.array.ptr[2].convert(nrows);

      // Check this is a reflection table
      if (type != "dials::af::reflection_table") {
        throw DIALS_ERROR("Expected dials::af::reflection_table, got something else");
      }

      // Check the version is what we expect
      if (version != 1) {
        throw DIALS_ERROR("Expected version 1, got something else");
      }

      // Resize the expected number of rows
      v.resize(nrows);

      // Get the column map
      const msgpack::object& map_object = o.via.array.ptr[3];

      // Ensure it is a column type
      if (map_object.type != msgpack::type::MAP) {
        throw DIALS_ERROR("msgpack type is not a map");
      }

      // Read the columns from the map
      if (map_object.via.map.size != 0) {
        msgpack::object_kv* first = map_object.via.map.ptr;
        msgpack::object_kv* last = first + map_object.via.map.size;
        for (msgpack::object_kv *it = first; it != last; ++it) {
          key_type key;
          mapped_type value;
          it->key.convert(key);
          it->val.convert(value);
          v[key] = value;
        }
      }
      return o;
    }
  };

} // namespace adaptor
} // MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS)
} // namespace msgpack

#endif
