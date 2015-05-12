
#ifndef DIALS_NEXUS_SERIALIZE_H
#define DIALS_NEXUS_SERIALIZE_H

#include <string>
#include <dials/error.h>

namespace dials { namespace nexus {

  template <typename T>
  struct serialize {

    template <typename Handle>
    static
    T load(const Handle &handle);

    template <typename Handle>
    static
    void dump(const T &obj, Handle &handle);

  };

  template <>
  struct serialize<std::string> {

    template <typename Handle>
    static
    std::string load(const Handle &dataset) {
      std::string result;
      H5::DataType datatype = dataset.getDataType();
      H5::DataSpace dataspace = dataset.getSpace();
      DIALS_ASSERT(datatype.isVariableStr());
      DIALS_ASSERT(dataspace.isSimple());
      int ndims = dataspace.getSimpleExtentNdims();
      DIALS_ASSERT(ndims == 1);
      hsize_t dims = 0;
      dataspace.getSimpleExtentDims(&dims);
      DIALS_ASSERT(dims == 1);
      dataset.read(result, datatype);
      return result;
    }

    template <typename Handle>
    static
    void dump(const std::string &obj, Handle &handle) {

    }

  };

  template <typename Handle>
  bool is_nx_class(const Handle &handle, std::string test_name) {
    if (handle.attrExists("NX_class")) {
      H5::Attribute attr = handle.openAttribute("NX_class");
      H5::DataType dtype = attr.getDataType();
      if (dtype.isVariableStr()) {
        std::string name;
        attr.read(dtype, name);
        if (name == test_name) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename Handle>
  bool is_nxmx_entry(const Handle &handle) {
    try {
      if ("NXmx" == serialize<std::string>::load(
            handle.openDataSet("definition"))) {
        return true;
      }
    } catch(H5::Exception) {
      // Do nothing
    } catch(std::exception) {
      // Do nothing
    }
    return false;
  }

  template <typename Handle>
  bool is_nx_entry(const Handle &handle) {
    if (handle.attrExists("NX_class")) {
      H5::Attribute attr = handle.openAttribute("NX_class");
      H5::DataType dtype = attr.getDataType();
      if (dtype.isVariableStr()) {
        std::string name;
        attr.read(dtype, name);
        if (name == "NXentry" || name == "NXsubentry") {
          return true;
        }
      }
    }
    return false;
  }

  template <typename Handle, typename Iterator>
  void find_nx_entries(Handle &handle, Iterator out) {
    for (std::size_t i = 0; i < handle.getNumObjs(); ++i) {
      if (handle.getObjTypeByIdx(i) == H5G_GROUP) {
        std::string name = handle.getObjnameByIdx(i);
        H5::Group group = handle.openGroup(name);
        if (is_nx_entry(group)) {
          *out++ = group;
          find_nx_entries(group, out);
        }
      }
    }
  }

}} // namespace dials::nexus

#endif // DIALS_NEXUS_SERIALIZE_H
