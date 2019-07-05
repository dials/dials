
#ifndef DIALS_NEXUS_SERIALIZE_H
#define DIALS_NEXUS_SERIALIZE_H

#include <string>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace nexus {

  using scitbx::vec2;
  using scitbx::vec3;

  std::string dataset_name(const H5::DataSet &ds) {
    size_t len = H5Iget_name(ds.getId(), NULL, 0);
    char buffer[len];
    H5Iget_name(ds.getId(), buffer, len + 1);
    std::string n = buffer;
    return n;
  }

  template <typename T>
  struct serialize {
    template <typename Handle>
    static T load(const Handle &handle);

    template <typename Handle>
    static void dump(const T &obj, Handle &handle);
  };

  template <>
  struct serialize<std::string> {
    template <typename Handle>
    static std::string load(const Handle &dataset) {
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
    static void dump(const std::string &obj, Handle &handle) {}
  };

  template <>
  struct serialize<bool> {
    template <typename Handle>
    static bool load(const Handle &dataset) {
      bool result;
      H5::DataType datatype = dataset.getDataType();
      H5::DataSpace dataspace = dataset.getSpace();
      DIALS_ASSERT(dataspace.isSimple());
      bool ndims = dataspace.getSimpleExtentNdims();
      DIALS_ASSERT(ndims == 1);
      hsize_t dims = 0;
      dataspace.getSimpleExtentDims(&dims);
      DIALS_ASSERT(dims == 1);
      dataset.read(&result, H5::PredType::NATIVE_HBOOL);
      return result;
    }

    template <typename Handle>
    static void dump(const bool &obj, Handle &handle) {}
  };

  template <>
  struct serialize<int> {
    template <typename Handle>
    static int load(const Handle &dataset) {
      int result;
      H5::DataType datatype = dataset.getDataType();
      H5::DataSpace dataspace = dataset.getSpace();
      DIALS_ASSERT(dataspace.isSimple());
      int ndims = dataspace.getSimpleExtentNdims();
      DIALS_ASSERT(ndims == 1);
      hsize_t dims = 0;
      dataspace.getSimpleExtentDims(&dims);
      DIALS_ASSERT(dims == 1);
      dataset.read(&result, H5::PredType::NATIVE_INT);
      return result;
    }

    template <typename Handle>
    static void dump(const int &obj, Handle &handle) {}
  };

  template <>
  struct serialize<double> {
    template <typename Handle>
    static double load(const Handle &dataset) {
      double result;
      H5::DataType datatype = dataset.getDataType();
      H5::DataSpace dataspace = dataset.getSpace();
      DIALS_ASSERT(dataspace.isSimple());
      int ndims = dataspace.getSimpleExtentNdims();
      DIALS_ASSERT(ndims == 1);
      hsize_t dims = 0;
      dataspace.getSimpleExtentDims(&dims);
      DIALS_ASSERT(dims == 1);
      dataset.read(&result, H5::PredType::NATIVE_DOUBLE);
      return result;
    }

    template <typename Handle>
    static void dump(const double &obj, Handle &handle) {}
  };

  template <>
  struct serialize<vec2<int> > {
    template <typename Handle>
    static vec2<int> load(const Handle &dataset) {
      vec2<int> result;
      H5::DataType datatype = dataset.getDataType();
      H5::DataSpace dataspace = dataset.getSpace();
      DIALS_ASSERT(dataspace.isSimple());
      vec2<int> ndims = dataspace.getSimpleExtentNdims();
      DIALS_ASSERT(ndims == 1);
      hsize_t dims = 0;
      dataspace.getSimpleExtentDims(&dims);
      DIALS_ASSERT(dims == 2);
      dataset.read(&result[0], H5::PredType::NATIVE_INT);
      return result;
    }

    template <typename Handle>
    static void dump(const vec2<int> &obj, Handle &handle) {}
  };

  template <>
  struct serialize<af::shared<double> > {
    template <typename Handle>
    static af::shared<double> load(const Handle &dataset) {
      af::shared<double> result;
      H5::DataType datatype = dataset.getDataType();
      H5::DataSpace dataspace = dataset.getSpace();
      DIALS_ASSERT(dataspace.isSimple());
      int ndims = dataspace.getSimpleExtentNdims();
      DIALS_ASSERT(ndims == 1);
      hsize_t dims = 0;
      dataspace.getSimpleExtentDims(&dims);
      result.resize(dims);
      dataset.read(&result[0], H5::PredType::NATIVE_DOUBLE);
      return result;
    }

    template <typename Handle>
    static void dump(const af::shared<double> &obj, Handle &handle) {}
  };

  template <>
  struct serialize<af::versa<double, af::c_grid<2> > > {
    template <typename Handle>
    static af::versa<double, af::c_grid<2> > load(const Handle &dataset) {
      af::versa<double, af::c_grid<2> > result;
      H5::DataType datatype = dataset.getDataType();
      H5::DataSpace dataspace = dataset.getSpace();
      DIALS_ASSERT(dataspace.isSimple());
      int ndims = dataspace.getSimpleExtentNdims();
      DIALS_ASSERT(ndims == 2);
      hsize_t dims[2] = {0, 0};
      dataspace.getSimpleExtentDims(dims);
      result.resize(af::c_grid<2>(dims[0], dims[1]));
      dataset.read(&result[0], H5::PredType::NATIVE_DOUBLE);
      return result;
    }

    template <typename Handle>
    static void dump(const af::versa<double, af::c_grid<2> > &obj, Handle &handle) {}
  };

  template <>
  struct serialize<af::versa<int, af::c_grid<2> > > {
    template <typename Handle>
    static af::versa<int, af::c_grid<2> > load(const Handle &dataset) {
      af::versa<int, af::c_grid<2> > result;
      H5::DataType datatype = dataset.getDataType();
      H5::DataSpace dataspace = dataset.getSpace();
      DIALS_ASSERT(dataspace.isSimple());
      int ndims = dataspace.getSimpleExtentNdims();
      DIALS_ASSERT(ndims == 2);
      hsize_t dims[2] = {0, 0};
      dataspace.getSimpleExtentDims(dims);
      result.resize(af::c_grid<2>(dims[0], dims[1]));
      dataset.read(&result[0], H5::PredType::NATIVE_INT);
      return result;
    }

    template <typename Handle>
    static void dump(const af::versa<int, af::c_grid<2> > &obj, Handle &handle) {}
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
      if ("NXmx" == serialize<std::string>::load(handle.openDataSet("definition"))) {
        return true;
      }
    } catch (H5::Exception) {
      // Do nothing
    } catch (std::exception) {
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

}}  // namespace dials::nexus

#endif  // DIALS_NEXUS_SERIALIZE_H
