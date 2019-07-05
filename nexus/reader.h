
#ifndef DIALS_NEXUS_READER_H
#define DIALS_NEXUS_READER_H

#include <H5Cpp.h>
#include <iostream>
#include <vector>
#include <dials/error.h>

namespace dials { namespace nexus {

  namespace detail {

    /* template <typename Object, typename Visitor> */
    /* void visitgroups(Object &obj, Visitor visitor) { */
    /*   for (std::size_t i = 0; i < obj.getNumObjs(); ++i) { */
    /*     std::string child_name = obj.getObjnameByIdx(i); */
    /*     H5G_obj_t   child_type = obj.getObjTypeByIdx(i); */
    /*     if (child_type == H5G_GROUP) { */
    /*       visitor(std::size_t i, obj.openGroup(child_name)); */
    /*     } */
    /*   } */
    /* } */

    /* template <typename Iterator> */
    /* struct find_nxmx_entries_visitor { */

    /*   Iterator out; */

    /*   find_nxmx_entries_visitor(Iterator out_) */
    /*     : out(out_) {} */

    /*   void operator()(std::size_t index, Group &group) const { */
    /*     if (group.attrExists("NX_class")) { */
    /*       H5::Attribute nx_class = group.openAttribute("NX_class"); */
    /*       H5::DataType dtype = nx_class.getDataType(); */
    /*       if (dtype.isVariableStr()) { */
    /*         std::string name; */
    /*         nx_class.read(dtype, name); */
    /*         if (name == "NXentry" || name == "NXsubentry") { */
    /*           *out++ = index; */
    /*           visitgroups(groups, find_nxmx_entries_visitor<Iterator>(out)); */
    /*         } */
    /*       } */
    /*     } */
    /*   } */

    /* }; */

    /* template <typename Object, typename OutputIterator> */
    /* void find_nxmx_entries(Object &obj, OuputIterator out) { */
    /*   visitgroups(obj, find_nxmx_entries_visitor<OutputIterator>(out)); */
    /* } */

    /* template <typename InputIterator, typename OutputIterator> */
    /* void find_nx_classes( */
    /*     const H5::CommonFG &obj, */
    /*     InputIterator first, */
    /*     InputIterator last, */
    /*     OutputIterator out) { */
    /*   for (std::size_t i = 0; i < obj.getNumObjs(); ++i) { */
    /*     std::string child_name = obj.getObjnameByIdx(i); */
    /*     H5G_obj_t   child_type = obj.getObjTypeByIdx(i); */
    /*     if (child_type == H5G_GROUP) { */
    /*       H5::Group group = obj.openGroup(child_name); */
    /*       if (group.attrExists("NX_class")) { */
    /*         H5::Attribute nx_class = group.openAttribute("NX_class"); */
    /*         H5::DataType dtype = nx_class.getDataType(); */
    /*         if (dtype.isVariableStr()) { */
    /*           std::string name; */
    /*           nx_class.read(dtype, name); */
    /*           for (InputIterator it = first; it != last; ++it) { */
    /*             if (name == *it) { */
    /*               *out++ = i; */
    /*             } */
    /*           } */
    /*         } */
    /*       } */
    /*     } */
    /*   } */
    /* } */

    /* template <typename OutputIterator> */
    /* void find_nxmx_entries(const H5::CommonFG &obj, OutputIterator out) { */
    /*   std::vector<std::string> names(2); */
    /*   names[0] = "NXentry"; */
    /*   names[1] = "NXsubentry"; */
    /*   std::vector<std::size_t> temp; */
    /*   find_nx_classes(obj, names.begin(), names.end(), std::back_inserter(temp)); */
    /* } */

    /* def find_nx_mx_entries(nx_file, entry): */
    /* ''' */
    /* Find NXmx entries */

    /* ''' */
    /* hits = [] */
    /* def visitor(name, obj): */
    /* if "NX_class" in obj.attrs.keys(): */
    /* if obj.attrs["NX_class"] in ["NXentry", "NXsubentry"]: */
    /*   if "definition" in obj.keys(): */
    /*     if obj["definition"].value == "NXmx": */
    /*       hits.append(obj) */
    /* nx_file[entry].visititems(visitor) */
    /* return hits */

  }  // namespace detail

  /* class NXinstrument { */
  /* public: */

  /* }; */

  /*     "depends_on" : { */
  /*       "minOccurs" : 1, */
  /*       "tests" : [ */
  /*         check_depends_on() */
  /*       ] */
  /*     }, */
  /*     "unit_cell" : { */
  /*       "minOccurs" : 1, */
  /*       "tests" : [ */
  /*         check_dset(dtype="float64", dims=2) */
  /*       ] */
  /*     }, */
  /*     "sample_orientation" : { */
  /*       "minOccurs" : 0, */
  /*       "tests" : [ */
  /*         check_dset(dtype="float64", shape=(3,)) */
  /*       ] */
  /*     }, */
  /*     "orientation_matrix" : { */
  /*       "minOccurs" : 1, */
  /*       "tests" : [ */
  /*         check_dset(dtype="float64", dims=3) */
  /*       ] */
  /*     }, */
  /* class NXsample { */
  /* public: */

  /*   std::string name() const { */
  /*     return ""; */
  /*   } */

  /*   std::string unit_cell_group() const { */
  /*     return ""; */
  /*   } */

  /*   NXbeam beam() const { */
  /*     return NXbeam(); */
  /*   } */

  /*   NXtransformations transformations() const { */
  /*     return NXtransformations(); */
  /*   } */

  /* }; */

  /* class NXdata { */
  /* public: */

  /* }; */

  /* class NXmx { */
  /* public: */

  /*   std::string title() const { */
  /*     return ""; */
  /*   } */

  /*   std::string start_time() const { */
  /*     return ""; */
  /*   } */

  /*   std::string end_time() const { */
  /*     return ""; */
  /*   } */

  /*   NXinstrument instrument() const { */
  /*     return NXinstrument(); */
  /*   } */

  /*   NXsample sample() const { */
  /*     return NXsample(); */
  /*   } */

  /*   NXdata data() const { */
  /*     return NXdata(); */
  /*   } */

  /* }; */

  class Reader {
    Object root_;

  public:
    Reader(const char *filename) : root_(H5::File(filename, H5F_ACC_RDONLY)) {}

    Object root() const {
      return root_;
    }
  };

}}  // namespace dials::nexus

#endif  // DIALS_NEXUS_READER_H
