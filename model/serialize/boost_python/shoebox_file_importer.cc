/*
 * shoebox_file_importer.cc
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
#include <dials/model/serialize/shoebox_file_importer.h>

namespace dials { namespace model { namespace serialize {
  namespace boost_python {

  using namespace boost::python;

  ShoeboxFileImporter* make_shoebox_file_importer(
      const std::string &filename,
      boost::python::tuple gmt,
      boost::python::tuple dmt,
      boost::python::tuple mmt) {

    // The input
    af::shared<ShoeboxFileImporter::gain_map_ref_type> gm;
    af::shared<ShoeboxFileImporter::dark_map_ref_type> dm;
    af::shared<ShoeboxFileImporter::mask_map_ref_type> mm;

    // Ensure all the same length
    std::size_t num = len(gmt);
    DIALS_ASSERT(len(dmt) == num);
    DIALS_ASSERT(len(mmt) == num);

    // Extract the stuff
    for (std::size_t i = 0; i < num; ++i) {
      gm.push_back(extract<ShoeboxFileImporter::gain_map_ref_type>(gmt[i]));
      dm.push_back(extract<ShoeboxFileImporter::dark_map_ref_type>(dmt[i]));
      mm.push_back(extract<ShoeboxFileImporter::mask_map_ref_type>(mmt[i]));
    }

    // Return the new importer
    return new ShoeboxFileImporter(filename,
      gm.const_ref(), dm.const_ref(), mm.const_ref());
  }

  /**
   * A proxy iterator to iterate over the shoeboxes
   */
  class ShoeboxIterator {
  public:

    typedef std::forward_iterator_tag iterator_category;
    typedef std::size_t base_iterator;
    typedef ptrdiff_t difference_type;
    typedef Shoebox<> value_type;
    typedef const value_type *pointer;
    typedef value_type reference;

    ShoeboxIterator(ShoeboxFileImporter &importer, base_iterator it)
      : importer_(importer),
        it_(it) {}

    reference operator*() {
      return importer_[it_];
    }

    ShoeboxIterator& operator++() {
      ++it_;
      return *this;
    }

    ShoeboxIterator operator++(int) {
      ShoeboxIterator result(*this);
      ++(*this);
      return result;
    }

    bool operator==(const ShoeboxIterator& rhs) const {
      return it_ == rhs.it_;
    }

    bool operator!=(const ShoeboxIterator& rhs) const {
      return !(*this == rhs);
    }

  private:
    ShoeboxFileImporter &importer_;
    base_iterator it_;
  };

  /**
   * Struct to help in creation of table proxy iterators
   */
  struct make_shoebox_iterator {
    static
    ShoeboxIterator begin(ShoeboxFileImporter &self) {
      return ShoeboxIterator(self, 0);
    }

    static
    ShoeboxIterator end(ShoeboxFileImporter &self) {
      return ShoeboxIterator(self, self.size());
    }

    static
    object range() {
      return boost::python::range(
        &make_shoebox_iterator::begin,
        &make_shoebox_iterator::end);
    }
  };

  void export_shoebox_file_importer()
  {
    // Typedefs of select function
    typedef af::shared< Shoebox<> > (ShoeboxFileImporter::*select_range)(
        std::size_t, std::size_t);
    typedef af::shared< Shoebox<> > (ShoeboxFileImporter::*select_many)(
        const af::const_ref<std::size_t>&);

    // Export the importer class
    class_<ShoeboxFileImporter, boost::noncopyable>(
        "ShoeboxFileImporter", no_init)
      .def(init<const std::string&>((
          arg("filename"))))
      .def("__init__", make_constructor(
        &make_shoebox_file_importer,
        default_call_policies(), (
          arg("filename"),
          arg("gain"),
          arg("dark"),
          arg("mask"))))
      .def("blob", &ShoeboxFileImporter::blob)
      .def("__len__", &ShoeboxFileImporter::size)
      .def("bboxes", &ShoeboxFileImporter::bboxes)
      .def("panels", &ShoeboxFileImporter::panels)
      .def("__getitem__", &ShoeboxFileImporter::operator[])
      .def("select", (select_range)&ShoeboxFileImporter::select)
      .def("select", (select_many)&ShoeboxFileImporter::select)
      .def("__iter__", make_shoebox_iterator::range());
      ;
  }

}}}} // namespace dials::model::serialize::boost_python
