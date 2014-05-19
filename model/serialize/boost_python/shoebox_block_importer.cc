/*
 * shoebox_block_importer.cc
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
#include <boost_adaptbx/std_pair_conversion.h>
#include <dials/model/serialize/shoebox_block_importer.h>

namespace dials { namespace model { namespace serialize {
  namespace boost_python {

  using namespace boost::python;

  ShoeboxBlockImporter* make_shoebox_block_importer(
      const std::string &filename,
      const af::const_ref<std::size_t> &blocks,
      const af::const_ref<double> &z,
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
    return new ShoeboxBlockImporter(filename, blocks, z,
      gm.const_ref(), dm.const_ref(), mm.const_ref());
  }

  /**
   * A proxy iterator to iterate over the shoeboxes
   */
  class ShoeboxBlockIterator {
  public:

    typedef std::forward_iterator_tag iterator_category;
    typedef std::size_t base_iterator;
    typedef ptrdiff_t difference_type;
    typedef af::shared<std::size_t> index_array_type;
    typedef af::shared< Shoebox<> > shoebox_array_type;
    typedef std::pair<index_array_type, shoebox_array_type> value_type;
    typedef const value_type *pointer;
    typedef value_type reference;

    ShoeboxBlockIterator(ShoeboxBlockImporter &importer, base_iterator it)
      : importer_(importer),
        it_(it) {}

    reference operator*() {
      return importer_[it_];
    }

    ShoeboxBlockIterator& operator++() {
      ++it_;
      return *this;
    }

    ShoeboxBlockIterator operator++(int) {
      ShoeboxBlockIterator result(*this);
      ++(*this);
      return result;
    }

    bool operator==(const ShoeboxBlockIterator& rhs) const {
      return it_ == rhs.it_;
    }

    bool operator!=(const ShoeboxBlockIterator& rhs) const {
      return !(*this == rhs);
    }

  private:
    ShoeboxBlockImporter &importer_;
    base_iterator it_;
  };

  /**
   * Struct to help in creation of table proxy iterators
   */
  struct make_shoebox_block_iterator {
    static
    ShoeboxBlockIterator begin(ShoeboxBlockImporter &self) {
      return ShoeboxBlockIterator(self, 0);
    }

    static
    ShoeboxBlockIterator end(ShoeboxBlockImporter &self) {
      return ShoeboxBlockIterator(self, self.size());
    }

    static
    object range() {
      return boost::python::range(
        &make_shoebox_block_iterator::begin,
        &make_shoebox_block_iterator::end);
    }
  };

  void export_shoebox_block_importer()
  {
    // Export the importer class
    class_<ShoeboxBlockImporter, boost::noncopyable>(
        "ShoeboxBlockImporter", no_init)
      .def(init<const std::string&,
                const af::const_ref<std::size_t>&,
                const af::const_ref<double>&>((
          arg("filename"),
          arg("blocks"),
          arg("z"))))
      .def("__init__", make_constructor(
        &make_shoebox_block_importer,
        default_call_policies(), (
          arg("filename"),
          arg("blocks"),
          arg("z"),
          arg("gain"),
          arg("dark"),
          arg("mask"))))
      .def("__len__", &ShoeboxBlockImporter::size)
      .def("__getitem__", &ShoeboxBlockImporter::operator[])
      .def("__iter__", make_shoebox_block_iterator::range());
      ;

    // Export a conversion for the pair
    boost_adaptbx::std_pair_conversions::to_tuple<
      af::shared<std::size_t>, af::shared< Shoebox<> > >();
  }

}}}} // namespace dials::model::serialize::boost_python
