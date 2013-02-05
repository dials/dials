#include <boost/python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include "multimap_indexing_suite.h"
#include <scitbx/stl/map_wrapper.h>
#include <boost/python/iterator.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/copy_non_const_reference.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/iterator.hpp>
#include <map>
#include <iostream>
#include <scitbx/boost_python/container_conversions.h>

namespace dials { namespace util {

namespace boost_python {

using namespace boost::python;

struct MyData {
  int value1;
  int value2;
  MyData() : value1(0), value2(0) {}
  MyData(int v1, int v2) : value1(v1), value2(v2) {}
};

struct MyData2 {
  int value1;
  int value2;
  MyData2() : value1(0), value2(0) {}
  MyData2(int v1, int v2) : value1(v1), value2(v2) {}
};

typedef std::multimap <int,MyData> MyMap;

void insert(MyMap &map, const MyMap::key_type &key, const MyMap::mapped_type &value) {
    map.insert(MyMap::value_type(key, value));
}

template<typename PairType>
struct PairToTupleConverter {
  static PyObject* convert(const PairType& pair) {
    return incref(make_tuple(pair.first, pair.second).ptr());
  }
};



//typedef std::pair <MyMap::iterator, MyMap::iterator> iterator_range_base;

//template <typename TIter>
//struct my_iterator_range {
//    typedef TIter iterator;
//    typedef std::pair <TIter, TIter> iterator_pair;
//    
//    my_iterator_range();
//    my_iterator_range(const iterator_pair &range_) {
//        this->range.first = range_.first;
//        this->range.second = range_.second;
//    }
//    MyMap::size_type size() const {
//        return std::distance(this->range.first, this->range.second);
//    }
//    iterator_pair range;
//    iterator& begin() { return this->range.first; }
//    iterator& end() { return this->range.second; }
//};


//typedef my_iterator_range <MyMap::iterator> MyIteratorRange;
typedef std::pair <MyMap::iterator, MyMap::iterator> MyIteratorRange;

//template <typename T1, typename T2, typename T3>
//class Triple {
//    typedef T1 first_type;
//    typedef T2 second_type;
//    typedef T3 third_type;
//    
//    first_type first;
//    second_type second;
//    third_type third;
//};

unsigned int iterator_range_size(MyIteratorRange &range) {
    return std::distance(range.first, range.second);
}

MyIteratorRange equal_range(MyMap &map, const MyMap::key_type &k) {
    MyIteratorRange range = map.equal_range(k);
    if (range.first == range.second) {
        PyErr_SetString(PyExc_KeyError, "Invalid key");
        throw_error_already_set();
    }
    return range;
}

void export_my_map() {




    MyMap::iterator (MyMap::*find)(const MyMap::key_type&) = &MyMap::find;
    MyMap::size_type (MyMap::*erase_by_key)(const MyMap::key_type&) = &MyMap::erase;
    
    to_python_converter<MyMap::value_type, PairToTupleConverter<MyMap::value_type> >();

    class_<MyMap::value_type>("MyMapValueType")
        .def(init <MyMap::value_type::first_type,
                   MyMap::value_type::second_type>())
        .def_readonly("first", &MyMap::value_type::first)
        .def_readwrite("second", &MyMap::value_type::second);
        
    class_<MyData> ("MyData")
        .def(init <int, int>())
        .def_readwrite("value1", &MyData::value1)
        .def_readwrite("value2", &MyData::value2);
        
    class_<MyData2> ("MyData2")
        .def(init <int, int>())
        .def_readwrite("value1", &MyData2::value1)
        .def_readwrite("value2", &MyData2::value2);
        
    class_<MyIteratorRange>("MyIteratorRange")
        .def("__iter__",
            range <return_internal_reference<> >(
                &MyIteratorRange::first, &MyIteratorRange::second))
        .def("__len__", &iterator_range_size);
   
         
    class_ <MyMap> ("MyMap")
        .def("__iter__", iterator<MyMap, return_internal_reference<> >())
        .def("__len__", &MyMap::size)
        .def("clear", &MyMap::clear)
        .def("find", find)
        .def("insert", &insert)
        .def("erase", erase_by_key)
        .def("count", &MyMap::count)
        .def("__getitem__", equal_range);
        //.def("__setitem__", set_item);
        
    typedef std::map <int, MyData2> TestMap;
        
    class_ <TestMap> ("TestMap")
        .def(map_indexing_suite <TestMap, true >());
}

}

}}
