#include <boost/python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <scitbx/stl/map_wrapper.h>
#include <boost/shared_ptr.hpp>
#include <map>
#include <iostream>

namespace dials { namespace util {

namespace boost_python {

using namespace boost::python;

template <class Key, class T, class Compare = std::less <Key> >
class SharedMap {

public:

  typedef typename std::map <Key, T, Compare> map_type;
  typedef typename map_type::key_type key_type;
  typedef typename map_type::mapped_type mapped_type;
  typedef typename map_type::value_type value_type;
  typedef typename map_type::key_compare key_compare;
  typedef typename map_type::value_compare value_compare;
  typedef typename map_type::allocator_type allocator_type;
  typedef typename map_type::reference reference;
  typedef typename map_type::const_reference const_reference;
  typedef typename map_type::pointer pointer;
  typedef typename map_type::const_pointer const_pointer;
  typedef typename map_type::iterator iterator;
  typedef typename map_type::const_iterator const_iterator;
  typedef typename map_type::reverse_iterator reverse_iterator;
  typedef typename map_type::const_reverse_iterator const_reverse_iterator;
  typedef typename map_type::difference_type difference_type;
  typedef typename map_type::size_type size_type;

  SharedMap() : map_(new map_type()) {}

  iterator begin() {
    return (*map_).begin();
  }
  
  const_iterator begin() const {
    return (*map_).begin();
  }

  iterator end() {
    return (*map_).end();
  }
  
  const_iterator end() const {
    return (*map_).end();
  }

  reverse_iterator rbegin() {
    return (*map_).rbegin();
  }
  
  const_reverse_iterator rbegin() const {
    return (*map_).rbegin();
  }
  
  reverse_iterator rend() {
    return (*map_).rend();
  }
  
  const_reverse_iterator rend() const {
    return (*map_).rend();
  }

  bool empty() const {
    return (*map_).empty();
  }

  size_type size() const {
    return (*map_).size();
  }
  
  size_type max_size() const {
    return (*map_).max_size();
  }
  
  mapped_type& operator[] (const key_type& k) {
    return (*map_)[k];
  }

  mapped_type& at(const key_type& k) {
    return (*map_).at(k);
  }

  const mapped_type& at(const key_type& k) const {
    return (*map_).at(k);
  }
  
  std::pair <iterator, bool> insert(const value_type& val) {
    return (*map_).insert(val);
  }

  iterator insert(iterator position, const value_type& val) {
    return (*map_).insert(position, val);
  }

  template <class InputIterator>
  void insert(InputIterator first, InputIterator last) {
    return (*map_).insert(first, last);
  }

  void erase(iterator position) {
    return (*map_).erase(position);
  }

  size_type erase(const key_type& k) {
    return (*map_).erase(k);
  }
  
  void erase(iterator first, iterator last) {
    return (*map_).erase(first, last);
  }

  void swap(map_type& x) {
    (*map_).swap(x);
  }
  
  void clear() {
    (*map_).clear();
  }
  
  key_compare key_comp() const {
    return (*map_).key_comp();
  }
  
  value_compare value_comp() const {
    return (*map_).value_comp();
  }

  iterator find(const key_type& k) {
    return (*map_).find(k);
  }

  const_iterator find(const key_type& k) const {
    return (*map_).find(k);
  }
  
  size_type count (const key_type& k) const {
    return (*map_).count(k);
  }
  
  iterator lower_bound(const key_type& k) {
    return (*map_).lower_bound(k);
  }
  
  const_iterator lower_bound(const key_type& k) const {
    return (*map_).lower_bound(k);
  }
  
  iterator upper_bound(const key_type& k) {
    return (*map_).upper_bound(k);
  }
  
  const_iterator upper_bound(const key_type& k) const {
    return (*map_).upper_bound(k);
  }
  
  std::pair <const_iterator, const_iterator> 
  equal_range(const key_type& k) const {
    return (*map_).equal_range(k);
  }
  
  std::pair <iterator, iterator> equal_range(const key_type& k) {
    return (*map_).equal_range(k);
  }
  
  allocator_type get_allocator() const {
    return (*map_).get_allocator();
  }
    
private:  

  boost::shared_ptr <map_type> map_;
};

template <class Key, class Value>
class MyMapBase : public SharedMap <Key, Value> {

public:

  typedef SharedMap <Key, Value> map_type;
  typedef typename SharedMap <Key, Value>::mapped_type mapped_type;
  typedef typename SharedMap <Key, Value>::key_type key_type;

private:

  mapped_type& at(const key_type& k) {
    return map_type::at(k);
  }

};

class MyMap : public MyMapBase <int, int> {

public:

  typedef MyMapBase <int, int> map_type;
  typedef map_type::const_iterator const_iterator;

  void print() const {
    for (const_iterator it = begin(); it != end(); ++it) {
      std::cout << it->first << " -> " << it->second << std::endl;
    }
  }
};

MyMap map_test(MyMap &map) {
 
  MyMap map2 = map;
  for (std::size_t i = 0; i < 10; ++i) {
    map2[i] = i;
  }
 
  return map2;
}

void map_print(const MyMap &map) {
  map.print();
}

void export_my_map() {
  class_<MyMap>("MyMap")
    .def(map_indexing_suite <MyMap> ());

  def("map_test", &map_test);
  def("map_print", &map_print);
}

}

}}
