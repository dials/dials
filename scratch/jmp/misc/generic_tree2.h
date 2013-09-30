
#ifndef DXTBX_MODEL_GENERIC_TREE_H
#define DXTBX_MODEL_GENERIC_TREE_H

#include <algorithm>
#include <iostream>
#include <list>
#include <boost/shared_ptr.hpp>

namespace dxtbx { namespace model {

  template <typename T>
  class generic_tree : public std::list< generic_tree<T> > , public T {
  public:

    typedef std::list< generic_tree<T> >  base_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::const_reference const_reference;
    typedef typename base_type::pointer pointer;
    typedef typename base_type::const_pointer const_pointer;
    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;
    typedef typename base_type::reverse_iterator reverse_iterator;
    typedef typename base_type::const_reverse_iterator const_reverse_iterator;

    generic_tree() {}

    generic_tree(int a) : T(a) {}

    generic_tree(const generic_tree &a)
      : base_type(a), T(a) {
      std::cout << "copy" << a.a << std::endl;
    }

  };

//  template <typename T>
//  class generic_tree : public generic_tree_node<T>::iterator {
//  public:

//    typedef typename generic_tree_node<T>::iterator base_type;
//    typedef typename generic_tree_node<T>::iterator iterator;
//    typedef typename generic_tree_node<T>::const_iterator const_iterator;
//    typedef typename base_type::value_type value_type;
//    typedef typename base_type::reference reference;
//    typedef typename base_type::const_reference const_reference;
//    typedef typename base_type::pointer pointer;
//    typedef typename base_type::const_pointer const_pointer;
//    typedef typename base_type::size_type size_type;
//    typedef typename base_type::node_type node_type;

//    generic_tree() {
//      base_type::node_ = root_node_.push_back(value_type()).node();
//    }
//
//    generic_tree(const value_type &value) {
//      base_type::node_ = root_node_.push_back(value).node();
//    }

//  private:
//    node_type root_node_;
//  };

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_GENERIC_TREE_H
