
#ifndef DXTBX_MODEL_GENERIC_TREE_H
#define DXTBX_MODEL_GENERIC_TREE_H

#include <algorithm>
#include <iostream>
#include <list>

namespace dxtbx { namespace model {

  template <typename T>
  class generic_tree_node : protected std::list< generic_tree_node<T> > {
  public:

    typedef T value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;

    typedef std::list< generic_tree_node<T> > base_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::value_type node_type;
    typedef typename base_type::reference node_reference;
    typedef typename base_type::const_reference const_node_reference;
    typedef typename base_type::pointer node_pointer;
    typedef typename base_type::const_pointer const_node_pointer;
    typedef typename base_type::iterator node_iterator;
    typedef typename base_type::const_iterator const_node_iterator;

    class iterator;
    typedef const iterator const_iterator;

    generic_tree_node() {}

    generic_tree_node(const_reference value)
      : value_(value) {}

    iterator begin() {
      return iterator(base_type::begin());
    }

    const_iterator begin() const {
      return const_iterator(base_type::begin());
    }

    iterator end() {
      return iterator(base_type::end());
    }

    const_iterator end() const {
      return const_iterator(base_type::end());
    }

    bool empty() const {
      return base_type::empty();
    }

    size_type size() const {
      return base_type::size();
    }

    reference front() {
      return base_type::front().value_;
    }

    const_reference front() const {
      return base_type::front().value_;
    }

    reference back() {
      return base_type::back().value_;
    }

    const_reference back() const {
      return base_type::back().value_;
    }

    iterator push_back(const_reference v) {
      base_type::push_back(node_type(v));
      return --end();
    }

    iterator push_front(const_reference v) {
      base_type::push_front(node_type(v));
      return begin();
    }

    void pop_back() {
      base_type::pop_back();
    }

    void pop_front() {
      base_type::pop_front();
    }

    iterator insert(iterator position, const_reference v) {
      return iterator(base_type::insert(position.node_, node_type(v)));
    }

    iterator erase(iterator position) {
      return iterator(base_type::erase(position.node_));
    }

    iterator erase(iterator first, iterator last) {
      return iterator(base_type::erase(first.node_, last.node_));
    }

    void clear() {
      base_type::clear();
    }

    class iterator {
    public:

      typedef value_type value_type;
      typedef value_type& reference;
      typedef const value_type& const_reference;
      typedef value_type* pointer;
      typedef const value_type* const_pointer;

      typedef size_type size_type;
      typedef node_type node_type;
      typedef node_reference node_reference;
      typedef const_node_reference const_node_reference;
      typedef node_pointer node_pointer;
      typedef const_node_pointer const_node_pointer;

      iterator() {}

      iterator(node_iterator iter)
        : node_(iter) {}

      bool operator!() {
        return !node_;
      }

      bool operator==(const iterator &other) const {
        return node_ == other.node_;
      }

      bool operator!=(const iterator &other) const {
        return node_ != other.node_;
      }

      reference operator*() {
        return node_->value_;
      }

      pointer operator->() {
        return &node_->value_;
      }

      const_reference operator*() const {
        return node_->value_;
      }

      const_pointer operator->() const {
        return &node_->value_;
      }

      iterator& operator++() {
        ++node_;
        return *this;
      }

      iterator operator++(int) {
        iterator it = *this;
        ++(*this);
        return it;
      }

      iterator& operator--() {
        --node_;
        return *this;
      }

      iterator operator--(int) {
        iterator it = *this;
        --(*this);
        return it;
      }

      iterator begin() {
        return node_->begin();
      }

      const_iterator begin() const {
        return node_->begin();
      }

      iterator end() {
        return node_->end();
      }

      const_iterator end() const {
        return node_->end();
      }

      bool empty() const {
        return node_->empty();
      }

      size_type size() const {
        return node_->size();
      }

      reference front() {
        return node_->front();
      }

      const_reference front() const {
        return node_->front();
      }

      reference back() {
        return node_->back();
      }

      const_reference back() const {
        return node_->back();
      }

      iterator push_back(const_reference v) {
        return node_->push_back(v);
      }

      iterator push_front(const_reference v) {
        return node_->push_front(v);
      }

      void pop_back() {
        node_->pop_back();
      }

      void pop_front() {
        node_->pop_front();
      }

      iterator insert(iterator position, const_reference v) {
        return node_->insert(position, v);
      }

      iterator erase(iterator position) {
        return node_->erase(position);
      }

      iterator erase(iterator first, iterator last) {
        return node_->erase(first, last);
      }

      void clear() {
        node_->clear();
      }

      node_iterator node() {
        return node_;
      }

      const_node_iterator node() const {
        return node_;
      }

      iterator next_sibling() {
        iterator it = *this;
        return ++it;
      }

      const_iterator next_sibling() const {
        iterator it = *this;
        return ++it;
      }

      iterator prev_sibling() {
        iterator it = *this;
        return --it;
      }

      const_iterator prev_sibling() const {
        iterator it = *this;
        return --it;
      }

      iterator first_child() {
        return begin();
      }

      const_iterator first_child() const {
        return begin();
      }

      iterator last_child() {
        return --end();
      }

//      const_iterator last_child() const {
//        return --end();
//      }

      bool leaf() const {
        return empty();
      }

    protected:

      node_iterator node_;
    };

  private:

    value_type value_;
  };

  template <typename T>
  class generic_tree : public generic_tree_node<T>::iterator {
  public:

    typedef typename generic_tree_node<T>::iterator base_type;
    typedef typename generic_tree_node<T>::iterator iterator;
    typedef typename generic_tree_node<T>::const_iterator const_iterator;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::const_reference const_reference;
    typedef typename base_type::pointer pointer;
    typedef typename base_type::const_pointer const_pointer;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::node_type node_type;

    generic_tree() {
      base_type::node_ = root_node_.push_back(value_type()).node();
    }

    generic_tree(const value_type &value) {
      base_type::node_ = root_node_.push_back(value).node();
    }

  private:
    node_type root_node_;
  };

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_GENERIC_TREE_H
