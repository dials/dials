/*
 * column_table.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_FRAMEWORK_TABLE_COLUMN_TABLE_H
#define DIALS_FRAMEWORK_TABLE_COLUMN_TABLE_H

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/signals2/signal.hpp>
#include <boost/signals2/shared_connection_block.hpp>
#include <boost/variant.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/remove_if.hpp>
#include <boost/mpl/transform.hpp>
#include <algorithm>
#include <vector>
#include <map>
#include <dials/error.h>

namespace dials { namespace framework {

  namespace bs2 = boost::signals2;

  class column_synchronizer_core {
  public:

    typedef std::size_t size_type;
    typedef bs2::connection connection;
    typedef bs2::shared_connection_block connection_block;
    typedef bs2::signal<void(size_type)> resize_signal_type;
    typedef bs2::signal<void(size_type,size_type)> insert_signal_type;
    typedef bs2::signal<void(size_type,size_type)> erase_signal_type;
    typedef resize_signal_type::slot_type resize_slot_type;
    typedef insert_signal_type::slot_type insert_slot_type;
    typedef erase_signal_type::slot_type erase_slot_type;

    column_synchronizer_core()
      : size_(0) {}

    column_synchronizer_core(size_type n)
      : size_(n) {}

    void resize(size_type n) {
      size_ = n;
      resize_(n);
    }

    void insert(size_type pos, size_type n) {
      size_ += n;
      insert_(pos, n);
    }

    void erase(size_type first, size_type last) {
      DIALS_ASSERT(last > first);
      DIALS_ASSERT(last <= size_);
      size_ -= (last - first);
      erase_(first, last);
    }

    size_type size() const {
      return size_;
    }

    connection connect_resize(resize_slot_type slot) {
      return resize_.connect(slot);
    }

    connection connect_insert(insert_slot_type slot) {
      return insert_.connect(slot);
    }

    connection connect_erase(erase_slot_type slot) {
      return erase_.connect(slot);
    }

  protected:
    size_type size_;
    resize_signal_type resize_;
    insert_signal_type insert_;
    erase_signal_type erase_;
  };


  class column_synchronizer {
  public:

    typedef column_synchronizer_core::size_type size_type;
    typedef column_synchronizer_core::resize_signal_type resize_signal_type;
    typedef column_synchronizer_core::insert_signal_type insert_signal_type;
    typedef column_synchronizer_core::erase_signal_type erase_signal_type;
    typedef column_synchronizer_core::resize_slot_type resize_slot_type;
    typedef column_synchronizer_core::insert_slot_type insert_slot_type;
    typedef column_synchronizer_core::erase_slot_type erase_slot_type;
    typedef column_synchronizer_core::connection connection;
    typedef column_synchronizer_core::connection_block connection_block;

    column_synchronizer()
      : core_(boost::make_shared<column_synchronizer_core>()) {}

    column_synchronizer(size_type n)
      : core_(boost::make_shared<column_synchronizer_core>(n)) {}

    void resize(size_type n) {
      core_->resize(n);
    }

    void insert(size_type pos, size_type n) {
      core_->insert(pos, n);
    }

    void erase(size_type first, size_type last) {
      core_->erase(first, last);
    }

    size_type size() const {
      return core_->size();
    }

    connection connect_resize(resize_slot_type slot) {
      return core_->connect_resize(slot);
    }

    connection connect_insert(insert_slot_type slot) {
      return core_->connect_insert(slot);
    }

    connection connect_erase(erase_slot_type slot) {
      return core_->connect_erase(slot);
    }

  protected:

    boost::shared_ptr<column_synchronizer_core> core_;
  };



  template <typename T>
  class column_data_core {
    struct resizer;
    struct inserter;
    struct eraser;

  public:

    typedef std::vector<T> storage_type;
    typedef typename storage_type::value_type value_type;
    typedef typename storage_type::reference reference;
    typedef typename storage_type::iterator iterator;
    typedef typename storage_type::const_iterator const_iterator;
    typedef typename storage_type::size_type size_type;
    typedef typename storage_type::difference_type difference_type;

    column_data_core()
      : sync_(0),
        resize_conn_(sync_.connect_resize(resizer(&storage_))),
        insert_conn_(sync_.connect_insert(inserter(&storage_))),
        erase_conn_(sync_.connect_erase(eraser(&storage_))) {}

    template <typename InputIterator>
    column_data_core(InputIterator first, InputIterator last)
      : storage_(first, last),
        sync_(storage_.size()),
        resize_conn_(sync_.connect_resize(resizer(&storage_))),
        insert_conn_(sync_.connect_insert(inserter(&storage_))),
        erase_conn_(sync_.connect_erase(eraser(&storage_))) {}

    column_data_core(column_synchronizer sync)
      : storage_(sync.size()),
        sync_(sync),
        resize_conn_(sync_.connect_resize(resizer(&storage_))),
        insert_conn_(sync_.connect_insert(inserter(&storage_))),
        erase_conn_(sync_.connect_erase(eraser(&storage_))) {}

    reference operator[](size_type index) {
      return storage_[index];
    }

    const_iterator begin() const {
      return storage_.begin();
    }

    const_iterator end() const {
      return storage_.end();
    }

    iterator begin() {
      return storage_.begin();
    }

    iterator end() {
      return storage_.end();
    }

    size_type size() const {
      return sync_.size();
    }

    void resize(size_type n) {
      sync_.resize(n);
    }

    void push_back(const value_type &v) {
      insert(end(), v);
    }

    iterator insert(iterator position, const value_type &v) {
      size_type index = std::distance(begin(), position);
      column_synchronizer::connection_block block(insert_conn_);
      iterator result = storage_.insert(position, v);
      sync_.insert(index, 1);
      DIALS_ASSERT(storage_.size() == sync_.size());
      return result;
    }

    void insert(iterator position, size_type n, const value_type &v) {
      size_type index = std::distance(begin(), position);
      column_synchronizer::connection_block block(insert_conn_);
      storage_.insert(position, n, v);
      sync_.insert(index, n);
      DIALS_ASSERT(storage_.size() == sync_.size());
    }

    template <typename InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last) {
      size_type n = std::distance(first, last);
      size_type index = std::distance(begin(), position);
      column_synchronizer::connection_block block(insert_conn_);
      storage_.insert(position, first, last);
      sync_.insert(index, n);
      DIALS_ASSERT(storage_.size() == sync_.size());
    }

    iterator erase(iterator first, iterator last) {
      size_type first_index = std::distance(begin(), first);
      size_type last_index = std::distance(begin(), last);
      column_synchronizer::connection_block block(erase_conn_);
      iterator result = storage_.erase(first, last);
      sync_.erase(first_index, last_index);
      DIALS_ASSERT(storage_.size() == sync_.size());
      return result;
    }

    iterator erase(iterator position) {
      size_type first_index = std::distance(begin(), position);
      size_type last_index = first_index + 1;
      column_synchronizer::connection_block block(erase_conn_);
      iterator result = storage_.erase(position);
      sync_.erase(first_index, last_index);
      DIALS_ASSERT(storage_.size() == sync_.size());
      return result;
    }

  private:

    struct resizer {
      storage_type *storage;
      resizer(storage_type *s)
        : storage(s) {}
      void operator()(size_type n) {
        storage->resize(n);
      }
    };

    struct inserter {
      storage_type *storage;
      inserter(storage_type *s)
        : storage(s) {}
      void operator()(size_type pos, size_type n) {
        storage->insert(storage->begin() + pos, n, value_type());
      }
    };

    struct eraser {
      storage_type *storage;
      eraser(storage_type *s)
        : storage(s) {}
      void operator()(size_type first, size_type last) {
        storage->erase(storage->begin()+first, storage->begin()+last);
      }
    };

    storage_type storage_;
    column_synchronizer sync_;
    column_synchronizer::connection resize_conn_;
    column_synchronizer::connection insert_conn_;
    column_synchronizer::connection erase_conn_;
  };


  template <typename T>
  class column_data {
  public:

    typedef column_data_core<T> core_type;
    typedef typename core_type::storage_type storage_type;
    typedef typename core_type::value_type value_type;
    typedef typename core_type::reference reference;
    typedef typename core_type::iterator iterator;
    typedef typename core_type::const_iterator const_iterator;
    typedef typename core_type::size_type size_type;
    typedef typename core_type::difference_type difference_type;

    column_data()
      : core_(boost::make_shared<core_type>()) {}

    template <typename InputIterator>
    column_data(InputIterator first, InputIterator last)
      : core_(boost::make_shared<core_type>(first, last)) {}

    column_data(column_synchronizer sync)
      : core_(boost::make_shared<core_type>(sync)) {}

    reference operator[](size_type index) {
      return (*core_)[index];
    }

    const_iterator begin() const {
      return core_->begin();
    }

    const_iterator end() const {
      return core_->end();
    }

    iterator begin() {
      return core_->begin();
    }

    iterator end() {
      return core_->end();
    }

    size_type size() const {
      return core_->size();
    }

    void resize(size_type n) {
      core_->resize(n);
    }

    void push_back(const value_type &v) {
      core_->push_back(v);
    }

    iterator insert(iterator position, const value_type &v) {
      return core_->insert(position, v);
    }

    void insert(iterator position, size_type n, const value_type &v) {
      core_->insert(position, n, v);
    }

    template <typename InputIterator>
    void insert(iterator pos, InputIterator first, InputIterator last) {
      core_->insert(pos, first, last);
    }

    iterator erase(iterator position) {
      return core_->erase(position);
    }

    iterator erase(iterator first, iterator last) {
      return core_->erase(first, last);
    }

  private:

    boost::shared_ptr<core_type> core_;
  };


  struct null_type {};

  template <typename T>
  struct is_null_type : public boost::mpl::bool_<false> {};

  template<>
  struct is_null_type<null_type> : public boost::mpl::bool_<true> {};


  template <typename T0,
            typename T1=null_type, typename T2=null_type, typename T3=null_type,
            typename T4=null_type, typename T5=null_type, typename T6=null_type,
            typename T7=null_type, typename T8=null_type, typename T9=null_type>
  class column_type_generator {
  private:

    template <typename T>
    struct create_column_type {
      typedef column_data<T> type;
    };

    // MPL List of all input types
    typedef boost::mpl::list<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9> all_types;

    // Remove any types if they are null
    typedef typename boost::mpl::remove_if<
      all_types,
      typename boost::mpl::lambda<
        is_null_type<boost::mpl::_1>
      >::type
    >::type valid_types;

    // Create a list of column_data<T> types
    typedef typename boost::mpl::transform<
      valid_types,
      typename boost::mpl::lambda<
        create_column_type<boost::mpl::_1>
      >::type
    >::type column_types;

  public:

    // Expose the variant type
    typedef typename boost::make_variant_over<column_types>::type type;
  };

  template <typename VarientType>
  class column_table {
  public:

    typedef std::map<std::string, VarientType> map_type;
    typedef typename map_type::key_type key_type;
    typedef typename map_type::mapped_type mapped_type;
    typedef typename map_type::value_type map_value_type;
    typedef typename map_type::iterator iterator;
    typedef typename map_type::const_iterator const_iterator;
    typedef typename map_type::size_type size_type;

  private:

    struct proxy {
      column_table *t_;
      key_type k_;

      proxy(column_table *t, key_type k)
        : t_(t), k_(k) {}

      template <typename T>
      operator column_data<T>() const{
        map_type& table = t_->table_;
        iterator it = table.lower_bound(k_);
        if (it == table.end() || table.key_comp()(k_, it->first)) {
          it = table.insert(it, map_value_type(k_,
            mapped_type(column_data<T>(t_->sync_))));
        }
        return boost::get< column_data<T> >(it->second);
      }

      operator mapped_type() const {
        return t_->get(k_);
      }
    };

  public:

    column_table()
      : sync_() {}

    proxy operator[](const key_type &key) {
      return proxy(this, key);
    }

    iterator begin() {
      return table_.begin();
    }

    iterator end() {
      return table_.end();
    }

    const_iterator begin() const {
      return table_.begin();
    }

    const_iterator end() const {
      return table_.end();
    }

    size_type nrows() const {
      return sync_.size();
    }

    size_type ncols() const {
      return table_.size();
    }

    size_type erase(const key_type &key) {
      return table_.erase(key);
    }

    bool empty() const {
      return table_.empty();
    }

    size_type size() const {
      return ncols() ;
    }

    void clear() {
      table_.clear();
    }

    size_type count(const key_type &key) const {
      return table_.count(key);
    }

    iterator find(const key_type &key) {
      return table_.find(key);
    }

    const_iterator find(const key_type &key) const {
      return table_.find(key);
    }

  private:

    mapped_type& get(const key_type &key) {
            return table_[key];
    }

    map_type table_;
    column_synchronizer sync_;
  };

}} // namespace dials::framework

#endif // DIALS_FRAMEWORK_TABLE_COLUMN_TABLE_H
