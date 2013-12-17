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

  /**
   * Base class for synchronizing the sizes of the different columns
   */
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

    /** Initialise the size to zero */
    column_synchronizer_core()
      : size_(0) {}

    /**
     * Initialise the columns to a given size
     * @param n The size of the columns
     */
    column_synchronizer_core(size_type n)
      : size_(n) {}

    /**
     * Resize all the columns to the given size
     * @param n The new size to make the columns.
     */
    void resize(size_type n) {
      size_ = n;
      resize_(n);
    }

    /**
     * Insert elements into each column
     * @param pos The position to insert at
     * @param n The number of elements to insert
     */
    void insert(size_type pos, size_type n) {
      DIALS_ASSERT(pos <= size_);
      size_ += n;
      insert_(pos, n);
    }

    /**
     * Erase some elements from the columns
     * @param pos The position to erase at
     * @param n The number of elements to erase
     */
    void erase(size_type pos, size_type n) {
      DIALS_ASSERT(pos + n <= size_);
      size_ -= n;
      erase_(pos, n);
    }

    /** @returns The size of the columns */
    size_type size() const {
      return size_;
    }

    /**
     * Connect a slot to the resize signal
     * @param slow The slow to connect.
     * @returns The connection
     */
    connection connect_resize(resize_slot_type slot) {
      return resize_.connect(slot);
    }

    /**
     * Connect a slot to the insert signal
     * @param slow The slow to connect.
     * @returns The connection
     */
    connection connect_insert(insert_slot_type slot) {
      return insert_.connect(slot);
    }

    /**
     * Connect a slot to the erase signal
     * @param slow The slow to connect.
     * @returns The connection
     */
    connection connect_erase(erase_slot_type slot) {
      return erase_.connect(slot);
    }

  protected:
    size_type size_;
    resize_signal_type resize_;
    insert_signal_type insert_;
    erase_signal_type erase_;
  };

  /**
   * Class for synchronizing the sizes of the different columns. This class
   * contains a shared pointer to the column_syncrhonizer_core class to ensure
   * that all the columns see the same object.
   */
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

    void erase(size_type pos, size_type n) {
      core_->erase(pos, n);
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


  /**
   * The base class for representing a column of data. The data is contained
   * in contiguous memory and the size of the internal storage is modified
   * automatically whenever another column with the same synchronization
   * object is modified.
   */
  template <typename T>
  class column_data_core {
    struct resizer;
    struct inserter;
    struct eraser;

  public:

    typedef std::vector<T> storage_type;
    typedef typename storage_type::value_type value_type;
    typedef typename storage_type::reference reference;
    typedef typename storage_type::const_reference const_reference;
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

    /**
     * Initialise the column data class with the synchronization object.
     * @param sync The syncrhonization object.
     */
    column_data_core(column_synchronizer sync)
      : storage_(sync.size()),
        sync_(sync),
        resize_conn_(sync_.connect_resize(resizer(&storage_))),
        insert_conn_(sync_.connect_insert(inserter(&storage_))),
        erase_conn_(sync_.connect_erase(eraser(&storage_))) {}

    /**
     * Access an element at the given index.
     * @param n The index of the element
     * @returns A reference to the element
     */
    reference operator[](size_type n) {
      return storage_[n];
    }

    /**
     * Access an element at the given index.
     * @param n The index of the element
     * @returns A const reference to the element
     */
    const_reference operator[](size_type n) const {
      return storage_[n];
    }

    /** @returns A const reference to the beginning of the column */
    const_iterator begin() const {
      return storage_.begin();
    }

    /** @returns A const reference to the end of the column */
    const_iterator end() const {
      return storage_.end();
    }

    /** @returns A reference to the beginning of the column */
    iterator begin() {
      return storage_.begin();
    }

    /** @returns A reference to the end of the column */
    iterator end() {
      return storage_.end();
    }

    /** @returns The size of the column */
    size_type size() const {
      return sync_.size();
    }

    /**
     * Resize the column
     * @param n The size to resize to.
     */
    void resize(size_type n) {
      sync_.resize(n);
    }

    /**
     * Append an element to the column
     * @param v The value of the element.
     */
    void push_back(const value_type &v) {
      insert(end(), v);
    }

    /**
     * Insert an element at the given position
     * @param position The position to insert
     * @param v The value of the element to insert
     * @returns An iterator to the new element
     */
    iterator insert(iterator position, const value_type &v) {
      size_type pos = std::distance(begin(), position);
      column_synchronizer::connection_block block(insert_conn_);
      iterator result = storage_.insert(position, v);
      sync_.insert(pos, 1);
      DIALS_ASSERT(storage_.size() == sync_.size());
      return result;
    }

    /**
     * Insert multiple elements at the given position
     * @param position The position to begin inserting
     * @param n The number of elements to insert
     * @param v The value of the element to insert
     */
    void insert(iterator position, size_type n, const value_type &v) {
      size_type pos = std::distance(begin(), position);
      column_synchronizer::connection_block block(insert_conn_);
      storage_.insert(position, n, v);
      sync_.insert(pos, n);
      DIALS_ASSERT(storage_.size() == sync_.size());
    }

    /**
     * Insert multiple elements at the given position
     * @param position The position to begin inserting
     * @param first The beginning of the sequence
     * @param last The end of the sequence
     */
    template <typename InputIterator>
    void insert(iterator position, InputIterator first, InputIterator last) {
      size_type n = std::distance(first, last);
      size_type pos = std::distance(begin(), position);
      column_synchronizer::connection_block block(insert_conn_);
      storage_.insert(position, first, last);
      sync_.insert(pos, n);
      DIALS_ASSERT(storage_.size() == sync_.size());
    }

    /**
     * Erase multiple elements from the column
     * @param first The beginning of the sequence
     * @param last The end of the sequence
     * @returns The iterator to the element following the last erased element
     */
    iterator erase(iterator first, iterator last) {
      size_type pos = std::distance(begin(), first);
      size_type n = std::distance(first, last);
      column_synchronizer::connection_block block(erase_conn_);
      iterator result = storage_.erase(first, last);
      sync_.erase(pos, n);
      DIALS_ASSERT(storage_.size() == sync_.size());
      return result;
    }

    /**
     * Erase an element from the column
     * @param position The beginning of the sequence
     * @returns The iterator to the element following the erased element
     */
    iterator erase(iterator position) {
      size_type pos = std::distance(begin(), position);
      size_type n = 1;
      column_synchronizer::connection_block block(erase_conn_);
      iterator result = storage_.erase(position);
      sync_.erase(pos, n);
      DIALS_ASSERT(storage_.size() == sync_.size());
      return result;
    }

  private:

    /**
     * A functor handling resizing the column on a signal from the synchronizer
     * class. This ensures that if another column is resized, this column
     * is resized to match it.
     */
    struct resizer {
      storage_type *storage;
      resizer(storage_type *s)
        : storage(s) {}
      void operator()(size_type n) {
        storage->resize(n);
      }
    };

    /**
     * A functor handling insertion in the column on a signal from the
     * synchronizer class. When elements are inserted into another column,
     * a default element is inserted into this column to ensure both columns
     * remain the same size.
     */
    struct inserter {
      storage_type *storage;
      inserter(storage_type *s)
        : storage(s) {}
      void operator()(size_type pos, size_type n) {
        storage->insert(storage->begin() + pos, n, value_type());
      }
    };

    /**
     * A functor handling erasing from the column on a signal from the
     * synchronizer class. When elements are erased from another column,
     * the corresponding element is erased from this column to ensure both
     * columns remain the same size.
     */
    struct eraser {
      storage_type *storage;
      eraser(storage_type *s)
        : storage(s) {}
      void operator()(size_type pos, size_type n) {
        iterator first = storage->begin()+pos;
        iterator last = first + n;
        storage->erase(first, last);
      }
    };

    storage_type storage_;
    column_synchronizer sync_;
    column_synchronizer::connection resize_conn_;
    column_synchronizer::connection insert_conn_;
    column_synchronizer::connection erase_conn_;
  };


  /**
   * Class representing a column of data. This class contains a shared pointer
   * to the column_data_core class. It also provides a thin wrapper to the
   * column_data_core interface.
   */
  template <typename T>
  class column_data {
  public:

    typedef column_data_core<T> core_type;
    typedef typename core_type::storage_type storage_type;
    typedef typename core_type::value_type value_type;
    typedef typename core_type::reference reference;
    typedef typename core_type::const_reference const_reference;
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

    reference operator[](size_type n) {
      return (*core_)[n];
    }

    const_reference operator[](size_type n) const {
      return (*core_)[n];
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


  /**
   * A class to represent a column-centric table. I.e. a table in which the
   * data is represented as a list of columns. It is created with a variant
   * type of column data. It can be instantiated as follows:
   *
   * typedef column_type_generator<int, double>::type column_types;
   * column_table<column_types> table;
   *
   * The columns can be accessed as values in a std::map as follows:
   *
   * column_data<int> col = table["column"];
   */
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
      operator column_data<T>() const {
        boost::shared_ptr<map_type> table = t_->table_;
        iterator it = table->lower_bound(k_);
        if (it == table->end() || table->key_comp()(k_, it->first)) {
          it = table->insert(it, map_value_type(k_,
            mapped_type(column_data<T>(t_->sync_))));
        }
        return boost::get< column_data<T> >(it->second);
      }

      operator mapped_type() const {
        return t_->get(k_);
      }
    };

  public:

    /** Initialise the table */
    column_table()
      : table_(boost::make_shared<map_type>()),
        sync_() {}

    /**
     * Initialise the table to a certain size
     * @param n The size to initialise to
     */
    column_table(size_type n)
      : table_(boost::make_shared<map_type>()),
        sync_(n) {}

    /**
     * Access a column by key
     * @param key The column name
     * @returns The proxy object to access the column data
     */
    proxy operator[](const key_type &key) {
      return proxy(this, key);
    }

    /** @returns An iterator to the beginning of the column map */
    iterator begin() {
      return table_->begin();
    }

    /** @returns An iterator to the end of the column map */
    iterator end() {
      return table_->end();
    }

    /** @returns A const iterator to the beginning of the column map */
    const_iterator begin() const {
      return table_->begin();
    }

    /** @returns A const iterator to the end of the column map */
    const_iterator end() const {
      return table_->end();
    }

    /** @returns The number of rows in the table */
    size_type nrows() const {
      return sync_.size();
    }

    /** @returns The number of columns in the table */
    size_type ncols() const {
      return table_->size();
    }

    /**
     * Erase a column from the table.
     * @param key The column name
     * @returns The number of columns removed
     */
    size_type erase(const key_type &key) {
      return table_->erase(key);
    }

    /** @returns Is the table empty */
    bool empty() const {
      return table_->empty();
    }

    /** @returns The number of columns in the table */
    size_type size() const {
      return ncols() ;
    }

    /** Clear the table */
    void clear() {
      table_->clear();
    }

    /** @returns The number of columns matching the key (0 or 1) */
    size_type count(const key_type &key) const {
      return table_->count(key);
    }

    /**
     * Find a column matching the key
     * @param key The column name
     * @returns An iterator to the column
     */
    iterator find(const key_type &key) {
      return table_->find(key);
    }

    /**
     * Find a column matching the key
     * @param key The column name
     * @returns A const iterator to the column
     */
    const_iterator find(const key_type &key) const {
      return table_->find(key);
    }

    /**
     * Resize the columns to the given size
     * @param n The size to resize to
     */
    void resize(size_type n) {
      sync_.resize(n);
    }

    /**
     * Insert an element at the given position in each column
     * @param pos The position to insert at
     */
    void insert(size_type pos) {
      insert(pos, 1);
    }

    /**
     * Insert a number of elements at the given position in each column
     * @param pos The position to insert at
     * @param n The number of elements to insert
     */
    void insert(size_type pos, size_type n) {
      sync_.insert(pos, n);
    }

    /**
     * Erase an element from all the columns
     * @param pos The position of the element
     */
    void erase(size_type pos) {
      erase(pos, 1);
    }

    /**
     * Erase a number of elements from all the columns
     * @param pos The position of the element
     * @param n The number of elements to erase
     */
    void erase(size_type pos, size_type n) {
      sync_.erase(pos, n);
    }

  private:

    /** @returns The element at the given key */
    mapped_type& get(const key_type &key) {
      return (*table_)[key];
    }

    boost::shared_ptr<map_type> table_;
    column_synchronizer sync_;
  };


  struct null_type {};

  template <typename T>
  struct is_null_type : public boost::mpl::bool_<false> {};

  template<>
  struct is_null_type<null_type> : public boost::mpl::bool_<true> {};


  /**
   * A class to help in generating a variant type for use in the column_table
   * class. Given a variadic template list of types it will convert the list
   * into list of column_data<T> types and then export a variant.
   *
   * For example:
   *  column_type_generator<int, double, std::string>::type
   *
   * Becomes:
   *
   *  variant<
   *    column_data<int>,
   *    column_data<double>,
   *    column_data<std::string>
   *  >
   */
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
      typename boost::mpl::lambda< is_null_type<boost::mpl::_1> >::type
    >::type valid_types;

    // Create a list of column_data<T> types
    typedef typename boost::mpl::transform<
      valid_types,
      typename boost::mpl::lambda< create_column_type<boost::mpl::_1> >::type
    >::type column_types;

  public:

    // Expose the variant type
    typedef typename boost::make_variant_over<column_types>::type type;
  };

}} // namespace dials::framework

#endif // DIALS_FRAMEWORK_TABLE_COLUMN_TABLE_H
