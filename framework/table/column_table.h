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
#include <dials/array_family/scitbx_shared_and_versa.h>

namespace dials { namespace framework {

  namespace bs2 = boost::signals2;

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

    /**
     * operator[] proxy to aid in returning and casting elements.
     */
    struct proxy {
      column_table *t_;
      key_type k_;

      proxy(column_table *t, key_type k)
        : t_(t), k_(k) {}

      /**
       * Cast the element to the desired column data type. If no element is
       * present, a new element with the desired type is created and returned.
       * Otherwise the current element is returned. If the types don't match,
       * an exception is raised.
       */
      template <typename T>
      operator af::shared<T>() const {
        size_type n = t_->sync_.size();
        boost::shared_ptr<map_type> table = t_->table_;
        iterator it = table->lower_bound(k_);
        if (it == table->end() || table->key_comp()(k_, it->first)) {
          it = table->insert(it, map_value_type(k_,
            mapped_type(af::shared<T>(n))));
        }
        return boost::get< af::shared<T> >(it->second);
      }

      /**
       * Return the mapped variant type at the given element directly.
       */
      operator mapped_type() const {
        return t_->get(k_);
      }
    };


    class column_synchronizer {
    public:

      typedef std::size_t size_type;

      struct size_visitor : boost::static_visitor<size_type> {

        template <typename T>
        size_type operator()(const T &v) const {
          return v.size();
        }
      };

      struct resize_visitor : boost::static_visitor<void> {
        size_type n_;
        resize_visitor(size_type n)
          : n_(n) {}
        template <typename T>
        void operator()(T &v) {
          v.resize(n_);
        }
      };

      struct insert_visitor : boost::static_visitor<void> {
        size_type pos, n;
        insert_visitor(size_type pos_, size_type n_)
          : pos(pos_), n(n_) {}
        template <typename T>
        void operator()(T &v) {
          v.insert(v.begin() + pos, n, typename T::value_type());
        }
      };

      struct erase_visitor : boost::static_visitor<void> {
        size_type pos, n;
        erase_visitor(size_type pos_, size_type n_)
          : pos(pos_), n(n_) {}
        template <typename T>
        void operator()(T &v) {
          iterator first = v.begin() + pos;
          iterator last = first + n;
          v.erase(first, last);
        }
      };

      /** Initialise the size to zero */
      column_synchronizer(column_table *t)
        : size_(0),
          table_(t) {}

      /**
       * Initialise the columns to a given size
       * @param n The size of the columns
       */
      column_synchronizer(size_type n, column_table *t)
        : size_(n),
          table_(t) {}

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

      /* @returns The size of the columns */
      size_type size() const {
        DIALS_ASSERT(is_consistent());
        return size_;
      }

      bool is_consistent() const {
        size_visitor visitor;
        for (iterator it = table_->begin(); it != table_->end(); ++it) {
          if (it->second.apply_visitor(visitor) != size_) {
            return false;
          }
        }
        return true;
      }

    protected:

      void resize_(size_type n) {
        resize_visitor visitor(n);
        for (iterator it = table_->begin(); it != table_->end(); ++it) {
          it->second.apply_visitor(visitor);
        }
      }

      void insert_(size_type pos, size_type n) {
        insert_visitor visitor(pos, n);
        for (iterator it = table_->begin(); it != table_->end(); ++it) {
          it->second.apply_visitor(visitor);
        }
      }

      void erase_(size_type pos, size_type n) {
        erase_visitor visitor(pos, n);
        for (iterator it = table_->begin(); it != table_->end(); ++it) {
          it->second.apply_visitor(visitor);
        }
      }

      size_type size_;
      column_table *table_;
    };

  public:

    /** Initialise the table */
    column_table()
      : table_(boost::make_shared<map_type>()),
        sync_(this) {}

    /**
     * Initialise the table to a certain size
     * @param n The size to initialise to
     */
    column_table(size_type n)
      : table_(boost::make_shared<map_type>()),
        sync_(n, this) {}

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
      return nrows() ;
    }

    bool is_consistent() const {
      return sync_.is_consistent();
    }

    /** Clear the table */
    void clear() {
      table_->clear();
      sync_.resize(0);
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
   * into list of af::shared<T> types and then export a variant.
   *
   * For example:
   *  column_type_generator<int, double, std::string>::type
   *
   * Becomes:
   *
   *  variant<
   *    af::shared<int>,
   *    af::shared<double>,
   *    af::shared<std::string>
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
      typedef af::shared<T> type;
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
