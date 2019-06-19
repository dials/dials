/*
 * qotree.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPATIAL_INDEXING_QOTREE_H
#define DIALS_ALGORITHMS_SPATIAL_INDEXING_QOTREE_H

#include <cstdlib>
#include <deque>
#include <list>
#include <cassert>
#include <stack>

namespace dials { namespace algorithms {

  /**
   * Template function to be specialized, calculate max depth
   * @tparam BoxType The box type
   * @return The maximum depth of the tree
   */
  template <typename BoxType>
  std::size_t maximum_depth(const BoxType &);

  /**
   * Template function to be specialized, subdivide the box.
   * @tparam BoxType The box type
   * @return The sub-divided box
   */
  template <typename BoxType>
  BoxType subdivide(const BoxType &, std::size_t);

  /**
   * Template structure containing functions to use in comparison of objects
   * @tparam ObjectTypeA The bounding box type
   * @tparam ObjectTypeB The object type
   */
  template <typename ObjectTypeA, typename ObjectTypeB>
  struct compare {
    /**
     * @return True/False object a contains object b
     */
    static bool contains(const ObjectTypeA &, const ObjectTypeB &);

    /**
     * @return True/False the objects collide
     */
    static bool collides(const ObjectTypeA &, const ObjectTypeB &);
  };

  /**
   * A class that implements the base behaviour for a quadtree and octree.
   * @tparam NDIM The dimensions (i.e. quad or oct)
   * @tparam BoxType The bounding box type to use
   * @tparam ObjectType The type of object to store
   */
  template <int NDIM, typename BoxType, typename ObjectType>
  class QOTree {
  public:
    // The number of children 2^NDIM
    static const int NCHILD = (1 << NDIM);

    // General typedefs
    typedef BoxType box_type;
    typedef ObjectType object_type;
    typedef std::size_t size_type;

    // Object list typedefs
    typedef std::list<object_type> object_list_type;
    typedef typename object_list_type::iterator object_iterator;
    typedef typename object_list_type::const_iterator const_object_iterator;

    // The node structure
    struct node_type {
      box_type box;              // The bounding box
      node_type *parent;         // The parent node
      node_type *child[NCHILD];  // The child nodes
      size_type bucket_size;     // The bucket size
      object_list_type bucket;   // The bucket list
      bool is_leaf;              // Is the node a leaf

      // Default node initialisation
      node_type() : parent(NULL), bucket_size(0), is_leaf(true) {}

      // Initialise the node with a box and parent
      node_type(const box_type &b, node_type *p)
          : box(b), parent(p), bucket_size(0), is_leaf(true) {}
    };

    // Node list typedefs
    typedef std::deque<node_type> node_list_type;
    typedef typename node_list_type::pointer node_pointer;
    typedef typename node_list_type::const_pointer const_node_pointer;
    typedef typename node_list_type::iterator node_iterator;
    typedef typename node_list_type::const_iterator const_node_iterator;

    /** Create the tree */
    QOTree(const box_type &box, size_type max_bucket_size = 10) {
      push_node(node_type(box, NULL));
      max_bucket_size_ = max_bucket_size;
      max_depth_ = maximum_depth(box);
    }

    /** Get the const iterator to the beginning of the node list */
    const_node_iterator node_begin() const {
      return node_list_.begin();
    }

    /** Get the const iterator to the end of the node list */
    const_node_iterator node_end() const {
      return node_list_.end();
    }

    /** Get the maximum bucket size */
    size_type max_bucket_size() const {
      return max_bucket_size_;
    }

    /** Get the maximum depth of the tree */
    size_type max_depth() const {
      return max_depth_;
    }

    /** Insert an element into the tree */
    bool insert(const object_type &v) {
      // Get the root node and ensure input is inside
      node_pointer root = &node_list_.front();
      if (compare<box_type, object_type>::contains(root->box, v)) {
        return insert(v, root);
      }

      // Unsuccessfully inserted
      return false;
    }

    /** Find all the objects in the given range */
    template <typename ObjectList>
    bool query_range(const box_type &range, ObjectList &elements) const {
      // Get the root node and ensure the range is contained within.
      // Then query the tree to get the list of elements that are contained
      // with the given range.
      const_node_pointer root = &node_list_.front();
      if (compare<box_type, box_type>::collides(root->box, range)) {
        return query_range(range, elements, root);
      }

      // Unsuccessful query
      return false;
    }

    /** Find all the objects colliding with the given object */
    template <typename ObjectList>
    bool query_collision(const object_type &v, ObjectList &elements) const {
      // Get the root node and ensure the object is contained within.
      // Then query the tree to get the list of elements that collide.
      const_node_pointer root = &node_list_.front();
      if (compare<box_type, object_type>::collides(root->box, v)) {
        return query_collision(v, elements, root);
      }

      // Unsuccessful query
      return false;
    }

  protected:
    /** Push the node and return a pointer to it */
    node_pointer push_node(const node_type &node) {
      node_list_.push_back(node);
      return &node_list_.back();
    }

    /** Choose the sub-division to that the object is contained within */
    template <typename T>
    size_type choose_subdivision(const T &v, const_node_pointer node) const {
      // Find the sub-division in which the element is contained. Return the
      // index of the sub-division. If no subdivision fully contains the
      // element, this will return NCHILD.
      size_type div;
      for (div = 0; div < NCHILD; ++div) {
        if (compare<box_type, T>::contains(node->child[div]->box, v)) {
          break;
        }
      }
      return div;
    }

    /** Move elements from the bucket of one node to another */
    void move_element(node_pointer from, node_pointer to, object_iterator it) {
      // Move the element from list to list without moving the element itself
      // in memory, increment the 'to' counter and decrement the 'from' counter.
      to->bucket.splice(to->bucket.end(), from->bucket, it);
      to->bucket_size++;
      from->bucket_size--;
    }

    /** Distribute all the elements in a nodes bucket to its children */
    void distribute_bucket(node_pointer node) {
      // Loop through all the elements in the node list. If the element is fully
      // contained in any of the child nodes, then transfer the element into the
      // child list, removing it from the parent list.
      object_iterator curr, next = node->bucket.begin();
      for (curr = next++; curr != node->bucket.end(); curr = next++) {
        size_type div;
        if ((div = choose_subdivision(*curr, node)) < NCHILD) {
          move_element(node, node->child[div], curr);
        }
      }
    }

    /** Sub-divide the node */
    void split_node(node_pointer node) {
      // Create all the child nodes
      node->is_leaf = false;
      for (std::size_t i = 0; i < NCHILD; ++i) {
        node->child[i] = push_node(node_type(subdivide(node->box, i), node));
      }

      // Distribute all parent elements onto the child nodes
      distribute_bucket(node);
    }

    /** Check the bucket size is < the max. If it isn't then split the node */
    bool ensure_bucket_size(node_pointer node) {
      // If the size of the bucket is >= max bucket size split the node
      if (node->bucket_size >= max_bucket_size_) {
        split_node(node);
        return false;
      }

      // No split
      return true;
    }

    /** Insert an element into the tree */
    bool insert(const object_type &v, node_pointer node) {
      // Iterate until we reach the maximum depth
      for (size_type depth = 0; depth < max_depth_; ++depth) {
        // If this is a leaf node, then check if the bucket has space for
        // this element. If so add it and return, otherwise split the node
        // into four quadrants and iterate again.
        if (node->is_leaf && ensure_bucket_size(node)) {
          break;
        }

        // If this item is not a leaf node, then check in which of the boxes
        // the new item is contained. If none of the child boxes fully contains
        // the item, add it to the list for the current node and return.
        size_type div = choose_subdivision(v, node);
        if (div < NCHILD) {
          node = node->child[div];
        } else {
          break;
        }
      }

      // Insert element into node bucket
      node->bucket.push_back(v);
      node->bucket_size++;

      // Successfully inserted
      return true;
    }

    /** Add objects contained in range in the current node to the list */
    template <typename ObjectList>
    void append_contained_objects(const box_type &range,
                                  ObjectList &elements,
                                  const_node_pointer node) const {
      // Loop through all objects in the node and add any that are contained
      // to the given object list.
      for (const_object_iterator it = node->bucket.begin(); it != node->bucket.end();
           ++it) {
        if (compare<box_type, object_type>::contains(range, *it)) {
          elements.push_back(*it);
        }
      }
    }

    /** Find all the objects in the given range */
    template <typename ObjectList>
    bool query_range(const box_type &range,
                     ObjectList &elements,
                     const_node_pointer node) const {
      // Find the node that fully contains the range.
      for (;;) {
        if (node->is_leaf) {
          break;
        }
        size_type div = choose_subdivision(range, node);
        if (div < NCHILD) {
          node = node->child[div];
        } else {
          break;
        }
      }

      // Create a stack and push the current node. Loop until the stack is
      // empty. Get the node from the stack and add any contained elements in
      // the node to the element list. Then, if the current node is not a leaf
      // node, search through all the children and push any that collide with
      // the object to the stack.
      std::stack<const_node_pointer> stack;
      stack.push(node);
      while (!stack.empty()) {
        const_node_pointer node = stack.top();
        stack.pop();
        append_contained_objects(range, elements, node);
        if (!node->is_leaf) {
          for (size_type i = 0; i < NCHILD; ++i) {
            if (compare<box_type, box_type>::collides(range, node->child[i]->box)) {
              stack.push(node->child[i]);
            }
          }
        }
      }

      // Successful query
      return true;
    }

    /** Add objects colliding with v in the current node to the list */
    template <typename ObjectList>
    void append_colliding_objects(const object_type &v,
                                  ObjectList &elements,
                                  const_node_pointer node) const {
      // Loop through all objects in the node and add any that collide
      // to the given object list.
      for (const_object_iterator it = node->bucket.begin(); it != node->bucket.end();
           ++it) {
        if (compare<object_type, object_type>::collides(*it, v)) {
          elements.push_back(*it);
        }
      }
    }

    /** Find all the objects colliding with the given object */
    template <typename ObjectList>
    bool query_collision(const object_type &v,
                         ObjectList &elements,
                         const_node_pointer node) const {
      // Create a stack and push the root node. Loop until the stack is empty.
      // Get the node from the stack and add any colliding elements in the node
      // to the element list. Then, if the current node is not a leaf node,
      // search through all the children and push any that collide with the
      // object to the stack.
      std::stack<const_node_pointer> stack;
      stack.push(node);
      while (!stack.empty()) {
        const_node_pointer node = stack.top();
        stack.pop();
        append_colliding_objects(v, elements, node);
        if (!node->is_leaf) {
          for (size_type i = 0; i < NCHILD; ++i) {
            if (compare<box_type, object_type>::collides(node->child[i]->box, v)) {
              stack.push(node->child[i]);
            }
          }
        }
      }

      // Successful query
      return true;
    }

    node_list_type node_list_;
    size_type max_bucket_size_;
    size_type max_depth_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPATIAL_INDEXING_QOTREE_H
