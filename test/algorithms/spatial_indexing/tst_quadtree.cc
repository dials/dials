
#include <iostream>
#include <string>
#include <dials/algorithms/spatial_indexing/quadtree.h>

using namespace dials::algorithms;

typedef Quadtree<Box> MyQuadtree;

/** Print out a box */
std::ostream &operator<<(std::ostream &os, const Box &box) {
  os << "(" << box.x0 << ", " << box.y0 << ", " << box.x1 << ", " << box.y1 << ")";
  return os;
}

/** The quadtree node bucket */
void print_bucket(const MyQuadtree::object_list_type &l, const std::string &indent) {
  MyQuadtree::const_object_iterator it;
  for (it = l.begin(); it != l.end(); ++it) {
    std::cout << indent << *it << "\n";
  }
}

/** Print the tree from the current node */
void print_tree(MyQuadtree::const_node_pointer node, const std::string &indent) {
  if (!node) return;
  std::cout << indent << "Box: " << node->box << "\n";
  std::cout << indent << "Nobj: " << node->bucket_size << "\n";
  print_bucket(node->bucket, indent);
  print_tree(node->child[0], indent + "  ");
  print_tree(node->child[1], indent + "  ");
  print_tree(node->child[2], indent + "  ");
  print_tree(node->child[3], indent + "  ");
}

/** Print the quadtree from the root */
void print_tree(MyQuadtree &tree) {
  print_tree(&*tree.node_begin(), "");
}

/** Generate a random box */
Box random_box(Box bounds, int min_size, int max_size) {
  int x = bounds.x0 + rand() % (bounds.x1 - bounds.x0);
  int y = bounds.y0 + rand() % (bounds.y1 - bounds.y0);
  int w = min_size + rand() % (max_size - min_size);
  int h = min_size + rand() % (max_size - min_size);
  return Box(x, y, x + w, y + h);
}

/** Test that no quadtree boxes have size < 1 */
bool test_node_box_sizes(const MyQuadtree &tree) {
  // Loop through all nodes and check the box size
  for (MyQuadtree::const_node_iterator it = tree.node_begin(); it != tree.node_end();
       ++it) {
    // If x or y size is less than 1 then fail the test
    if ((it->box.x1 - it->box.x0) < 1 || (it->box.y1 - it->box.y0) < 1) {
      return false;
    }
  }

  // Test passed
  return true;
}

/** Test that all the elements have been put in the correct box */
bool test_elements_in_correct_box(const MyQuadtree &tree) {
  typedef compare<Box, Box> check;

  // Loop through all the nodes in the tree
  for (MyQuadtree::const_node_iterator it = tree.node_begin(); it != tree.node_end();
       ++it) {
    // Loop through all the items contained in the node
    std::size_t n_obj = 0;
    for (MyQuadtree::const_object_iterator vit = it->bucket.begin();
         vit != it->bucket.end();
         ++vit) {
      // If the item does not belong in this node then fail the test
      if (!check::contains(it->box, *vit)) {
        return false;
      }

      // If this is not a leaf node and the item could belong in one of the
      // nodes children then fail the test
      if (!it->is_leaf) {
        for (std::size_t i = 0; i < 4; ++i) {
          if (check::contains(it->child[i]->box, *vit)) {
            return false;
          }
        }
      }
      n_obj++;
    }

    // If the number of objects in the list does not match the size of the
    // bucket then fail the test.
    if (n_obj != it->bucket_size) {
      return false;
    }
  }

  // Test passed
  return true;
}

/** Test the total number of elements in the tree */
bool test_num_elements(const MyQuadtree &tree, std::size_t n_obj_expected) {
  // Loop through all the nodes in the tree and sum the number of elements
  // in the whole tree
  std::size_t n_obj = 0;
  for (MyQuadtree::const_node_iterator it = tree.node_begin(); it != tree.node_end();
       ++it) {
    n_obj += it->bucket_size;
  }

  // Check if the number of elements is the same as expected.
  if (n_obj != n_obj_expected) {
    return false;
  }

  // Test passed
  return true;
}

/** By brute-force find the number of elements in the range */
void query_range_brute_force(MyQuadtree &tree,
                             const Box &box,
                             std::size_t &n_elements) {
  typedef compare<Box, Box> check;

  // Init n_elements to zero
  n_elements = 0;

  // Loop through all the tree nodes
  for (MyQuadtree::const_node_iterator node = tree.node_begin();
       node != tree.node_end();
       ++node) {
    // Loop through all the values in the node
    for (MyQuadtree::const_object_iterator value = node->bucket.begin();
         value != node->bucket.end();
         ++value) {
      // If the value is within the box then add to the counter
      if (check::contains(box, *value)) {
        n_elements++;
      }
    }
  }
}

/** By brute-force find the number of elements that collide */
void query_collision_brute_force(MyQuadtree &tree,
                                 const Box &box,
                                 std::size_t &n_elements) {
  typedef compare<Box, Box> check;

  // Init n_elements to zero
  n_elements = 0;

  // Loop through all the tree nodes
  for (MyQuadtree::const_node_iterator node = tree.node_begin();
       node != tree.node_end();
       ++node) {
    // Loop through all the values in the node
    for (MyQuadtree::const_object_iterator value = node->bucket.begin();
         value != node->bucket.end();
         ++value) {
      // If the value collides with the box then add to the counter
      if (check::collides(box, *value)) {
        n_elements++;
      }
    }
  }
}

/** Perform a selection of range queries and check the results */
bool test_query_range(MyQuadtree &tree, const std::list<Box> &range) {
  // Loop through all the range items
  std::list<Box> elements;
  for (std::list<Box>::const_iterator it = range.begin(); it != range.end(); ++it) {
    // Clear elements
    elements.clear();

    // Query the quadtree with a range and get the elements in the range
    if (tree.query_range(*it, elements)) {
      // Find the number of elements by brute forace
      std::size_t n_brute_force_elements;
      query_range_brute_force(tree, *it, n_brute_force_elements);

      // Calculate the number of elements found by the quadtree
      std::size_t n_elements = std::distance(elements.begin(), elements.end());

      // If the number of elements does not match those obtained through
      // brute force then fail the test.
      if (n_elements != n_brute_force_elements) {
        return false;
      }
    }
  }

  // Test passed
  return true;
}

/** Perform a selection of collision queries and check the results */
bool test_query_collides(MyQuadtree &tree, const std::list<Box> &v) {
  std::list<Box> elements;

  // Loop through all the given boxes
  for (std::list<Box>::const_iterator it = v.begin(); it != v.end(); ++it) {
    // Clear the elements
    elements.clear();

    // Query the quad tree for elements that collide with the box
    if (tree.query_collision(*it, elements)) {
      // Find the number of elements by brute forace
      std::size_t n_brute_force_elements;
      query_collision_brute_force(tree, *it, n_brute_force_elements);

      // Calculate the number of elements found
      std::size_t n_elements = std::distance(elements.begin(), elements.end());

      // If the number of elements is not the same as brute force then
      // fail the test.
      if (n_elements != n_brute_force_elements) {
        return false;
      }
    }
  }

  // Test passed
  return true;
}

std::string pass_fail(bool pass) {
  return pass ? "OK" : "Fail";
}

int main(int argc, char const *argv[]) {
  MyQuadtree tree(Box(0, 0, 512, 512));
  // clock_t st = clock();
  std::size_t num = 10000;
  std::size_t num_valid_obj = 0;
  for (std::size_t i = 0; i < num; ++i) {
    if (tree.insert(random_box(Box(0, 0, 512, 512), 3, 8))) {
      num_valid_obj++;
    }
  }
  // std::cout << "Time: " << ((float)(clock() - st)) / CLOCKS_PER_SEC << "\n";

  std::list<Box> range;
  for (std::size_t i = 0; i < 10000; ++i) {
    range.push_back(random_box(Box(0, 0, 512, 512), 10, 100));
  }

  //  std::list <Box> elements;
  //  st = clock();
  //  for (std::list<Box>::iterator it = range.begin(); it != range.end(); ++it) {
  //    tree.query_collision(*it, elements);
  //  }
  //  std::cout << "Time: " << ((float)(clock() - st)) / CLOCKS_PER_SEC << "\n";

  std::cout << pass_fail(test_node_box_sizes(tree)) << std::endl;
  std::cout << pass_fail(test_elements_in_correct_box(tree)) << std::endl;
  std::cout << pass_fail(test_num_elements(tree, num_valid_obj)) << std::endl;
  std::cout << pass_fail(test_query_range(tree, range)) << std::endl;
  std::cout << pass_fail(test_query_collides(tree, range)) << std::endl;

  // print_tree(tree);

  return 0;
}
