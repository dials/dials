/*
 * find_overlapping_reflections.cc
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
#include <boost/python/iterator.hpp>
#include <boost_adaptbx/std_pair_conversion.h>
#include <dials/algorithms/integration/find_overlapping_reflections.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  typedef AdjacencyList::vertex_descriptor vertex_descriptor;
  typedef AdjacencyList::vertex_iterator vertex_iterator;
  typedef AdjacencyList::edge_iterator edge_iterator;
  typedef AdjacencyList::adjacency_iterator adjacency_iterator;
  typedef AdjacencyList::edge_descriptor edge_descriptor;
  typedef std::pair<vertex_iterator, vertex_iterator> vertex_iterator_range;
  typedef std::pair<edge_iterator, edge_iterator> edge_iterator_range;
  typedef std::pair<adjacency_iterator, adjacency_iterator> 
    adjacency_iterator_range;
    
  void adjacency_list_add_edge(boost::shared_ptr<AdjacencyList> &list, 
      vertex_descriptor i, vertex_descriptor j) {
    add_edge(i, j, *list);
  }
  
  void adjacency_list_remove_edge(boost::shared_ptr<AdjacencyList> &list,
      vertex_descriptor i, vertex_descriptor j) {
    remove_edge(i, j, *list);
  }

  void adjacency_list_clear_vertex(boost::shared_ptr<AdjacencyList> &list,
      vertex_descriptor v) {
    clear_vertex(v, *list);
  }

  void adjacency_list_clear(boost::shared_ptr<AdjacencyList> &list) {
    list->clear();
  }

  vertex_iterator_range adjacency_list_vertices(
      boost::shared_ptr<AdjacencyList> &list) {
    return vertices(*list);
  }

  edge_iterator_range adjacency_list_edges(
      boost::shared_ptr<AdjacencyList> &list) {
    return edges(*list);
  }

  adjacency_iterator_range adjacency_list_adjacent_vertices(
      boost::shared_ptr<AdjacencyList> &list, vertex_descriptor v) {
    return adjacent_vertices(v, *list);
  }
  
  std::size_t adjacency_list_num_vertices(
      boost::shared_ptr<AdjacencyList> &list) {
    return num_vertices(*list);
  }
  std::size_t adjacency_list_num_edges(
      boost::shared_ptr<AdjacencyList> &list) {
    return num_edges(*list);
  }
  
  std::pair<vertex_descriptor, vertex_descriptor> adjacency_list_edge_vertices(
      boost::shared_ptr<AdjacencyList> &list, edge_descriptor edge) {
    return std::pair<vertex_descriptor, vertex_descriptor>(
      source(edge, *list), target(edge, *list)); 
  }

  void export_find_overlapping_reflections()
  {
    boost_adaptbx::std_pair_conversions::
      to_and_from_tuple<vertex_descriptor, vertex_descriptor>();
  
    class_<vertex_iterator_range>("VertexIteratorRange", no_init)
      .def("__iter__",
        boost::python::range(
          &vertex_iterator_range::first, 
          &vertex_iterator_range::second));

    class_<edge_descriptor>("EdgeDescriptor", no_init);

    class_<edge_iterator_range>("EdgeIteratorRange", no_init)
      .def("__iter__",
        boost::python::range(
          &edge_iterator_range::first, 
          &edge_iterator_range::second));
        
    class_<adjacency_iterator_range>("AdjacencyIteratorRange", no_init)
      .def("__iter__",
        boost::python::range(
          &adjacency_iterator_range::first, 
          &adjacency_iterator_range::second));        
        
    class_<boost::shared_ptr<AdjacencyList> >("AdjacencyList")
      .def("add_edge", &adjacency_list_add_edge)
      .def("remove_edge", &adjacency_list_remove_edge)
      .def("clear_vertex", &adjacency_list_clear_vertex)
      .def("clear", &adjacency_list_clear)
      .def("num_vertices", &adjacency_list_num_vertices)
      .def("num_edges", &adjacency_list_num_edges)
      .def("vertices", &adjacency_list_vertices)
      .def("edges", &adjacency_list_edges)
      .def("adjacent_vertices", &adjacency_list_adjacent_vertices)
      .def("edge_vertices", &adjacency_list_edge_vertices)
      .def("__len__", &adjacency_list_num_edges);
      
    def("find_overlapping_reflections", 
      &find_overlapping_reflections, (
        arg("reflection_list")));
  }

}}} // namespace = dials::algorithms::boost_python
