

#include <boost/detail/iterator.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <iostream>
#include <cstdlib>

using namespace boost::multi_index;

struct Data  
{  
  int key;  
  int number;  
  int value;
  Data(int key_,int number_,int value_)
    : key(key_),
      number(number_),
      value(value_){}  
};  


/* tags for accessing the corresponding indices*/  
struct key{};  
struct number{};  

typedef multi_index_container<  
  Data,  
  indexed_by<
    random_access <>,
    ordered_non_unique<  
      tag <key>, BOOST_MULTI_INDEX_MEMBER(Data, int, key)> >
  > map_type;  

//typedef multi_index_container<  
//  Data,  
//  indexed_by<
//    ordered_non_unique<  
//      tag <key>,    BOOST_MULTI_INDEX_MEMBER(Data, int, key)>,  
//    ordered_non_unique<  
//      tag <number>, BOOST_MULTI_INDEX_MEMBER(Data, int, number)> >  
//> map_type;  


int main()  
{  
  map_type data_map;  

  for (int k = 0; k < 10; ++k) {
    for (int n = 0; n < 5; ++n) {
      int v = std::rand();
      data_map.push_back(Data(k, n, v));
    }
  }
  
  int size = data_map.size();
  std::cout << "Size: " << size << std::endl;
  
  for (std::size_t i = 0; i < data_map.size(); ++i) {
    std::cout << 
      "Key: " << data_map[i].key << ", " << 
      "Number: " << data_map[i].number << ", " << 
      "Value: " << data_map[i].value << std::endl;
  }
  
  
//  for (int k = 0; k < 10; ++k) {
//    for (int n = 0; n < 5; ++n) {
//      int v = std::rand();
//      data_map.insert(Data(k, n, v));
//    }
//  }
  
  
  //data_map.get <0> ().find(10);
//  for (int k = 0; k < 10; ++k) {
//    for (int n = 0; n < 5; ++n) {
//      int count1 = data_map.get <0> ().count(k);
//      int count2 = data_map.get <1> ().count(n);
//      std::cout << "Count: " << count1 << ", " << count2 << std::endl;
//    }
//  }
  
  return 0; 
}

//struct data {
//  int key;
//  int value;
//};

//int main (int argc, char const* argv[])
//{
//  typedef multi_index_container <
//    data, 
//    indexed_by<    
//      random_access <>,  // keep insertion order
//      ordered_non_unique <member <data, int, &data::key> >
//    > 
//  > map_type;
//  
//  map_type map;
//  
//  
//  //typedef map_type::index <int>::type map_type_by_key;
//  //map_type_by_key::iterator it = 
//  
//  data a = { 10, 20 };
//  
//  map.insert(a);  
//  return 0;
//}
