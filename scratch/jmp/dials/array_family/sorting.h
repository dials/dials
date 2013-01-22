#ifndef DIALS_ARRAY_FAMILY_SORTING_H
#define DIALS_ARRAY_FAMILY_SORTING_H

#include <scitbx/array_family/flex_types.h>
#include "../error.h"

namespace dials { namespace array_family {

/**
 * A class to aid comparison based on an index array
 * @tparam IndexIterator The index array iterator type
 * @tparam DataIterator The data array iterator type
 * @tparam Compare The comparison function type
 */
template <typename IndexIterator, typename DataIterator, typename Compare>
class CompareIndex {
public:

    typedef typename std::iterator_traits <IndexIterator>::value_type IndexType;
    typedef typename std::iterator_traits <DataIterator>::value_type DataType;

    /**
     * Initialise and keep hold of the data iterator and comparison function
     * @param data The data iterator
     * @param compare The comparison function
     */      
    CompareIndex(const DataIterator data, Compare compare) 
        : data_(data), 
          compare_(compare) {}
    
    /**
     * Compare the two elements at the given indices a and b
     * @param a Index A
     * @param b Index B
     */
    bool operator()(const IndexType &a, const IndexType &b) {
        return compare_(data_[a], data_[b]);
    }

private:

    const DataIterator data_;
    Compare compare_;
};

/**
 * Perform a partial sort of the data using an index list. The data itself will
 * not be sorted but the list of indices will be sorted in data order. This
 * variant of the function uses the std::greater function object for comparisons.
 * @param index_first The first element of the index array
 * @param index_middle The middle element of the index array to sort about
 * @param index_last The last element in the index array
 * @param data_first The first element in the data array
 */ 
template <typename IndexIterator, typename DataIterator>
void partial_sort_indices(IndexIterator index_first, 
                          IndexIterator index_middle,
                          IndexIterator index_last,
                          const DataIterator data_first)
{
    typedef typename std::iterator_traits <IndexIterator>::value_type IndexType;
    typedef typename std::iterator_traits <DataIterator>::value_type DataType;
    std::partial_sort(index_first, index_middle, index_last, 
                      CompareIndex <IndexIterator, 
                                    DataIterator, 
                                    std::greater <DataType> > (
                                        data_first, std::greater <DataType>()));
}

/**
 * Perform a partial sort of the data using an index list. The data itself will
 * not be sorted but the list of indices will be sorted in data order. This
 * variant of the function uses the given compare function for comparisons.
 * @param index_first The first element of the index array
 * @param index_middle The middle element of the index array to sort about
 * @param index_last The last element in the index array
 * @param data_first The first element in the data array
 * @param compare The comparison function
 */ 
template <typename IndexIterator, typename DataIterator, typename Compare>
void partial_sort_indices(IndexIterator index_first, 
                          IndexIterator index_middle,
                          IndexIterator index_last,
                          const DataIterator data_first,
                          Compare compare)
{
    std::partial_sort(index_first, index_middle, index_last, 
                      CompareIndex <IndexIterator, 
                                    DataIterator, 
                                    Compare> (data_first, compare));
}

/**
 * Partial sort a flex array based on the array indices
 * @param data The array to sort
 * @param n The number of elements to sort
 * @returns An array of sorted indices
 */
template <typename ArrayType>
scitbx::af::flex_int partial_sort_indices_flex(const ArrayType &data, int n) {
    DIALS_ASSERT(0 < n && n <= data.size());
    scitbx::af::flex_int result(data.size());
    for (int i = 0; i < result.size(); ++i) result[i] = i;
    partial_sort_indices(result.begin(), 
                         result.begin() + n, 
                         result.end(), 
                         data.begin());
    return result;
}

}} // namespace dials::array_family

#endif // DIALS_ARRAY_FAMILY_SORTING_H
