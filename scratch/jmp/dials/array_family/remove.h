#ifndef DIALS_ARRAY_FAMILY_REMOVE_H
#define DIALS_ARRAY_FAMILY_REMOVE_H

#include <scitbx/array_family/flex_types.h>
#include "../error.h"

namespace dials { namespace af {

/**
 * Remove the elements from the list and return if the element is true
 * @param in The input array
 * @param yesno True/False remove the element
 * @returns The array with the selected elements removed
 */
template <typename ArrayType>
ArrayType remove_if(const ArrayType &in, const scitbx::af::flex_bool &yesno) {
    DIALS_ASSERT(in.size() == yesno.size());
    int size = 0;
    for (int i = 0; i < yesno.size(); ++i) {
        if (!yesno[i]) size++;
    }
    ArrayType result(size);
    for (int i = 0, j = 0; i < in.size(); ++i) {
        if (!yesno[i]) {
            result[j++] = in[i];
        }
    }
    return result;
}

/**
 * Remove the elements from the list and return if the element is false
 * @param in The input array
 * @param yesno True/False don't remove the element
 * @returns The array with the selected elements removed
 */
template <typename ArrayType>
ArrayType remove_if_not(const ArrayType &in, const scitbx::af::flex_bool &yesno) {
    DIALS_ASSERT(in.size() == yesno.size());
    int size = 0;
    for (int i = 0; i < yesno.size(); ++i) {
        if (yesno[i]) size++;
    }
    ArrayType result(size);
    for (int i = 0, j = 0; i < in.size(); ++i) {
        if (yesno[i]) {
            result[j++] = in[i];
        }
    }
    return result;
}

}} // namespace dials::af

#endif // DIALS_ARRAY_FAMILY_REMOVE_H
