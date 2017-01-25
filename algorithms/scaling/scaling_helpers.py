#!/usr/bin/env cctbx.python
#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Auxiliary functions for the scaling package"""

from __future__ import absolute_import, division

def products_omitting_one_item(items):
  """For n items, efficiently calculate every product of n-1 items. That is
  calculate every possible product that omits a single item. For example,
  given the list of items [a, b, c, d], calculate and return the products

  1) b*c*d
  2) a*c*d
  3) a*b*d
  4) a*b*c.

  The aim here is to avoid calculating sub-products such as c*d more than once.
  """

  items = list(items)
  n = len(items)

  # The function is not defined for 1-elt lists
  assert n > 1

  # cumulatively calculate n-1 products from the front
  p = items[0]
  front_prod = [p]
  for i in xrange(1, n-1):
    p *= items[i]
    front_prod.append(p)

  # cumulatively calculate n-1 products from the back
  items.reverse()
  p = items[0]
  back_prod = [p]
  for i in xrange(1, n-1):
    p *= items[i]
    back_prod.append(p)

  back_prod.reverse()
  results = [back_prod[0]]
  for i in xrange(1, n-1):
    results.append(back_prod[i] * front_prod[i - 1])
  results.append(front_prod[n - 2])

  return results

