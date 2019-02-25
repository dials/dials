"""
Tests for the reindex_to_reference module.
"""
from __future__ import absolute_import, division, print_function
import pytest
from dials.util import Sorry

from cctbx import sgtbx

from dials.algorithms.symmetry.cosym._generate_test_data import generate_test_data
from dials.algorithms.symmetry.reindex_to_reference import \
  determine_reindex_operator_against_reference

def test_determine_reindex_operator_against_reference():
  """Test that the correct reindex operator is returned by the function."""

  # create a test dataset of random intensities
  data, _ = generate_test_data(
    space_group=sgtbx.space_group_info(symbol='P4').group(), sample_size=1)

  # reindexing operator is a,-b,-c
  op = 'a,-b,-c'
  # reindex the data into a new array, so that the function should determine
  # that the same change of basis operator should be applied to give consistent
  # indexing back to the original data.
  reindexed_data = data[0].change_basis(op)

  cb_op = determine_reindex_operator_against_reference(data[0], reindexed_data)
  assert cb_op.as_abc() == 'a,-b,-c'

  # Repeat but with no reindexing
  cb_op = determine_reindex_operator_against_reference(data[0], data[0])
  assert cb_op.as_abc() == 'a,b,c'

  #Test that a Sorry is raised if inconsistent indexing
  data_2, _ = generate_test_data(
    space_group=sgtbx.space_group_info(symbol='P1').group(), sample_size=1)
  with pytest.raises(Sorry):
    cb_op = determine_reindex_operator_against_reference(data[0], data_2[0])

  # Test case for a simple space group with no ambiguity
  data, _ = generate_test_data(
    space_group=sgtbx.space_group_info(symbol='P1').group(), sample_size=1)
  cb_op = determine_reindex_operator_against_reference(data[0], data[0])
  assert cb_op.as_abc() == 'a,b,c'
