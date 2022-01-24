"""Shared functions for tests."""
from __future__ import annotations

from unittest.mock import Mock

from dials.array_family import flex


def mock_reflection_file_object(id_=0, identifier=True):
    """Create a mock reflection_file_object."""
    fileobj = Mock()
    r = flex.reflection_table()
    r["id"] = flex.int([-1, id_, id_])
    if identifier:
        r.experiment_identifiers()[id_] = str(id_)
    fileobj.data = r
    return fileobj


def mock_two_reflection_file_object(ids=(0, 2)):
    """Create a mock reflection_file_object with two datasets."""
    fileobj = Mock()
    r = flex.reflection_table()
    r["id"] = flex.int([-1, ids[0], ids[0], ids[1], ids[1]])
    r.experiment_identifiers()[ids[0]] = str(ids[0])
    r.experiment_identifiers()[ids[1]] = str(ids[1])
    fileobj.data = r
    return fileobj
