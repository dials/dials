from __future__ import annotations

import boost_adaptbx.boost.python as bp
from iotbx.mtz import extract_from_symmetry_lib

from dials_util_ext import *  # noqa: F403; lgtm
from dials_util_ext import GemmiMtzObject


class CrystalView:
    def __init__(self, crystal_name, project_name, unit_cell_parameters, mtz_object):

        self.datasets = []
        self.crystal_name = crystal_name
        self.project_name = project_name
        self.unit_cell_parameters = unit_cell_parameters
        self.mtz_object = mtz_object

        return

    def add_dataset(self, dataset_name, wavelength):

        dataset = DatasetView(dataset_name, wavelength, self)
        self.mtz_object.add_dataset(
            self.project_name,
            self.crystal_name,
            dataset_name,
            self.unit_cell_parameters,
            wavelength,
        )
        self.datasets.append(dataset)

        return dataset


class DatasetView:
    def __init__(self, dataset_name, wavelength, crystal_view):

        self.dataset_name = dataset_name
        self.wavelength = wavelength
        self.crystal_view = crystal_view

        return

    def add_column(self, column_name, column_type):

        self.crystal_view.mtz_object.add_column(column_name, column_type)
        column = ColumnView(column_name, column_type, self)
        return column


class ColumnView:
    def __init__(self, column_name, column_type, dataset_view):

        self.column_name = column_name
        self.column_type = column_type
        self.dataset_view = dataset_view

    def set_values(self, array):

        # Do not set dummy data for these columns, as they will be set by
        # replace_original_index_miller_indices instead.
        if self.column_name not in ("H", "K", "L", "M_ISYM"):
            self.dataset_view.crystal_view.mtz_object.add_column_data(array)


@bp.inject_into(GemmiMtzObject)
class _:

    crystals = []

    def set_space_group_info(self, space_group_info, symbol=None):
        if symbol is None:
            symbol = extract_from_symmetry_lib.ccp4_symbol(
                space_group_info=space_group_info, lib_name="symop.lib"
            )
            if symbol is None:
                symbol = "No.%d" % space_group_info.type().number()
        self.set_space_group_by_name(symbol)
        # FIXME is this enough? Note that the equivalent method in iotbx.mtz.object
        # has all this extra code:
        #
        # group = space_group_info.group()
        # self.set_space_group_name(name=symbol)
        # self.set_space_group_number(number=space_group_info.type().number())
        # self.set_point_group_name(name=group.point_group_type())
        # self.set_lattice_centring_type(
        #  symbol=group.conventional_centring_type_symbol())
        # if (self.lattice_centring_type() == "\0"):
        #  self.set_lattice_centring_type(symbol="?")
        # self.set_space_group(space_group=space_group_info.group())
        # return self

    def add_crystal(self, crystal_name, project_name, unit_cell_parameters):
        """Cache information about a crystal to add with the dataset later"""

        crystal = CrystalView(crystal_name, project_name, unit_cell_parameters, self)
        self.crystals.append(crystal)
        return crystal

    def adjust_column_array_sizes(self, nref):
        pass


__all__ = (  # noqa: F405
    "ResolutionMaskGenerator",
    "add_dials_batches",
    "dials_u_to_mosflm",
    "ostream",
    "scale_down_array",
    "streambuf",
    "GemmiMtzObject",
)
