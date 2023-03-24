from __future__ import annotations

import boost_adaptbx.boost.python as bp
from iotbx.mtz import extract_from_symmetry_lib

from dials_util_ext import *  # noqa: F403; lgtm
from dials_util_ext import GemmiMtzObject


@bp.inject_into(GemmiMtzObject)
class _:
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


__all__ = (  # noqa: F405
    "ResolutionMaskGenerator",
    "add_dials_batches",
    "dials_u_to_mosflm",
    "ostream",
    "scale_down_array",
    "streambuf",
    "GemmiMtzObject",
)
