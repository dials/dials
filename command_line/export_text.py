# LIBTBX_SET_DISPATCHER_NAME dev.dials.export_text
from __future__ import absolute_import, division, print_function

from dials.util.export_text import export_text
from dials.array_family import flex

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        raise RuntimeError("%s integrated.mpack")

    integrated_data = flex.reflection_table.from_msgpack_file(sys.argv[1])
    export_text(integrated_data)
