# LIBTBX_SET_DISPATCHER_NAME libtbx.precommit

from __future__ import absolute_import, division, print_function

import dials.precommitbx.installer


def run(_=None):
    dials.precommitbx.installer.main()


if __name__ == "__main__":
    run()
