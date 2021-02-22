# LIBTBX_SET_DISPATCHER_NAME libtbx.precommit

import dials.precommitbx.installer
import dials.util


@dials.util.show_mail_handle_errors()
def run(_=None):
    dials.precommitbx.installer.main()


if __name__ == "__main__":
    run()
