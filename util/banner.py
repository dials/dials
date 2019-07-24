from __future__ import absolute_import, division, print_function

import dials.util.log
import warnings

dials.util.log.print_banner()
warnings.warn(
    "dials.util.banner is deprecated. use dials.util.log.print_banner()",
    DeprecationWarning,
    stacklevel=2,
)
