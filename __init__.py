import logging
import sys

if sys.version_info.major == 2:
    sys.exit("Python 2 is no longer supported")

import pathlib

_dials = pathlib.Path(__file__).parents[2]

exit(
    ("=" * 80)
    + """

Your dials repository is still tracking 'master',
but the main dials branch has been renamed to 'main'.

Please go into your dials repository at %s and run the following commands:
  git branch -m master main
  git fetch origin
  git branch -u origin/main main
  git pull --rebase

For more information please see https://github.com/dials/dials/issues/1546
"""
    % _dials
)


logging.getLogger("dials").addHandler(logging.NullHandler())

# Intercept easy_mp exceptions to extract stack traces before they are lost at
# the libtbx process boundary/the easy_mp API. In the case of a subprocess
# crash we print the subprocess stack trace, which will be most useful for
# debugging parallelized sections of DIALS code.
import libtbx.scheduling.stacktrace as _lss


def _stacktrace_tracer(error, trace, intercepted_call=_lss.set_last_exception):
    """Intercepts and prints ephemeral stacktraces."""
    if error and trace:
        logging.getLogger("dials").error(
            "\n\neasy_mp crash detected; subprocess trace: ----\n%s%s\n%s\n\n",
            "".join(trace),
            error,
            "-" * 46,
        )
    return intercepted_call(error, trace)


if _lss.set_last_exception.__doc__ != _stacktrace_tracer.__doc__:
    # ensure function is only redirected once
    _lss.set_last_exception = _stacktrace_tracer
