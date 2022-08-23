from __future__ import annotations

import sys
import warnings

from dials.command_line.split import *  # noqa: F403

try:
    import colorama

    YELLOW = colorama.Fore.YELLOW
    NC = colorama.Style.RESET_ALL
except ImportError:
    YELLOW = ""
    NC = ""
else:
    colorama.init()


if __name__ == "__main__":
    print(
        f"\n{YELLOW}Warning: dials.split_experiments is now deprecated and will be removed in the future. Please run dials.split instead.\n\n{NC}",
        file=sys.stderr,
    )
    run()  # noqa: F405
else:
    warnings.warn(
        "dials.command_line.split_experiments is deprecated. Please import dials.command_line.split instead.",
        UserWarning,
        stacklevel=1,
    )
