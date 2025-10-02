from __future__ import annotations


def version():
    import os
    import sys

    import dials
    from dials.util.version import dials_version

    print(dials_version())
    print(
        f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    )
    print(f"Installed in: {os.path.split(dials.__file__)[0]}")


def run(args=None):
    import dials.util

    with dials.util.show_mail_handle_errors():
        version()


if __name__ == "__main__":
    run()
