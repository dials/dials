from __future__ import annotations


def version():
    import os
    import sys

    import dials
    from dials.util.version import dials_version

    print(dials_version())
    print("Python {0.major}.{0.minor}.{0.micro}".format(sys.version_info))
    print(f"Installed in: {os.path.split(dials.__file__)[0]}")


def run(args=None):
    import dials.util

    with dials.util.show_mail_handle_errors():
        version()


if __name__ == "__main__":
    run()
