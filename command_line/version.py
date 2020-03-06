from __future__ import absolute_import, division, print_function


def version():
    from dials.util.version import dials_version
    import dials
    import os
    import sys

    print(dials_version())
    print("Python {0.major}.{0.minor}.{0.micro}".format(sys.version_info))
    print("Installed in: %s" % os.path.split(dials.__file__)[0])


version()
