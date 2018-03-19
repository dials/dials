#!/usr/bin/python
#
# this is used by the build scripts - it should not normally need to be
# called in a full DIALS installation (although this will probably work too)

"""
Echo the current DIALS version and build number.
"""

from __future__ import absolute_import, division, print_function

import datetime

# XXX TO CHANGE THE DIALS VERSION NAME, EDIT BASE_VERSION
BASE_VERSION = "dev"
# this is when the first DIALS installer was created (UTC time)
start_date = datetime.date(2014, 7, 10)
today = datetime.datetime.utcnow().date()
delta = today - start_date
n_days = abs(delta.days) + 1
print("%s-%d" % (BASE_VERSION, n_days))
