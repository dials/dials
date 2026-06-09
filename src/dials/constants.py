from __future__ import annotations

import math

# constants removed from elsewhere - EPS which is the "standard small number"
EPS = 1e-7

# FULL_PARTIALITY which corresponds to the greatest partiality we can
# reasonably expect integrating to +/- 3 sigma
FULL_PARTIALITY = math.erf(3 / math.sqrt(2))
