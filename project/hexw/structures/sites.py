"""The seven transition-metal sites of the W-phase (P6_3/mmc) in the unified W1..W7 notation (literature/fit/site_notation.md).

Counts are per formula unit (Z = 2 in the hexagonal cell). Spin direction: Gorter model, confirmed by the signs of the
refined moments in L01 (Fig. 4), the mCIF of L02 and the collinear ferrimagnet of L05 (4e, 4f1, 4f3 antiparallel to
6g, 4f2, 12k, 2d in Islam's labels = W1, W2, W6 down; W3, W4, W5, W7 up).
"""
from __future__ import annotations

import numpy as np

from hexw.registry.schemas import SITE_IDS, SITE_NOTATION

SITES: tuple[str, ...] = tuple(SITE_IDS)
assert SITES == ("W1", "W2", "W3", "W4", "W5", "W6", "W7")

N_FU = np.array([SITE_NOTATION[s]["atoms_per_fu"] for s in SITES], dtype=float)       # 2 2 3 2 6 2 1
N_CELL = 2 * N_FU                                                                      # hexagonal cell, Z = 2
SPIN_SIGN = np.array([+1.0 if SITE_NOTATION[s]["gorter_spin"] == "up" else -1.0 for s in SITES])
WYCKOFF = tuple(SITE_NOTATION[s]["wyckoff"] for s in SITES)
ENV = tuple(SITE_NOTATION[s]["env"] for s in SITES)
assert N_FU.sum() == 18.0

# spin-only ionic moments of Fe (g S muB, g = 2): used by the Gorter model (hexw.params.moments) and by the weights of the Fe2+
# placement scenarios (hexw.structures.ensemble)
SPIN_FE3 = 5.0
SPIN_FE2 = 4.0
