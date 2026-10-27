"""map_literature_constants rejects NaN / inf (stage 4B)."""
import math

import pytest

from hexw.theory.magnet import map_literature_constants


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf"), math.nan])
def test_non_finite_constant_rejected(bad):
    with pytest.raises(ValueError, match="not finite"):
        map_literature_constants({"K1": bad}, {"K1": "K1"})
    with pytest.raises(ValueError, match="not finite"):
        map_literature_constants({"K1": -1e5, "K3": bad}, {"K1": "K1", "K3": "K4"})


def test_finite_still_ok():
    assert map_literature_constants({"K1": -1e5, "K3": 8e3}, {"K1": "K1", "K3": "K4"}) == {"K1": -1e5, "K4": 8e3}
