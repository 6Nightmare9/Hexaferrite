"""hexw: literature-based W-hexaferrite model with frozen predictions and independent analysis."""
__version__ = "0.2.0"

# The data-isolation audit hook must be active in EVERY process that imports hexw (also `python -m hexw.<x>`).
from hexw.io import paths  # noqa: E402,F401
