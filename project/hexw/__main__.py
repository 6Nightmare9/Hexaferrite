import os
import sys
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(key, "2")
from hexw.cli import main

if __name__ == "__main__":
    sys.exit(main())
