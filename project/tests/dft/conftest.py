"""Tests of hexw.dft. Tests marked `dft` start the real pw.x (minutes of CPU): they run only with HEXW_RUN_DFT=1
(`HEXW_RUN_DFT=1 python -m pytest tests/dft -m dft`), never in the ordinary test run."""
import json
import os
from pathlib import Path

import pytest

from hexw.dft import settings as ST


def pytest_configure(config):
    config.addinivalue_line("markers", "dft: runs the real Quantum ESPRESSO pw.x (slow); needs HEXW_RUN_DFT=1")


def pytest_collection_modifyitems(config, items):
    if os.environ.get("HEXW_RUN_DFT") == "1":
        return
    skip = pytest.mark.skip(reason="real pw.x run: set HEXW_RUN_DFT=1 and select with -m dft")
    for it in items:
        if it.get_closest_marker("dft") is not None:
            it.add_marker(skip)


@pytest.fixture(scope="session")
def cfg():
    return ST.load_config()


@pytest.fixture
def tmp_cfg(cfg, tmp_path):
    """Config whose runs directory is a temporary directory (absolute runs_root)."""
    return cfg.with_overrides(runs_root=str(tmp_path / "runs"))
