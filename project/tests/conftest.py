import pytest

from _helpers import run_as  # noqa: F401  (re-export for old imports)
from hexw.io import paths


@pytest.fixture
def root(tmp_path):
    (tmp_path / "raw" / "theory_side").mkdir(parents=True)
    (tmp_path / "raw" / "original_zip").mkdir(parents=True)
    (tmp_path / "vault" / "lab_validation").mkdir(parents=True)
    (tmp_path / "vault" / "legacy").mkdir(parents=True)
    (tmp_path / "registry").mkdir()
    (tmp_path / "raw" / "theory_side" / "lit.txt").write_text("literature")
    (tmp_path / "vault" / "lab_validation" / "d.spc").write_bytes(b"lab")
    (tmp_path / "raw" / "original_zip" / "a.zip").write_bytes(b"zip")
    with paths.test_root(tmp_path):      # explicit, test-only; grants no vault access by itself
        yield tmp_path
