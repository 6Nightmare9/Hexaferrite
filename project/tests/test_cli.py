import io

import pytest

from hexw import cli
from hexw.cli import STUBS
from hexw.io import paths


@pytest.mark.parametrize("cmd", [c for c in STUBS if c != "explicit-calibrate" and c not in cli.IMPLEMENTED])
def test_stub_exits_nonzero_with_message(cmd, capsys):
    rc = cli.main([cmd])
    err = capsys.readouterr().err
    assert rc == cli.EXIT_NOT_IMPLEMENTED != 0
    assert "этап не реализован" in err


def test_explicit_calibrate_requires_flag(capsys):
    assert cli.main(["explicit-calibrate"]) == cli.EXIT_CONFIRMATION_REQUIRED
    assert "подтверждение" in capsys.readouterr().err
    assert cli.main(["explicit-calibrate", cli.CONFIRM_FLAG]) == cli.EXIT_NOT_IMPLEMENTED


def test_all_tz_commands_present():
    need = {"literature", "structures", "dft", "spin-model", "predict", "freeze-model", "analyze-exp",
            "validate", "explicit-calibrate", "report", "resume"}
    assert need <= set(STUBS)


def test_audit_passes_on_real_data():
    """Integration: real audit on the populated project root (read-only, no writes)."""
    root = paths._default_root()
    if not (root / "registry" / "zip_entries.json").exists():
        pytest.skip("real project root not populated")
    old = paths.get_config().root
    try:
        from hexw.audit import run_audit
        buf = io.StringIO()
        assert run_audit(root, write_manifest=False, out=buf) == 0, buf.getvalue()
    finally:
        paths.configure(old)
