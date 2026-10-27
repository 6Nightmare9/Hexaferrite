import pytest

from _helpers import run_as
from hexw.io import paths

READ = "from hexw.io import paths\ndef go(p):\n    return paths.open_data(p, 'rb', purpose='metadata').read()\n"
READ_NOPURPOSE = "from hexw.io import paths\ndef go(p):\n    return paths.open_data(p, 'rb').read()\n"
RAW_OPEN = "def go(p):\n    return open(p, 'rb').read()\n"


def vault_file(root):
    return root / "vault" / "lab_validation" / "d.spc"


@pytest.mark.parametrize("mod", ["hexw.theory.x", "hexw.theory.sub.y", "hexw.params.z"])
def test_theory_and_params_cannot_open_vault(root, mod):
    s = run_as(mod, READ)
    with pytest.raises(paths.TheoryAccessError):
        s["go"](vault_file(root))


@pytest.mark.parametrize("mod", ["hexw.theory.x", "hexw.params.z"])
def test_builtin_open_blocked_by_audit_hook(root, mod):
    s = run_as(mod, RAW_OPEN)
    with pytest.raises(paths.TheoryAccessError):
        s["go"](vault_file(root))


def test_theory_cannot_read_original_zip_with_lab_bytes(root):
    s = run_as("hexw.theory.x", RAW_OPEN)
    with pytest.raises(paths.TheoryAccessError):
        s["go"](root / "raw" / "original_zip" / "a.zip")


def test_theory_can_read_theory_side(root):
    s = run_as("hexw.theory.x", READ)
    assert s["go"](root / "raw" / "theory_side" / "lit.txt") == b"literature"


def test_dotdot_does_not_bypass(root):
    s = run_as("hexw.theory.x", READ)
    sneaky = root / "raw" / "theory_side" / ".." / ".." / "vault" / "lab_validation" / "d.spc"
    with pytest.raises(paths.TheoryAccessError):
        s["go"](sneaky)


@pytest.mark.parametrize("mod", ["hexw.experiment.x", "hexw.validate.y"])
def test_experiment_validate_need_freeze(root, mod):
    s = run_as(mod, READ_NOPURPOSE)
    with pytest.raises(paths.FreezeRequiredError):
        s["go"](vault_file(root))
    (root / "registry" / "FROZEN_MODEL.json").write_text('{"model_version": "mv-test"}')
    assert s["go"](vault_file(root)) == b"lab"


def test_theory_frame_below_experiment_is_still_denied(root):
    """Confused deputy: experiment helper called from theory code must not open the vault."""
    (root / "registry" / "FROZEN_MODEL.json").write_text("{}")
    exp = run_as("hexw.experiment.helper", READ_NOPURPOSE)
    th = run_as("hexw.theory.caller", "def go(f, p):\n    return f(p)\n")
    with pytest.raises(paths.TheoryAccessError):
        th["go"](exp["go"], vault_file(root))


def test_unrelated_module_has_no_vault_role(root):
    s = run_as("some.other.module", READ)
    with pytest.raises(paths.AccessError):
        s["go"](vault_file(root))


def test_audit_and_ingest_only_for_integrity_or_metadata(root):
    ok = run_as("hexw.ingest.x", READ)
    assert ok["go"](vault_file(root)) == b"lab"
    bad = run_as("hexw.ingest.x", READ.replace("purpose='metadata'", "purpose='analysis'"))
    with pytest.raises(paths.AccessError):
        bad["go"](vault_file(root))


def test_write_to_raw_and_vault_forbidden(root):
    for p in (vault_file(root), root / "raw" / "theory_side" / "lit.txt"):
        with pytest.raises(paths.ReadOnlyZoneError):
            paths.open_data(p, "wb")
        with pytest.raises(paths.ReadOnlyZoneError):
            open(p, "ab")
