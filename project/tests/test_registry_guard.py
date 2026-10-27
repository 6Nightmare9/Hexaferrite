"""Protection of the registry against being overwritten by build_registry / umnik (synthetic root only)."""
import os
import subprocess

import pytest
import yaml

from hexw.ingest import build_registry, umnik
from hexw.ingest import registry_guard as G


def _git(reg, *a):
    subprocess.run(["git", "-C", str(reg), "-c", "user.email=t@t", "-c", "user.name=t", *a], check=True,
                   capture_output=True)


def _init_clean_repo(reg):
    _git(reg, "init", "-q")
    _git(reg, "add", "-A")
    _git(reg, "commit", "-q", "-m", "x")


def test_build_registry_refuses_existing_files(root):
    reg = root / "registry"
    (reg / "samples.yaml").write_text("samples: [MANUAL]\n")
    with pytest.raises(G.RegistryOverwriteError, match="already exist"):
        build_registry.build()
    assert (reg / "samples.yaml").read_text() == "samples: [MANUAL]\n"


def test_umnik_refuses_existing_files(root):
    reg = root / "registry"
    (reg / "prior_use.yaml").write_text("manual: 1\n")
    with pytest.raises(G.RegistryOverwriteError):
        umnik.build()
    assert (reg / "prior_use.yaml").read_text() == "manual: 1\n"


def test_force_needs_git_repo_and_clean_status(root):
    reg = root / "registry"
    (reg / "samples.yaml").write_text("s: 1\n")
    with pytest.raises(G.RegistryOverwriteError, match="not a git repository"):
        G.check_write_allowed(reg, ["samples.yaml"], force=True)
    _init_clean_repo(reg)
    assert G.check_write_allowed(reg, ["samples.yaml"], force=True) == ["samples.yaml"]
    (reg / "samples.yaml").write_text("s: 2   # manual edit\n")          # uncommitted manual edit
    with pytest.raises(G.RegistryOverwriteError, match="not clean"):
        G.check_write_allowed(reg, ["samples.yaml"], force=True)
    with pytest.raises(G.RegistryOverwriteError):                        # without force always refused
        G.check_write_allowed(reg, ["samples.yaml"], force=False)


def test_no_existing_files_needs_no_force(root):
    assert G.check_write_allowed(root / "registry", ["samples.yaml"]) == []


def test_cli_exit_code_on_refusal(root, capsys):
    (root / "registry" / "samples.yaml").write_text("s: 1\n")
    assert build_registry.main([]) == 2
    assert "Refusing" in capsys.readouterr().out


EXISTING = """# header kept
schema_version: 1
independent_lab_validation_available: false
independent_lab_validation_note: set by hand
entries:
- observation_id: OBS-1
  prior_use_status: used_prior
  validation_mode: retrospective
  user_confirmation: {by: user, date: 2026-10-25}
"""


def _gen(oid, status="no_prior_use_found"):
    return {"observation_id": oid, "prior_use_status": status, "validation_mode": "independent_candidate",
            "user_confirmation": None}


def test_prior_use_merge_keeps_manual_fields_byte_for_byte():
    text, rep = G.merge_prior_use_text(EXISTING, [_gen("OBS-1"), _gen("OBS-2")])
    assert text.startswith(EXISTING)                                     # nothing manual changed
    doc = yaml.safe_load(text)
    assert doc["independent_lab_validation_available"] is False
    assert doc["entries"][0]["user_confirmation"]["by"] == "user"
    assert [e["observation_id"] for e in doc["entries"]] == ["OBS-1", "OBS-2"]
    assert rep["added"] == ["OBS-2"]
    assert rep["status_differs_from_generated_kept_manual"] == ["OBS-1"]  # the manual status wins, reported


def test_prior_use_merge_nothing_new_is_noop():
    text, rep = G.merge_prior_use_text(EXISTING, [_gen("OBS-1")])
    assert text == EXISTING and rep["added"] == []


def test_prior_use_merge_refuses_when_entries_not_last():
    bad = "entries:\n- observation_id: OBS-1\nextra: 1\n"
    with pytest.raises(G.RegistryOverwriteError):
        G.merge_prior_use_text(bad, [_gen("OBS-2")])


def test_write_keeps_read_only_attribute(tmp_path):
    p = tmp_path / "f.yaml"
    p.write_text("a")
    os.chmod(p, 0o444)
    G.write_text_keep_mode(p, "b")
    assert p.read_text() == "b" and not os.access(p, os.W_OK)
    os.chmod(p, 0o666)


# ---- stage 4A: style of `entries:` (flow-style lists used to turn into invalid YAML) ----
FLOW_EMPTY = "# header\nschema_version: 1\nentries: []   # nothing yet\n"


@pytest.mark.parametrize("nl", ["\n", "\r\n"])
def test_prior_use_merge_into_empty_flow_list(nl):
    existing = FLOW_EMPTY.replace("\n", nl)
    text, rep = G.merge_prior_use_text(existing, [_gen("OBS-1"), _gen("OBS-2")])
    doc = yaml.safe_load(text)
    assert [e["observation_id"] for e in doc["entries"]] == ["OBS-1", "OBS-2"]
    assert doc["schema_version"] == 1 and rep["added"] == ["OBS-1", "OBS-2"]
    assert "# nothing yet" in text and "[]" not in text                    # the hand comment survives
    assert text.startswith("# header" + nl + "schema_version: 1" + nl + "entries:")
    assert ("\r\n" in text) == (nl == "\r\n") and (nl == "\r\n" or "\r" not in text)


def test_prior_use_merge_refuses_non_empty_flow_list():
    flow = "schema_version: 1\nentries: [{observation_id: OBS-1, prior_use_status: used_prior}]\n"
    with pytest.raises(G.RegistryOverwriteError, match="flow"):
        G.merge_prior_use_text(flow, [_gen("OBS-1"), _gen("OBS-2")])
    multi = "schema_version: 1\nentries: [\n  {observation_id: OBS-1}\n]\n"
    with pytest.raises(G.RegistryOverwriteError, match="flow"):
        G.merge_prior_use_text(multi, [_gen("OBS-2")])


def test_prior_use_merge_into_indented_block_list():
    indented = "schema_version: 1\nentries:\n  - observation_id: OBS-1\n    prior_use_status: used_prior\n"
    text, rep = G.merge_prior_use_text(indented, [_gen("OBS-1"), _gen("OBS-2")])
    assert text.startswith(indented)
    assert [e["observation_id"] for e in yaml.safe_load(text)["entries"]] == ["OBS-1", "OBS-2"]


def test_prior_use_merge_verifies_the_result(monkeypatch):
    # a style the detector does not understand must not produce a silently corrupt file
    weird = "schema_version: 1\nentries: &a\n- observation_id: OBS-1\n"
    with pytest.raises(G.RegistryOverwriteError):
        G.merge_prior_use_text(weird, [_gen("OBS-2")])


def test_write_prior_use_checked_rolls_back(tmp_path):
    p = tmp_path / "prior_use.yaml"
    p.write_bytes(EXISTING.encode())
    good, _ = G.merge_prior_use_text(EXISTING, [_gen("OBS-2")])
    G.write_prior_use_checked(p, good)
    assert yaml.safe_load(p.read_text())["entries"][-1]["observation_id"] == "OBS-2"
    before = p.read_bytes()
    with pytest.raises(G.RegistryOverwriteError, match="rolled back"):
        G.write_prior_use_checked(p, before.decode() + "\n  broken: [unclosed\n")
    assert p.read_bytes() == before                                         # old bytes restored
    q = tmp_path / "new.yaml"
    with pytest.raises(G.RegistryOverwriteError):
        G.write_prior_use_checked(q, "a: [unclosed\n")
    assert not q.exists()
