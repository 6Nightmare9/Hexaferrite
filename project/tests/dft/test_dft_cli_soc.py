import json

import pytest
import yaml

from hexw import cli as top_cli
from hexw.dft import cli as dcli
from hexw.dft import runner as R
from hexw.dft import soc


@pytest.fixture
def cfg_file(cfg, tmp_path):
    raw = json.loads(json.dumps(cfg.raw))
    raw["runs_root"] = str(tmp_path / "runs")
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump(raw), encoding="utf-8")
    return str(p)


def test_dft_is_a_real_command_in_the_main_cli(capsys, cfg_file):
    rc = top_cli.main(["dft", "--config", cfg_file, "plan", "--system", "hematite"])
    out = capsys.readouterr()
    assert rc == 0 and "J classes: 3" in out.out and "этап не реализован" not in out.err
    assert "H = -sum_{i<j} J_ij S_i.S_j" in out.out and "[pilot]" in out.out


def test_plan_json_and_generate_validate(capsys, cfg_file, tmp_path):
    assert dcli.main(["--config", cfg_file, "plan", "--system", "hematite", "--json"]) == 0
    plan = json.loads(capsys.readouterr().out)
    assert plan["n_configs"] == 5 and plan["conditioning_fit_set"]["rank"] == 4 and plan["runtime"]["n_jobs"] == 5
    assert plan["settings"]["hubbard"]["status"]["Co"] == "assumption"
    assert dcli.main(["--config", cfg_file, "generate", "--system", "hematite", "--soc"]) == 0
    out = capsys.readouterr().out
    assert "inputs validated" in out and "not_feasible_on_laptop" in out
    d = tmp_path / "runs" / "campaign_v2" / "hematite_alpha-Fe2O3"
    assert (d / "c00_ground" / "pw.in").is_file() and (d / "soc").is_dir()
    assert dcli.main(["--config", cfg_file, "status", "--system", "hematite"]) == 0
    assert "pending" in capsys.readouterr().out


def test_run_requires_the_runtime_flag(capsys, cfg_file, tmp_path, monkeypatch):
    raw = yaml.safe_load(open(cfg_file, encoding="utf-8"))
    raw["runtime"]["confirm_above_cpu_hours"] = 0.001                # make the threshold small: the estimate exceeds it
    open(cfg_file, "w", encoding="utf-8").write(yaml.safe_dump(raw))
    assert dcli.main(["--config", cfg_file, "generate", "--system", "hematite"]) == 0
    capsys.readouterr()
    called = []
    monkeypatch.setattr(R, "run_queue", lambda *a, **k: called.append(1) or {})
    rc = dcli.main(["--config", cfg_file, "run", "--system", "hematite"])
    err = capsys.readouterr()
    assert rc == dcli.EXIT_CONFIRM == 2 and "--i-accept-runtime" in err.err and not called and "CPU-h" in err.out
    assert dcli.main(["--config", cfg_file, "run", "--system", "hematite", "--i-accept-runtime"]) == 0
    assert called == [1]


def test_run_before_generate_is_an_error(capsys, cfg_file):
    assert dcli.main(["--config", cfg_file, "run", "--system", "hematite"]) == 1
    assert "generate" in capsys.readouterr().err


def test_fit_before_results_is_an_error(capsys, cfg_file):
    assert dcli.main(["--config", cfg_file, "generate", "--system", "hematite"]) == 0
    capsys.readouterr()
    assert dcli.main(["--config", cfg_file, "fit-j", "--system", "hematite"]) == 1
    assert "usable fit configurations" in capsys.readouterr().err


# ---- SOC --------------------------------------------------------------------------------------------------------------------
def test_soc_inputs_are_generated_but_never_runnable(cfg, tmp_path):
    from hexw.dft.hematite import hematite_system
    s = hematite_system()
    res = soc.generate_soc_inputs(s, cfg, tmp_path)
    assert res["state"] == "not_feasible_on_laptop" and "FEASIBILITY" in res["reference"] and len(res["jobs"]) == 6
    assert res["pseudo_fully_relativistic"] is False
    txt = (tmp_path / "soc" / "soc_t90_p0" / "pw.in").read_text()
    assert "noncolin=.true." in txt and "lspinorb=.true." in txt and "angle1(1)=90" in txt and "NEVER RUN" in txt
    assert "angle1(2)=90" in txt and "angle2(2)=180" in txt                        # antiparallel sublattice
    for j in res["jobs"]:
        assert R.read_status(__import__("pathlib").Path(j))["state"] == "not_feasible_on_laptop"
    with pytest.raises(soc.SocNotFeasible, match="not_feasible_on_laptop"):
        soc.run_soc(s, cfg)
    assert soc.soc_status(cfg)["state"] == "not_feasible_on_laptop" and cfg.raw["soc"]["allow_run"] is False


def test_soc_jobs_are_skipped_by_the_queue(cfg, tmp_path, monkeypatch):
    from hexw.dft.hematite import hematite_system
    res = soc.generate_soc_inputs(hematite_system(), cfg, tmp_path)
    called = []
    monkeypatch.setattr(R, "run_external", lambda *a, **k: called.append(1))
    from pathlib import Path
    out = R.run_queue(cfg, [Path(j) for j in res["jobs"]], timeout_s=5, force=True, retry_failed=True)
    assert set(out.values()) == {"skip:not_feasible_on_laptop"} and not called


def test_kpoints_tag_and_start_overrides(capsys, cfg_file, tmp_path):
    """k test options: a variant directory per tag, the k mesh and the Fe2+ start mode follow the options (new config hash)."""
    assert dcli.main(["--config", cfg_file, "generate", "--system", "hematite", "--kpoints", "1,1,2", "--tag", "k112",
                      "--fe2-start", "none"]) == 0
    d = tmp_path / "runs" / "campaign_v2" / "hematite_alpha-Fe2O3_k112"
    txt = (d / "c00_ground" / "pw.in").read_text()
    assert "K_POINTS automatic\n1 1 2 0 0 0" in txt
    assert dcli.main(["--config", cfg_file, "status", "--system", "hematite", "--kpoints", "1,1,2", "--tag", "k112"]) == 0
    assert "pending" in capsys.readouterr().out
    base = json.loads((d / "design.json").read_text())["config_sha256"]
    assert dcli.main(["--config", cfg_file, "generate", "--system", "hematite", "--tag", "plain"]) == 0
    other = json.loads((tmp_path / "runs" / "campaign_v2" / "hematite_alpha-Fe2O3_plain" / "design.json").read_text())["config_sha256"]
    assert base != other
    assert dcli.main(["--config", cfg_file, "plan", "--system", "hematite", "--kpoints", "1,1"]) == 1          # needs three integers
