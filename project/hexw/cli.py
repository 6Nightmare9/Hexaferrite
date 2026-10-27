"""CLI for material forecasts, frozen comparison, independent analysis, DFT adapters and local UI."""
from __future__ import annotations

import argparse
import sys

STUBS = {
    "serve": "локальный интерфейс",
    "literature": "построение литературного корпуса (этап 2)",
    "structures": "генерация структурных кандидатов (этап 4)",
    "dft": "DFT-параметризация (этап 5)",
    "spin-model": "спиновая модель / ФМР (этап 4-6)",
    "predict": "предварительные прогнозы для серии (этап 6)",
    "freeze-model": "фиксация модели и прогнозов (этап 7)",
    "analyze-exp": "независимый анализ экспериментов (этап 8)",
    "validate": "сравнение теории и эксперимента (этап 9)",
    "explicit-calibrate": "калибровка по лабораторным данным (этап 10)",
    "report": "итоговый отчёт",
    "resume": "возобновление с контрольной точки",
}

IMPLEMENTED = {"dft", "predict", "resume", "report", "freeze-model", "analyze-exp", "validate", "spin-model", "serve"}
EXIT_NOT_IMPLEMENTED = 3
EXIT_CONFIRMATION_REQUIRED = 2
CONFIRM_FLAG = "--i-confirm-explicit-calibration"


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="hexw", description="HexW: W-hexaferrite literature model and independent analysis")
    sub = ap.add_subparsers(dest="command", required=True)
    a = sub.add_parser("audit", help="integrity + leak-guard checks and summary")
    a.add_argument("--root", help="project root; must equal the configured root (set HEXW_ROOT to change it)")
    a.add_argument("--write-manifest", action="store_true",
                   help="write a RunManifest (checksums, git HEADs) to registry/runs; by default audit writes nothing")
    for name, desc in STUBS.items():
        if name in IMPLEMENTED:
            p = sub.add_parser(name, help=desc, add_help=False)       # real parser: hexw.dft.cli
            p.add_argument("args", nargs=argparse.REMAINDER)
            continue
        p = sub.add_parser(name, help=f"[не реализовано] {desc}")
        p.add_argument("args", nargs="*")
        if name == "explicit-calibrate":
            p.add_argument(CONFIRM_FLAG, action="store_true", dest="confirm",
                           help="required: explicit confirmation that lab data may be used for calibration")
    return ap


def main(argv=None) -> int:
    args = sys.argv[1:] if argv is None else list(argv)
    if args and args[0] == "dft":                      # its own parser (hexw.dft.cli); options such as --config follow the name
        from hexw.dft.cli import main as dft_main
        return dft_main(args[1:])
    if args and args[0] in {"predict", "resume", "report"}:
        return forecast_main(args)
    if args and args[0] in {"freeze-model", "analyze-exp", "validate"}:
        return analysis_main(args)
    if args and args[0]=="spin-model":
        return spin_main(args[1:])
    if args and args[0]=="serve":
        from hexw.ui.server import main as serve_main
        return serve_main(args[1:])
    ap = build_parser()
    ns = ap.parse_args(args)
    if ns.command == "audit":
        from hexw.audit import run_audit
        try:
            return run_audit(ns.root, write_manifest=ns.write_manifest)
        except PermissionError as e:
            print(f"hexw audit: {e}", file=sys.stderr)
            return EXIT_CONFIRMATION_REQUIRED
    if ns.command == "explicit-calibrate" and not ns.confirm:
        print(f"hexw explicit-calibrate: требуется явное подтверждение флагом {CONFIRM_FLAG}. "
              "Калибровка по лабораторным данным создаёт новую версию модели и лишает данные статуса "
              "слепой проверки. Ничего не выполнено.", file=sys.stderr)
        return EXIT_CONFIRMATION_REQUIRED
    print(f"hexw {ns.command}: этап не реализован ({STUBS[ns.command]}). Результатов нет.", file=sys.stderr)
    return EXIT_NOT_IMPLEMENTED


def forecast_main(args):
    import json
    from pathlib import Path
    from hexw.workflow import friendly_resources
    friendly_resources()
    ap = argparse.ArgumentParser(prog="hexw " + args[0])
    if args[0] == "predict":
        ap.add_argument("--config", required=True)
        ap.add_argument("--out")
    else:
        ap.add_argument("--run", required=True)
    ns = ap.parse_args(args[1:])
    if args[0] == "report":
        print((Path(ns.run) / "summary.json").read_text(encoding="utf-8"))
        return 0
    from hexw.theory.request import ForecastRequest
    from hexw.theory.pipeline import run_forecast
    try:
        if args[0] == "resume":
            request = ForecastRequest(**json.loads((Path(ns.run) / "request.json").read_text(encoding="utf-8")))
            out = ns.run
        else:
            request = ForecastRequest.load(ns.config)
            out = ns.out
        result = run_forecast(request, out, progress=lambda s: print(s, flush=True))
        print(str(result), flush=True)
        return 0
    except (ValueError, OSError, RuntimeError) as e:
        print(f"forecast: {e}", file=sys.stderr)
        return 1


def analysis_main(args):
    import json
    from hexw.workflow import friendly_resources
    friendly_resources()
    ap=argparse.ArgumentParser(prog="hexw "+args[0])
    if args[0] in {"freeze-model", "validate"}:
        ap.add_argument("--run",required=True)
    if args[0]=="validate":
        ap.add_argument("--targets",required=True)
    if args[0]=="analyze-exp":
        ap.add_argument("--config",required=True)
    ns=ap.parse_args(args[1:])
    try:
        if args[0]=="freeze-model":
            from hexw.registry.freeze import freeze_run
            print(json.dumps(freeze_run(ns.run),ensure_ascii=False,indent=2))
        elif args[0]=="analyze-exp":
            from hexw.experiment.service import analyze
            print(analyze(ns.config))
        else:
            from hexw.validate.service import compare
            print(compare(ns.run,ns.targets))
        return 0
    except (ValueError,RuntimeError,OSError,KeyError) as e:
        print(f"{args[0]}: {e}",file=sys.stderr)
        return 1


def spin_main(args):
    import json
    from pathlib import Path
    import numpy as np
    from hexw.io import paths
    from hexw.workflow import friendly_resources,atomic_json,digest,utc
    friendly_resources()
    ap=argparse.ArgumentParser(prog="hexw spin-model")
    ap.add_argument("--config",required=True)
    ns=ap.parse_args(args)
    from hexw.spin.multisublattice import SpinSystem,equilibrium,response
    try:
        cfg=json.loads(Path(ns.config).read_text(encoding="utf-8"))
        system=SpinSystem(**cfg["system"])
        if "unit_magnetizations" in cfg:
            m=np.asarray(cfg["unit_magnetizations"]); info={"equilibrium":"user-specified; stability checked"}
        else:
            m,info=equilibrium(system,cfg["B_T"],starts=cfg.get("starts",12))
        result=response(system,m,cfg["B_T"],cfg["frequencies_Hz"])
        out=paths.get_config().root/"runs"/"spin"/digest(cfg)[:20]
        out.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(out/"response.npz",**{k:v for k,v in result.items() if isinstance(v,np.ndarray)})
        atomic_json(out/"result.json",{"config":cfg,"equilibrium":info,"unit_magnetizations":m,
            "undamped_modes_Hz":result["undamped_modes_Hz"],"status":"explicit-parameter scenario; not a formula-derived exchange prediction","created_utc":utc()})
        print(out); return 0
    except (ValueError,RuntimeError,OSError,KeyError) as e:
        print(f"spin-model: {e}",file=sys.stderr); return 1
