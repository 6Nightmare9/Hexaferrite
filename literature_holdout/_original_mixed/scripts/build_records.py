"""Assemble records.yaml / records.jsonl, validate against hexw.registry.schemas (read-only import of project/).
Run: D:/Work/Hexx/.venv/Scripts/python.exe -I literature/scripts/build_records.py   (cwd = D:/Work/Hexx)"""
import sys, json, collections, datetime
sys.dont_write_bytecode = True
sys.path.insert(0, r"D:/Work/Hexx/literature/scripts")
sys.path.insert(0, r"D:/Work/Hexx/project")           # read-only import of the schema module
import yaml
import rc
import data_morch, data_nakai_islam, data_tier3          # noqa: F401  (populate rc.REG)
from hexw.registry.schemas import ParameterEstimate, LiteratureRecord, Claim, from_dict, to_plain
from hexw.units import check_unit

OUT = r"D:/Work/Hexx/literature"
REG = rc.REG

# ---- validate estimates
bad = 0
for r in REG:
    if r["estimate"]:
        try:
            pe = from_dict(ParameterEstimate, r["estimate"])
            assert pe.value == r["estimate"]["value"]
        except Exception as e:                       # noqa
            bad += 1
            print("INVALID", r["rec_id"], r["name"], e)
print("estimates:", sum(1 for r in REG if r["estimate"]), "invalid:", bad)

# ---- literature records with claims (schema-valid)
by_lit = collections.defaultdict(list)
for r in REG:
    by_lit[r["lit_id"]].append(r)
lits = []
for lid, rs in by_lit.items():
    info = rc.LIT[lid]
    claims = []
    for r in rs:
        e = r["estimate"]
        c = dict(text=f"[{r['rec_id']}] {r['name']}: {r['ext']['verbatim']}", locator=r["ext"]["locator"], extraction_status="verified" if r["ext"]["status"] != "derived_by_us" else "derived_by_us")
        if e:
            c.update(value=e["value"], unit=e["unit"], quantity_kind=e["quantity_kind"])
        claims.append(c)
    lr = dict(lit_id=lid, citation=info["citation"], doi=info["doi"], year=info["year"], path=info.get("path"), access=info["access"], role=None, claims=claims, license=info.get("license"))
    obj = from_dict(LiteratureRecord, lr)                  # raises if a numeric claim lacks unit/kind
    lits.append(to_plain(obj))

# ---- emit
meta = dict(
    generated_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    stage="Stage 2 part 1 (literature extraction)", n_records=len(REG), n_schema_valid_estimates=sum(1 for r in REG if r["estimate"]),
    n_pending_kind=sum(1 for r in REG if r["pending"] and r["pending"]["reason"].startswith("schema")),
    n_claims_only=sum(1 for r in REG if not r["estimate"] and not r["pending"]),
    note=("literature_records: strictly hexw.registry.schemas.LiteratureRecord/Claim (loadable with from_dict). parameter_records: 'estimate' is a ParameterEstimate dict "
          "(loadable); where the schema has no suitable Kind/unit the numeric value sits in 'pending' with proposed_quantity_kind; 'ext' carries the extra fields demanded by TZ section 4/13 "
          "(sample id, synthesis, measurement conditions, definition, status, verbatim, transformation, measurement_group, split, tier, figure_only). source_role of estimates = proposed split "
          "(literature_fit / literature_holdout); records with split excluded_* / methodology_only carry no estimate. NOTHING here is final: the split is a PROPOSAL (see split_proposal.md)."),
    status_vocabulary=["measured", "calculated", "digitized", "assumed_by_authors", "derived_by_us"],
)
doc = dict(meta=meta, literature_records=lits, parameter_records=REG)
with open(OUT + "/records.yaml", "w", encoding="utf-8") as fh:
    yaml.safe_dump(doc, fh, allow_unicode=True, sort_keys=False, width=200, default_flow_style=False)
with open(OUT + "/records.jsonl", "w", encoding="utf-8") as fh:
    for r in REG:
        fh.write(json.dumps(r, ensure_ascii=False) + "\n")

# ---- statistics (for the report)
cnt_src = collections.Counter(r["lit_id"] for r in REG)
cnt_split = collections.Counter(r["ext"]["split"] for r in REG)
cnt_tier = collections.Counter(r["ext"]["tier"] for r in REG)
cnt_fig = sum(1 for r in REG if r["ext"]["figure_only"])
print("total", len(REG)); print("by source", dict(cnt_src)); print("by split", dict(cnt_split)); print("by tier", dict(cnt_tier)); print("figure_only", cnt_fig)
prop = collections.Counter(r["name"].split("_")[0] + "_" + (r["name"].split("_")[1] if "_" in r["name"] else "") for r in REG)
print("top names", prop.most_common(25))
json.dump(dict(by_source=cnt_src, by_split=cnt_split, by_tier={str(k): v for k, v in cnt_tier.items()}, figure_only=cnt_fig), open(OUT + "/scripts/_stats.json", "w"))
