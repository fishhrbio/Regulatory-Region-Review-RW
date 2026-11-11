"""Rebuild the audited tables, figures and standalone explorer from raw BEDs."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import platform
from pathlib import Path
from urllib.request import urlopen
from . import AUTHOR, __version__
from .core import read_track, intersect_tracks, summarize

BASE_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/"
NAMES = {"GM12878": "wgEncodeBroadHmmGm12878HMM.bed.gz", "K562": "wgEncodeBroadHmmK562HMM.bed.gz"}
STATE_NAMES = ["Active promoter", "Weak promoter", "Inactive/poised promoter", "Strong enhancer 4",
               "Strong enhancer 5", "Weak/poised enhancer 6", "Weak/poised enhancer 7", "Insulator",
               "Transcriptional transition", "Transcriptional elongation", "Weak transcribed",
               "Polycomb-repressed", "Heterochromatin / low signal", "Repetitive/CNV 14", "Repetitive/CNV 15"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ensure_inputs(root, download):
    manifest = json.loads((root/"INPUT_MANIFEST.json").read_text())
    for entry in manifest["files"]:
        path = root / entry["path"]
        if not path.exists() and download:
            path.parent.mkdir(parents=True, exist_ok=True)
            # Do not install an unverified download as a raw source.
            payload = urlopen(BASE_URL+path.name, timeout=120).read()
            if hashlib.sha256(payload).hexdigest() != entry["sha256"]:
                raise ValueError(f"Downloaded checksum mismatch: {path.name}")
            path.write_bytes(payload)
        if not path.is_file() or path.stat().st_size != entry["bytes"] or digest(path) != entry["sha256"]:
            raise ValueError(f"Missing or changed source: {path}. Use --download for absent inputs.")
    return manifest


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(root, config_path, output, download):
    root, output = Path(root), Path(output)
    config = json.loads(Path(config_path).read_text())
    if config["assembly"] != "hg19" or config["chromosome"] != "chr22":
        raise ValueError("This release is validated for hg19 chr22 only")
    for states in [config["primary_excluded_states"], *config["classes"].values()]:
        if not states or len(states) != len(set(states)) or any(type(x) is not int or not 1 <= x <= 15 for x in states):
            raise ValueError("State definitions must be nonempty unique integer subsets of 1..15")
    # Public chart labels and narrative encode these biological definitions.
    # Fail rather than produce mislabeled exports after an unsupported edit.
    supported_classes = {"active_promoter":[1], "strong_enhancer":[4,5], "broad_enhancer":[4,5,6,7]}
    if config["classes"] != supported_classes or sorted(config["primary_excluded_states"]) != [14,15]:
        raise ValueError("This release supports the supplied classes and 14/15 mask. Update labels, narrative and validation before changing definitions.")
    manifest = ensure_inputs(root, download)
    tracks, qc = {}, {}
    for cell, filename in NAMES.items():
        tracks[cell], qc[cell] = read_track(root/"data/raw"/filename, config["chromosome"])
    atoms = intersect_tracks(tracks["GM12878"],tracks["K562"])
    common_bp = sum(a.bp for a in atoms)
    for cell in NAMES:
        qc[cell]["selected_bp_outside_common"] = qc[cell]["selected_bp"]-common_bp
    result = {"project":"Waite Chromatin Atlas", "version": __version__, "credit": AUTHOR,
              "config":config, "input_manifest":manifest, "qc":qc,
              "common_covered_bp":common_bp, "common_atom_count":len(atoms),
              "state_names":STATE_NAMES,
              "scenarios": {"primary":summarize(atoms,config["classes"],config["primary_excluded_states"]),
                            "retain_repetitive":summarize(atoms,config["classes"],[])}}
    processed = output/"data/processed"
    processed.mkdir(parents=True,exist_ok=True)
    (processed/"results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    write_csv(processed/"qc.csv",[{"cell":c, **{k:v for k,v in q.items() if k!="observed_states"}} for c,q in qc.items()])
    coverage, overlaps, matrix_rows = [],[],[]
    for scenario, data in result["scenarios"].items():
        coverage += [{"scenario":scenario, **r,"eligible_bp":data["eligible_bp"]} for r in data["coverage"]]
        overlaps += [{"scenario":scenario,"class":name,**v} for name,v in data["overlap"].items()]
        matrix_rows += [{"scenario":scenario,"gm_state":i+1,"k_state":j+1,"bp":v}
                        for i,row in enumerate(data["state_matrix_bp"]) for j,v in enumerate(row)]
    write_csv(processed/"state_coverage.csv",coverage)
    write_csv(processed/"class_overlap.csv",overlaps)
    write_csv(processed/"paired_state_matrix.csv",matrix_rows)
    # Fixed gzip timestamp and no filename in header make processed data reproducible.
    with (processed/"paired_intervals.tsv.gz").open("wb") as raw:
        with gzip.GzipFile(fileobj=raw,mode="wb",mtime=0,filename="") as gz:
            with io.TextIOWrapper(gz,encoding="utf-8",newline="\n") as text:
                text.write("chrom\tstart\tend\tgm_state\tk_state\tprimary_eligible\n")
                excluded = set(config["primary_excluded_states"])
                for a in atoms:
                    text.write(f"chr22\t{a.start}\t{a.end}\t{a.gm}\t{a.k}\t{int(a.gm not in excluded and a.k not in excluded)}\n")
    from .plots import make_figures
    make_figures(result,output/"figures")
    web_data = {**result,"atoms":[[a.start,a.end,a.gm,a.k] for a in atoms]}
    template = (root/"web/template.html").read_text(encoding="utf-8")
    app = (root/"web/app.js").read_text(encoding="utf-8")
    page = template.replace("__DATA__",json.dumps(web_data,separators=(",",":"))).replace("__APP__",app)
    (output/"docs").mkdir(parents=True,exist_ok=True)
    (output/"docs/index.html").write_text(page,encoding="utf-8")
    (output/"docs/.nojekyll").touch()
    import matplotlib
    import numpy
    environment = {"python":platform.python_version(),"platform":platform.platform(),
                   "matplotlib":matplotlib.__version__,"numpy":numpy.__version__,
                   "config_sha256":digest(Path(config_path)),"pipeline_version":__version__}
    (processed/"environment.json").write_text(json.dumps(environment,indent=2)+"\n")
    print(json.dumps({"eligible_bp":result["scenarios"]["primary"]["eligible_bp"],
                      "overlap":result["scenarios"]["primary"]["overlap"]},indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description="Rebuild Waite Chromatin Atlas from verified hg19 inputs")
    parser.add_argument("--root",type=Path,default=Path("."))
    parser.add_argument("--config",type=Path)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--download",action="store_true",help="Download absent sources with checksum validation")
    args = parser.parse_args()
    run(args.root,args.config or args.root/"config.json",args.output or args.root,args.download)
