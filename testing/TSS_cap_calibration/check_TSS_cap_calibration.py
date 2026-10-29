#!/usr/bin/env python3
# encoding: utf-8

"""Execution test for the library-level TSS cap calibration (v0.47.0). See README.md.

Runs the standalone step (util/calibrate_TSS_cap_signature.py) on saved scan tables --
the per-peak and off-peak 5' G-run counts a bam scan produces -- and checks the
decision, the config keys it hands to LRAA, and the report: .json (decision and
statistics), .hist.tsv (per-peak capped-fraction histogram) and .svg (its plot). The
runs are kept under __test_runs/, with a summary.md showing every case's plot.

Cases (data/, see README.md for where each scan came from):
  BT474_on            bulk PacBio: the signature, baseline 1 G -> "on"
  SIRV_undetermined   SIRV spike-ins: 38 peaks, too few to judge -> "undetermined",
                      no config keys (LRAA keeps "auto")
  sim_off             Arabidopsis simulation: one mode at 0 -> "off"
  BT474_floor         BT474 with the peak floor raised past its peak count by
                      --config_update: the caller's settings reach the decision
  BT474_merge         --merge: the output config keeps the caller's keys
"""

import os
import sys
import json
import shutil
import subprocess
import xml.dom.minidom

HERE = os.path.dirname(os.path.realpath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
CALIBRATE = os.path.join(REPO, "util", "calibrate_TSS_cap_signature.py")
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "__test_runs")

failures = []


def check(cond, msg):
    print(("  ok    " if cond else "  FAIL  ") + msg)
    if not cond:
        failures.append(msg)


def calibrate(name, table, config=None, merge=False):
    d = os.path.join(OUT, name)
    os.makedirs(d)
    cmd = [sys.executable, CALIBRATE, "--scan_table", os.path.join(DATA, table), "--output_prefix", name]
    if config is not None:
        cmd += ["--config_update", json.dumps(config)]
    if merge:
        cmd += ["--merge"]
    with open(os.path.join(d, "run.log"), "wt") as log:
        rc = subprocess.call(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT)
    check(rc == 0, "exit 0")
    return os.path.join(d, name)


def check_report(prefix, decision, baseline):
    files = {ext: "{}.TSS_cap_calibration.{}".format(prefix, ext) for ext in ("json", "hist.tsv", "svg", "config.json")}
    for p in files.values():
        check(os.path.exists(p) and os.path.getsize(p) > 0, "wrote " + os.path.basename(p))
    rep = json.load(open(files["json"]))
    st, cu = rep["stats"], rep["config_update"]
    check(st["decision"] == decision, "decision {} (got {})".format(decision, st["decision"]))
    check(st["baseline_G_run"] == baseline, "baseline G run {} (got {})".format(baseline, st["baseline_G_run"]))
    check(sum(st["histogram"]) == st["peaks_judged"], "histogram covers every judged peak")
    rows = open(files["hist.tsv"]).read().strip().split("\n")[1:]
    check(sum(int(r.split("\t")[2]) for r in rows) == st["peaks_judged"], "hist.tsv covers every judged peak")
    try:
        xml.dom.minidom.parse(files["svg"])
        check(True, "svg is well-formed")
    except Exception as e:
        check(False, "svg is well-formed ({})".format(e))
    if decision in ("on", "off"):
        check(cu == {"TSS_cap_signature": decision, "TSS_cap_G_baseline": baseline},
              "config_update sets TSS_cap_signature={} and TSS_cap_G_baseline={}".format(decision, baseline))
    else:
        check(cu == {}, "undetermined: no config keys (LRAA keeps auto)")
    return rep, json.load(open(files["config.json"]))


def main():
    shutil.rmtree(OUT, ignore_errors=True)
    os.makedirs(OUT)
    summary = []

    print("BT474_on: bulk PacBio")
    p = calibrate("BT474_on", "BT474_PacBio.scan.tsv.gz")
    rep, cfg = check_report(p, "on", 1)
    check(cfg == {"TSS_cap_signature": "on", "TSS_cap_G_baseline": 1, "TSS_cap_prepass": False},
          "config.json: the calibration, and no re-calibration downstream")
    st = rep["stats"]
    check(st["capped_mode_frac"] > 0.8 and st["mode_separation"] > 0.5, "capped mode dominant and well separated")
    check(st["valley"] is not None and st["valley"] >= st["per_site_threshold"], "valley at or above the per-site threshold")
    summary.append(("BT474_on", rep))

    print("SIRV_undetermined: SIRV spike-ins")
    p = calibrate("SIRV_undetermined", "SIRV.scan.tsv.gz")
    rep, cfg = check_report(p, "undetermined", 1)
    check("too few peaks" in rep["stats"].get("reason", ""), "reason: too few peaks")
    check(cfg == {"TSS_cap_prepass": False}, "config.json: only TSS_cap_prepass false (runs keep auto per graph)")
    summary.append(("SIRV_undetermined", rep))

    print("sim_off: Arabidopsis simulation")
    p = calibrate("sim_off", "Arabidopsis_sim.Chr4_6M.scan.tsv.gz")
    rep, cfg = check_report(p, "off", 0)
    check(rep["stats"]["mode_separation"] < 0.2, "modes not separated")
    summary.append(("sim_off", rep))

    print("BT474_floor: --config_update reaches the decision")
    p = calibrate("BT474_floor", "BT474_PacBio.scan.tsv.gz", config={"TSS_cap_prepass_min_peaks": 100000})
    rep, cfg = check_report(p, "undetermined", 1)
    summary.append(("BT474_floor", rep))

    print("BT474_merge: --merge keeps the caller's keys")
    p = calibrate("BT474_merge", "BT474_PacBio.scan.tsv.gz", config={"require_cap_for_new_TSS": True}, merge=True)
    rep, cfg = check_report(p, "on", 1)
    check(cfg == {"require_cap_for_new_TSS": True, "TSS_cap_signature": "on", "TSS_cap_G_baseline": 1,
                  "TSS_cap_prepass": False}, "config.json: caller's key plus the calibration")

    with open(os.path.join(OUT, "summary.md"), "wt") as ofh:
        ofh.write("# TSS cap calibration: execution test\n\n")
        ofh.write("| case | decision | baseline G | peaks judged | uncapped mode median | capped mode median "
                  "| capped mode share | separation | valley | reason |\n")
        ofh.write("|---|---|---|---|---|---|---|---|---|---|\n")
        keys = ("decision", "baseline_G_run", "peaks_judged", "uncapped_mode_median", "capped_mode_median",
                "capped_mode_frac", "mode_separation", "valley", "reason")
        for name, rep in summary:
            s = rep["stats"]
            ofh.write("| " + " | ".join([name] + ["" if s.get(k) is None else str(s[k]) for k in keys]) + " |\n")
        ofh.write("\n")
        for name, _ in summary:
            ofh.write("## {}\n\n![]({}/{}.TSS_cap_calibration.svg)\n\n".format(name, name, name))
    print("\nreports and plots: {}".format(os.path.join(OUT, "summary.md")))

    if failures:
        print("\n{} check(s) FAILED".format(len(failures)))
        sys.exit(1)
    print("\nall checks passed")


if __name__ == "__main__":
    main()
