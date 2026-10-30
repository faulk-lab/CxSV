#!/usr/bin/env python3
"""Per-node runner of the full panel build: works through this node's shard.

    python3 run_node.py BUILD NODE          (started by launch.sh, under setsid)

Runs on the node itself with the system python3 (standard library only).
Reads <BUILD>/nodes/<NODE>/runner.json (re-read every cycle, so `slots` can
be changed while it runs) and the node's rows of <BUILD>/plan/shard_plan.tsv,
in plan order (longest first). Per sample:

  1. stage  copy alignment + index from the NAS to stage_dir (node-local
            disk), one copy at a time per node (max_concurrent_stage), so the
            NAS is read once per sample instead of once per caller;
  2. run    one Snakemake instance (workflow/panel_build/Snakefile) inside
            the pinned image, --cores cores_per_sample, in the node's own
            working directory <BUILD>/nodes/<NODE>/work (--nolock: every
            instance writes only its own sample directory);
  3. done   when the sample's cost record exists
            (<results>/samples/<cohort>/<sample>/panel_build/cost.tsv); the
            staged copy is deleted. A failed sample is retried up to
            max_attempts, then left `failed` for a person to look at.

Staging may be throttled per node with runner.json stage_bwlimit_mb_s (0 =
unlimited; the NAS serves ~630 MB/s in total to all readers, lab work included).
Concurrency: at most `slots` samples run; a new one starts only while the
node's MemAvailable exceeds reserve_ram_gb + that sample's peak RSS, the 1-min
load average is below threads - reserve_threads, and (if set) inside
active_hours. That is how other lab work on a shared node keeps its share.

Control files in <BUILD>/nodes/<NODE>/ (scripts/panel_build/control.sh):
  RETRY   COHORT/SAMPLE per line: give gave_up samples another max_attempts
  PAUSE   start nothing new; running samples finish (drain). Remove to resume.
  STOP    terminate running samples now (Snakemake marks their outputs
          incomplete; the next run redoes them) and exit.
Resumable: re-running the runner skips done samples, re-adopts instances that
are still alive (pid files in run/), and redoes the rest; Snakemake's
--rerun-incomplete in the same node working directory cleans partial outputs.

Ledger: <BUILD>/nodes/<NODE>/ledger.tsv, append-only, written only by this
runner (one file per node, never shared): time, epoch, node, cohort,
sample_id, event, attempt, seconds, detail. Events: start_runner, staging,
staged, started, adopted, done, failed, retry, gave_up, paused, resumed,
stopped, exit_runner.
"""
from __future__ import annotations

import csv
import datetime as dt
import fcntl
import json
import os
import shutil
import signal
import socket
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import shardlib  # noqa: E402

POLL_S = 20
LEDGER_COLUMNS = ["time", "epoch", "node", "cohort", "sample_id", "event", "attempt", "seconds", "detail"]


class Runner:
    def __init__(self, build: str, node: str):
        self.build, self.node = build, node
        self.nd = f"{build}/nodes/{node}"
        self.cfg = {}
        self.reload()
        self.ledger_path = f"{self.nd}/ledger.tsv"
        for d in ("logs", "run", "overlays", "work"):
            os.makedirs(f"{self.nd}/{d}", exist_ok=True)
        self.plan = [r for r in shardlib.read_sheet(f"{build}/plan/shard_plan.tsv")[1] if r["node"] == node]
        self.plan.sort(key=lambda r: int(r["order"]))
        self.attempts = {}       # key -> attempts used
        self.running = {}        # key -> {"pid", "t0", "proc" or None}
        self.staging = {}        # key -> {"proc", "t0"}
        self.staged = set()
        self.failed = set()
        self.stage_failures = {}
        self.paused = False
        self.disk_short = False

    # ── config / ledger ────────────────────────────────────────────────────
    def reload(self):
        with open(f"{self.nd}/runner.json") as fh:
            self.cfg = json.load(fh)

    def log(self, key, event, attempt="", seconds="", detail=""):
        now = time.time()
        row = [dt.datetime.now().astimezone().isoformat(timespec="seconds"), f"{now:.0f}", self.node,
               key[0] if key else "", key[1] if key else "", event, str(attempt),
               f"{seconds:.0f}" if isinstance(seconds, float) else str(seconds), str(detail).replace("\t", " ")]
        new = not os.path.exists(self.ledger_path)
        with open(self.ledger_path, "a") as fh:
            if new:
                fh.write("\t".join(LEDGER_COLUMNS) + "\n")
            fh.write("\t".join(row) + "\n")
            fh.flush()
            os.fsync(fh.fileno())
        print("\t".join(row), flush=True)

    def previous_attempts(self):
        if not os.path.exists(self.ledger_path):
            return
        with open(self.ledger_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["event"] == "started" and r["sample_id"]:
                    k = (r["cohort"], r["sample_id"])
                    self.attempts[k] = max(self.attempts.get(k, 0), int(r["attempt"] or 1))
                if r["event"] == "gave_up":
                    self.failed.add((r["cohort"], r["sample_id"]))
                if r["event"] == "retry" and r["detail"] == "RETRY file":
                    k = (r["cohort"], r["sample_id"])
                    self.failed.discard(k)
                    self.attempts[k] = 0

    # ── paths ──────────────────────────────────────────────────────────────
    def sample_dir(self, key):
        return f"{self.cfg['results_root']}/samples/{key[0]}/{key[1]}"

    def is_done(self, key):
        return os.path.exists(f"{self.sample_dir(key)}/panel_build/cost.tsv")

    def tag(self, key):
        return f"{key[0]}__{key[1]}"

    def staged_paths(self, row):
        st = shardlib.staged_row(self.sheet_row(row), self.cfg["stage_dir"], row["cohort"])
        return st["alignment"], st["index"]

    def sheet_row(self, row):
        if not hasattr(self, "_sheets"):
            self._sheets = {}
        c = row["cohort"]
        if c not in self._sheets:
            hdr, rows = shardlib.read_sheet(f"{self.nd}/sheets/{c}.tsv")
            self._sheets[c] = (hdr, {r["sample_id"]: r for r in rows})
        return self._sheets[c][1][row["sample_id"]]

    # ── node state ─────────────────────────────────────────────────────────
    @staticmethod
    def mem_available_gb():
        with open("/proc/meminfo") as fh:
            for ln in fh:
                if ln.startswith("MemAvailable:"):
                    return int(ln.split()[1]) / 1048576
        return 0.0

    def in_active_hours(self):
        win = (self.cfg.get("active_hours") or "").strip()
        if not win:
            return True
        a, b = (int(x) for x in win.split("-"))
        h = dt.datetime.now().hour
        return a <= h < b if a < b else (h >= a or h < b)

    def stage_bytes_held(self):
        tot = 0
        for k in list(self.staged) + list(self.staging):
            r = self.row_of(k)
            tot += int(r["bytes"])
        return tot

    def row_of(self, key):
        for r in self.plan:
            if (r["cohort"], r["sample_id"]) == key:
                return r
        raise KeyError(key)

    # ── actions ────────────────────────────────────────────────────────────
    def start_stage(self, row):
        key = (row["cohort"], row["sample_id"])
        src = self.sheet_row(row)
        aln, idx = self.staged_paths(row)
        src_idx = src.get("index") or src["alignment"] + (".crai" if src["alignment"].endswith(".cram") else ".bai")
        cmd = (f"source '{HERE}/lib.sh' && cxsv_stage '{src['alignment']}' '{aln}' >/dev/null && "
               f"cxsv_stage '{src_idx}' '{idx}' >/dev/null")
        logf = open(f"{self.nd}/logs/{self.tag(key)}.stage.log", "a")
        env = dict(os.environ, CXSV_STAGE_BWLIMIT_MB=str(int(self.cfg.get("stage_bwlimit_mb_s", 0) or 0)))
        p = subprocess.Popen(["bash", "-c", cmd], stdout=logf, stderr=subprocess.STDOUT,
                             env=env, start_new_session=True)
        self.staging[key] = {"proc": p, "t0": time.time()}
        self.log(key, "staging", detail=f"{int(row['bytes']) / 1e9:.1f} GB")

    def write_overlay(self, row):
        key = (row["cohort"], row["sample_id"])
        st_row = shardlib.staged_row(self.sheet_row(row), self.cfg["stage_dir"], row["cohort"])
        hdr, _ = self._sheets[row["cohort"]]
        sheet = f"{self.nd}/overlays/{self.tag(key)}.tsv"
        shardlib.write_sheet(sheet, hdr, [st_row])
        node_ov = json.load(open(f"{self.nd}/node.yaml"))
        node_ov["panel_cohorts"] = {row["cohort"]: {"samples": sheet,
                                                    "platform": self.cfg["site_cohort_platform"][row["cohort"]]}}
        ov = f"{self.nd}/overlays/{self.tag(key)}.yaml"
        shardlib.dump_json(ov, node_ov)
        return ov

    def start_sample(self, row):
        key = (row["cohort"], row["sample_id"])
        attempt = self.attempts.get(key, 0) + 1
        self.attempts[key] = attempt
        ov = self.write_overlay(row)
        tmp = f"{self.cfg['tmp_dir']}/{self.tag(key)}"
        env = dict(os.environ, CXSV_SIF=self.cfg["sif"], CXSV_EXTRA_BINDS=self.cfg["stage_dir"])
        cmd = (f"source '{HERE}/lib.sh' && cxsv_snakemake '{self.nd}/work' '{tmp}' '{self.cfg['code']}' "
               f"{int(self.cfg['cores_per_sample'])} '{ov}'")
        logf = open(f"{self.nd}/logs/{self.tag(key)}.attempt{attempt}.log", "w")
        p = subprocess.Popen(["bash", "-c", cmd], stdout=logf, stderr=subprocess.STDOUT,
                             env=env, start_new_session=True)
        with open(f"{self.nd}/run/{self.tag(key)}.pid", "w") as fh:
            fh.write(f"{p.pid}\t{time.time():.0f}\t{attempt}\n")
        self.running[key] = {"pid": p.pid, "t0": time.time(), "proc": p}
        self.log(key, "started", attempt, detail=f"pid {p.pid}")

    def unstage(self, row):
        for f in self.staged_paths(row):
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
        self.staged.discard((row["cohort"], row["sample_id"]))

    def adopt(self):
        """Re-attach to instances a previous runner started that are still alive."""
        for fn in os.listdir(f"{self.nd}/run"):
            if not fn.endswith(".pid"):
                continue
            cohort, sample = fn[:-4].split("__", 1)
            key = (cohort, sample)
            pid, t0, attempt = open(f"{self.nd}/run/{fn}").read().split()
            alive = False
            try:
                os.kill(int(pid), 0)
                alive = sample in open(f"/proc/{pid}/cmdline").read().replace("\0", " ")
            except (OSError, FileNotFoundError):
                pass
            if alive:
                self.running[key] = {"pid": int(pid), "t0": float(t0), "proc": None}
                self.staged.add(key)
                self.log(key, "adopted", attempt, detail=f"pid {pid}")
            else:
                os.remove(f"{self.nd}/run/{fn}")

    def reap(self):
        for key, st in list(self.staging.items()):
            rc = st["proc"].poll()
            if rc is None:
                continue
            del self.staging[key]
            if rc == 0:
                self.staged.add(key)
                self.log(key, "staged", seconds=time.time() - st["t0"])
            else:
                self.stage_failures[key] = self.stage_failures.get(key, 0) + 1
                if self.stage_failures[key] < int(self.cfg.get("max_attempts", 2)):
                    self.log(key, "retry", detail=f"staging exit {rc}")
                else:
                    self.log(key, "gave_up", detail=f"staging exit {rc}; see logs/")
                    self.failed.add(key)
        for key, st in list(self.running.items()):
            if st["proc"] is not None:
                rc = st["proc"].poll()
                if rc is None:
                    continue
            else:
                try:
                    os.kill(st["pid"], 0)
                    continue
                except OSError:
                    rc = None        # adopted: exit code unknown, judged by outputs
            del self.running[key]
            try:
                os.remove(f"{self.nd}/run/{self.tag(key)}.pid")
            except FileNotFoundError:
                pass
            row = self.row_of(key)
            took = time.time() - st["t0"]
            if self.is_done(key):
                self.log(key, "done", self.attempts.get(key, ""), took, "")
                self.unstage(row)
            elif self.attempts.get(key, 0) < int(self.cfg.get("max_attempts", 2)):
                self.log(key, "retry", self.attempts.get(key, ""), took, f"exit {rc}")
            else:
                self.log(key, "gave_up", self.attempts.get(key, ""), took, f"exit {rc}; see logs/")
                self.failed.add(key)
                self.unstage(row)

    def stop_all(self):
        for key, st in list(self.running.items()) + list(self.staging.items()):
            pid = st["proc"].pid if st.get("proc") else st["pid"]
            try:
                os.killpg(pid, signal.SIGTERM)
            except OSError:
                pass
        deadline = time.time() + 120
        while time.time() < deadline and any(
                (s.get("proc") and s["proc"].poll() is None) for s in list(self.running.values()) + list(self.staging.values())):
            time.sleep(2)
        for key in list(self.running):
            self.log(key, "stopped", self.attempts.get(key, ""))

    # ── main loop ──────────────────────────────────────────────────────────
    def pending(self):
        busy = set(self.running) | set(self.staging) | self.failed
        return [r for r in self.plan if (r["cohort"], r["sample_id"]) not in busy
                and not self.is_done((r["cohort"], r["sample_id"]))]

    def can_start(self, row):
        c = self.cfg
        if self.mem_available_gb() < float(c["reserve_ram_gb"]) + float(row["est_peak_rss_gb"]):
            return False
        if os.getloadavg()[0] > int(c["threads"]) - int(c["reserve_threads"]):
            return False
        return True

    def take_retries(self):
        """RETRY file (control.sh retry): COHORT/SAMPLE per line -> back to pending."""
        f = f"{self.nd}/RETRY"
        if not os.path.exists(f):
            return
        for ln in open(f).read().split():
            key = tuple(ln.strip().split("/", 1))
            if key in self.failed:
                self.failed.discard(key)
                self.attempts[key] = 0
                self.stage_failures.pop(key, None)
                self.log(key, "retry", detail="RETRY file")
        os.remove(f)

    def cycle(self):
        self.reload()
        self.take_retries()
        self.reap()
        if os.path.exists(f"{self.nd}/STOP"):
            self.stop_all()
            return False
        pause = os.path.exists(f"{self.nd}/PAUSE") or not self.in_active_hours()
        if pause != self.paused:
            self.paused = pause
            self.log(None, "paused" if pause else "resumed",
                     detail="PAUSE file" if os.path.exists(f"{self.nd}/PAUSE") else "active_hours")
        pend = self.pending()
        if not pend and not self.running and not self.staging:
            return False
        if pause:
            return True
        # start staged samples first, in plan order
        for row in pend:
            key = (row["cohort"], row["sample_id"])
            if len(self.running) >= int(self.cfg["slots"]):
                break
            if key in self.staged and self.can_start(row):
                self.start_sample(row)
        # stage ahead: keep up to (free slots + 1) staged, within the disk cap
        want = int(self.cfg["slots"]) - len(self.running) + 1
        ready = len([k for k in self.staged if k not in self.running])
        for row in self.pending():
            key = (row["cohort"], row["sample_id"])
            if key in self.staged or key in self.staging:
                continue
            if len(self.staging) >= int(self.cfg.get("max_concurrent_stage", 1)) or ready >= want:
                break
            if (self.stage_bytes_held() + int(row["bytes"])) / 1e9 > float(self.cfg["stage_max_gb"]):
                break
            free = shutil.disk_usage(self.cfg["stage_dir"]).free
            short = free < int(row["bytes"]) * 1.1 + 50e9
            if short != self.disk_short:
                self.disk_short = short
                self.log(key, "paused" if short else "resumed", detail=f"stage disk has {free / 1e9:.0f} GB free")
            if short:
                break
            self.start_stage(row)
            ready += 1
        return True

    def main(self):
        os.makedirs(self.cfg["stage_dir"], exist_ok=True)
        os.makedirs(self.cfg["tmp_dir"], exist_ok=True)
        lock_path = f"{self.cfg['stage_dir']}/.runner.lock"      # node-local: flock is reliable there
        lock = open(lock_path, "w")
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            sys.exit(f"another run_node.py holds {lock_path} on {socket.gethostname()}")
        self.previous_attempts()
        self.log(None, "start_runner", detail=f"host {socket.gethostname()} pid {os.getpid()} "
                                              f"plan {len(self.plan)} samples")
        self.adopt()
        # copies left by an earlier runner are reused if complete (cxsv_stage checks sizes)
        for row in self.plan:
            key = (row["cohort"], row["sample_id"])
            aln, idx = self.staged_paths(row)
            src = self.sheet_row(row)
            if key not in self.running and os.path.exists(aln) and os.path.exists(idx) \
                    and os.path.getsize(aln) == os.path.getsize(src["alignment"]) and not self.is_done(key):
                self.staged.add(key)
        while self.cycle():
            time.sleep(POLL_S)
        done = sum(self.is_done((r["cohort"], r["sample_id"])) for r in self.plan)
        self.log(None, "exit_runner", detail=f"{done}/{len(self.plan)} done, {len(self.failed)} failed")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    Runner(sys.argv[1], sys.argv[2]).main()
