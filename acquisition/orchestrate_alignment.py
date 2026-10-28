#!/usr/bin/env python3
"""Distribute CHM13 alignment of the 1kGP-ONT cohort across the lab nodes.

Run from freya, which can reach the others over the passwordless ssh mesh:

    python3 orchestrate_alignment.py --manifest .../manifest_PRJEB89727_odp.tsv

Why a dispatcher rather than GNU parallel: the nodes are not interchangeable.
athena logs in as a different user, shadowfax has 61 GB of RAM against 125 on the
others and needs a smaller minimap2 minibatch, and work can only start on a
sample whose download has already been *verified* - which changes while this is
running, because the download is still going. A small queue loop expresses that;
a static job list does not.

Safety properties worth keeping if this is edited:
* Only samples marked ``verified`` in the download ledger are eligible, so a
  half-transferred .sra is never aligned.
* A sample already carrying an indexed CRAM is skipped, making the whole run
  resumable after an interruption.
* Failures are recorded and retried a bounded number of times rather than
  silently dropped, and a sample that fails repeatedly is reported at the end.
"""

from __future__ import annotations

import argparse
import csv
import os
import signal
import socket
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

#: Per-node tuning. threads/batch reflect real RAM, not core count: shadowfax has
#: 61 GB, so a 10g minibatch alongside an 8.5 GB index and the sort buffer is too
#: close to the edge.
NODES = {
    "freya":     dict(slots=1, threads_mm=24, threads_sort=8, sort_mem="2G", batch="4g",
                      home="/home/minknow", local=True, env="sratools"),
    "shadowfax": dict(slots=1, threads_mm=20, threads_sort=6, sort_mem="1G", batch="4g",
                      home="/home/minknow", local=False, env="sratools"),
    # athena returned to service 2026-10-22 after a failed DIMM was removed. It
    # now reports 94 GB rather than the 125 GB it had before, so the old tuning
    # was NOT carried over unexamined; 4g keeps peak RSS well inside 94 GB.
    "athena":    dict(slots=1, threads_mm=24, threads_sort=8, sort_mem="2G", batch="4g",
                      home="/home/tia", local=False, env="sratools"),
    # betsy shipped minimap2 2.27; the pinned cxsvalign env holds 2.24 so its
    # alignments match the rest of the panel. It also has a history of DRAM
    # faults, which is a further reason the read-count integrity gate must pass
    # before any .sra it processed is deleted.
    "betsy":     dict(slots=1, threads_mm=24, threads_sort=8, sort_mem="2G", batch="4g",
                      home="/home/minknow", local=False, env="cxsvalign"),
    # princess is the odd one out: 12 cores and only 29 GB RAM, against 61-126 GB
    # elsewhere, and it shipped minimap2 2.31 so it runs from the pinned
    # cxsvalign env like betsy.
    #
    # Settings are MEASURED, not projected. -K 1g with 10 threads was tried first
    # and put minimap2 at 19.0 GB RSS with the machine 24/29 GB used and 3 GB into
    # swap - running, but with no headroom for the sort phase. The 9.0 GB index is
    # irreducible fixed cost here, so the batch carries the cut. Note minimap2
    # memory tracks -K and thread count, NOT sample depth, so these hold for the
    # deepest samples too. Expect roughly 3x freya's wall time per sample.
    "princess":  dict(slots=1, threads_mm=8, threads_sort=3, sort_mem="512M", batch="400m",
                      home="/home/minknow", local=False, env="cxsvalign"),
    # keurig: 32 cores / 60 GB, i.e. shadowfax-class memory, so it takes
    # shadowfax's envelope rather than freya's. It shipped minimap2 2.27 and so
    # runs from the pinned cxsvalign env. Two things had to be fixed before it
    # could join: its NAS shares were unmounted (fstab was right but needed a
    # root `mount -a`), and freya had no host key for it - installed only after
    # confirming the fingerprint read off keurig itself matched what freya's
    # ssh-keyscan returned, so it was verified rather than trust-on-first-use.
    "keurig":    dict(slots=1, threads_mm=20, threads_sort=6, sort_mem="1G", batch="2g",
                      home="/home/minknow", local=False, env="cxsvalign"),
}

ROOT = "/mnt/sargasso1/CxSV"
SCRIPT = f"{ROOT}/code/align_sra_to_chm13.sh"
SRA_DIR = f"{ROOT}/datasets/1kGP_ONT_1019_sra"
OUT_DIR = f"{ROOT}/datasets/1kGP_ONT_1019_CHM13"
LOG_DIR = f"{ROOT}/datasets/logs"
CLAIM_DIR = f"{ROOT}/datasets/logs/claims"

_print_lock = threading.Lock()
_ledger_lock = threading.Lock()

#: Set by SIGTERM/SIGINT. Workers finish the alignment they hold, then exit.
_stop = threading.Event()


def _handle_signal(signum, _frame):
    """Stop scheduling new work but let running alignments finish.

    Killing the orchestrator outright orphans its ssh children: the alignments
    survive and complete, but no ledger row is written for them. Draining is
    cleaner, so SIGTERM and SIGINT are handled rather than fatal.
    """
    _stop.set()
    log(f"signal {signum} received - draining: no new samples will be dispatched")


def drain_requested(node: str, log_dir: str) -> bool:
    """True when <log_dir>/drain.<node> exists.

    Lets one node be taken out of the pool for a reboot with a single `touch`,
    without restarting the run and without wasting the sample it is mid-way
    through. Remove the file to let the node take work again.
    """
    from pathlib import Path as _P
    return _P(log_dir, f"drain.{node}").exists()


def log(msg: str) -> None:
    with _print_lock:
        print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def record(path: Path, row: list) -> None:
    with _ledger_lock:
        new = not path.exists()
        with open(path, "a", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            if new:
                w.writerow(["ts", "sample", "accession", "node", "status",
                            "wall_s", "bytes", "note"])
            w.writerow(row)


def verified_downloads(ledger: Path) -> set[str]:
    """Accessions whose .sra transfer the downloader marked verified."""
    if not ledger.exists():
        return set()
    out = set()
    with open(ledger, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] == "verified":
                out.add(r["filename"].replace(".sra", ""))
    return out


def claim(acc: str, claim_dir: str) -> bool:
    """Atomically claim ``acc`` for this process. True if we got it.

    O_CREAT|O_EXCL is atomic over SMB/CIFS (it maps to a FILE_CREATE
    disposition), so two orchestrator instances - or a restart overlapping its
    predecessor - cannot both take the same sample. Without this, two instances
    sort the queue identically, pick the same top sample, and write the same
    ``.cram.tmp`` path.
    """
    os.makedirs(claim_dir, exist_ok=True)
    path = os.path.join(claim_dir, f"{acc}.claim")
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        return False
    with os.fdopen(fd, "w") as fh:
        fh.write(f"{os.getpid()} {datetime.now(timezone.utc).isoformat()}\n")
    return True


def release_claim(acc: str, claim_dir: str) -> None:
    """Drop a claim so a failed sample can be retried by anyone."""
    try:
        os.unlink(os.path.join(claim_dir, f"{acc}.claim"))
    except FileNotFoundError:
        pass


def reap_stale_claims(claim_dir: str, nodes: dict, max_age_s: int = 4 * 3600) -> int:
    """Release claims that no longer correspond to work in progress.

    A claim whose sample has no CRAM and whose aligning process has gone - a
    node rebooted mid-sample, a killed orchestrator - would otherwise strand
    that sample forever: claimed, never aligned, and silently missing from the
    panel. Only claims older than ``max_age_s`` are considered, so a claim held
    by a legitimately long-running alignment is never pulled out from under it.

    A claim is kept if EITHER the CRAM exists (work finished) or the accession
    is in some node's live process list (work ongoing). ``nodes`` must therefore
    be the FULL fleet, not just this instance's dispatch set - a claim held by
    work on an un-probed node would otherwise be reaped while still live.
    """
    from pathlib import Path as _P
    now = time.time()
    running: set[str] = set()
    for node, cfg in nodes.items():
        # Match the accession from either the minimap2 read group (mapping
        # phase) or the align script's own command line (which persists through
        # the sort/index/QC phases after minimap2 has exited).
        probe = ("{ ps -eo args --no-headers | grep -- '-x map-ont' | "
                 "grep -oE 'ID:ERR[0-9]+' | cut -d: -f2; "
                 "ps -eo args --no-headers | grep '[a]lign_sra_to_chm13.sh' | "
                 "grep -oE 'ERR[0-9]+' ; } 2>/dev/null | sort -u || true")
        cmd = ["bash", "-lc", probe] if cfg["local"] else [
            "ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", node, probe]
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout
            running.update(x.strip() for x in out.split() if x.strip())
        except Exception:
            # Treat an unreachable node as possibly-busy: do not reap on its behalf.
            return 0
    freed = 0
    for cf in _P(claim_dir).glob("*.claim"):
        acc = cf.stem
        if acc in running:
            continue
        if now - cf.stat().st_mtime < max_age_s:
            continue
        if list(_P(OUT_DIR).glob(f"*.{acc}.*.cram")):
            continue
        # A .tmp means the alignment is past minimap2 and inside samtools sort /
        # index / QC, which on a 15 GB CRAM runs for many minutes. The process
        # probe above only sees minimap2, so without this a sample in its sort
        # phase looks idle and its claim gets reaped out from under live work.
        if list(_P(OUT_DIR).glob(f"*.{acc}.*.cram.tmp")):
            continue
        cf.unlink()
        freed += 1
        log(f"  reaped stale claim for {acc} (no CRAM, no running process)")
    return freed


def node_busy(node: str, cfg: dict) -> bool:
    """True while an alignment is already running on ``node``.

    A restart can overlap alignments orphaned by its predecessor; dispatching
    into that would put two 24-thread jobs on one machine. Checking the node
    itself is more reliable than tracking what we think we started.
    """
    # samtools is included deliberately: an alignment spends its sort, index and
    # flagstat phases with no minimap2 process alive, and dispatching into that
    # window would put two jobs on the node.
    probe = "pgrep -c 'minimap2|fasterq-dump|samtools' || true"
    cmd = ["bash", "-lc", probe] if cfg["local"] else [
        "ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", node, probe]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip()
        return int(out.split()[0] or 0) > 0
    except Exception:
        return False


def already_aligned(acc: str, sample: str) -> bool:
    """True if a final CRAM exists for this sample.

    Deliberately keyed on the CRAM, not its index. An interrupted run can leave
    a complete CRAM whose .crai never got written; treating that as "not done"
    re-aligned good samples for hours and raced the verifier on the same file.
    A missing index is a repair job, not a reason to redo the alignment.
    """
    cram = Path(OUT_DIR) / f"{sample}.{acc}.hs1.minimap2.cram"
    return cram.exists() and cram.stat().st_size > 0


def run_one(node: str, cfg: dict, sample: str, acc: str, ledger: Path) -> bool:
    env = (
        f"REF_DIR={cfg['home']}/refs/hs1 "
        f"SCRATCH={cfg['home']}/scratch/cxsv_align "
        f"OUT_DIR={OUT_DIR} LOG_DIR={LOG_DIR}/align "
        f"THREADS_MM={cfg['threads_mm']} THREADS_SORT={cfg['threads_sort']} "
        f"SORT_MEM={cfg['sort_mem']} MM_BATCH={cfg['batch']} "
        f"CONDA_ENV={cfg['env']}"
    )
    inner = f"{env} bash {SCRIPT} {SRA_DIR}/{acc}.sra {sample}"
    cmd = ["bash", "-lc", inner] if cfg["local"] else [
        "ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", node, inner]

    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    wall = int(time.time() - t0)
    tail = (p.stdout or "").strip().split("\n")[-1][:200]

    if p.returncode == 0:
        size = ""
        for tok in tail.split():
            if tok.isdigit() and len(tok) > 8:
                size = tok
        record(ledger, [datetime.now(timezone.utc).isoformat(), sample, acc, node,
                        "aligned", wall, size, tail])
        log(f"  OK   {sample}/{acc} on {node} in {wall}s")
        return True

    err = (p.stderr or "").strip().split("\n")[-1][:200]

    # 126 = found but not executable, 127 = not found, or a stale NAS handle on
    # the script itself: problems with the node, not with the sample.
    if p.returncode in (126, 127) or "Stale file handle" in (p.stderr or ""):
        record(ledger, [datetime.now(timezone.utc).isoformat(), sample, acc, node,
                        "node_error", wall, "", f"rc={p.returncode} {err}"])
        return "node_error"

    # rc 255 is ssh's own exit code - the transport died, which says nothing
    # about whether the alignment finished. Observed in practice: a "Broken
    # pipe" drop after minimap2 and the sort had completed, leaving a valid
    # 9.3 GB CRAM that the ledger called a failure. Check the output before
    # believing the return code.
    if p.returncode == 255:
        cram = Path(OUT_DIR) / f"{sample}.{acc}.hs1.minimap2.cram"
        if cram.exists() and cram.stat().st_size > 0:
            verdict = repair_after_transport_drop(node, cfg, sample, acc, cram)
            record(ledger, [datetime.now(timezone.utc).isoformat(), sample, acc, node,
                            verdict, wall, cram.stat().st_size,
                            f"ssh transport drop; output {verdict}"])
            log(f"  {verdict.upper()} {sample}/{acc} on {node} (ssh dropped, output checked)")
            if verdict == "unverified":
                # Keep the claim: a CRAM exists and nothing showed it is bad.
                # Re-aligning would only overwrite it. Leave it for review.
                return True
            if verdict == "bad":
                # The check ran and the CRAM failed it. already_aligned() keys
                # on the CRAM existing, so it must leave the panel or the sample
                # would count as done forever. Moved, not deleted.
                qdir = Path(OUT_DIR).parent / "quarantine"
                qdir.mkdir(exist_ok=True)
                for f in (cram, Path(f"{cram}.crai")):
                    if f.exists():
                        f.rename(qdir / f.name)
                log(f"  quarantined failed CRAM {cram.name}; sample will be re-aligned")
                return False
            return verdict == "aligned"

    record(ledger, [datetime.now(timezone.utc).isoformat(), sample, acc, node,
                    "failed", wall, "", f"rc={p.returncode} {err}"])
    log(f"  FAIL {sample}/{acc} on {node} rc={p.returncode}: {err}")
    return False


INTEGRITY_DIR = f"{ROOT}/datasets/logs/align/integrity"


def write_integrity_record(sample: str, acc: str, sra_n: str, cram_n: str,
                           verdict: str, size: int) -> None:
    """Write one sample's integrity record atomically.

    Never append to a shared file here: several hosts write integrity records
    concurrently, and O_APPEND is not atomic across SMB clients. The shared
    align_integrity.tsv was corrupted exactly that way (NUL-filled rows, two
    rows lost) once six nodes were appending to it. Write to a temp name and
    rename(); rename is atomic on SMB, so concurrent writers cannot clobber.
    """
    d = Path(INTEGRITY_DIR)
    d.mkdir(parents=True, exist_ok=True)
    tmp = d / f".{acc}.{os.getpid()}.tmp"
    tmp.write_text(f"{sample}\t{acc}\t{sra_n}\t{cram_n}\t{verdict}\t{size}\n")
    os.replace(tmp, d / f"{acc}.tsv")


def repair_after_transport_drop(node: str, cfg: dict, sample: str, acc: str,
                                cram: Path) -> str:
    """Verify (and finish) an alignment whose ssh session died late.

    Returns one of three verdicts - the distinction matters:
      "aligned"      CRAM passes quickcheck AND its primary-alignment count
                     equals the read count the .sra reports.
      "bad"          the check RAN and the CRAM failed it.
      "unverified"   the check could not be completed. This is NOT evidence of
                     a bad CRAM and must never trigger re-alignment.

    The check runs on the machine running the dispatcher, never over ssh to
    the node that just dropped. The CRAM, the .sra and the reference are all on
    the shared NAS, so the check does not need that node at all. The first
    version ran it over a second ssh to the same flaky link; when that also
    dropped, four good CRAMs were reported as failures and re-queued.
    """
    ref = None
    for cand in ("/home/minknow/refs/hs1/hs1.fa", f"{ROOT}/refs/hs1.fa"):
        if Path(cand).exists():
            ref = cand
            break
    if ref is None:
        return "unverified"
    sra = f"{SRA_DIR}/{acc}.sra"
    env = "source ~/miniforge3/etc/profile.d/conda.sh && conda activate sratools && "

    def sh(cmd: str, timeout: int = 3600):
        return subprocess.run(["bash", "-lc", env + cmd], capture_output=True,
                              text=True, timeout=timeout)
    try:
        if sh(f"samtools quickcheck {cram}").returncode != 0:
            return "bad"
        crai = Path(f"{cram}.crai")
        if not crai.exists() or crai.stat().st_size == 0:
            if sh(f"samtools index -@ 4 {cram}").returncode != 0:
                return "unverified"
        p1 = sh(f"sra-stat --quick {sra} | awk -F'|' '{{print $3}}' | awk -F: '{{print $1}}'")
        p2 = sh(f"samtools view -c -F 0x900 --reference {ref} {cram}")
        sra_n, cram_n = p1.stdout.strip(), p2.stdout.strip()
        if p1.returncode != 0 or p2.returncode != 0 or not sra_n.isdigit() or not cram_n.isdigit():
            return "unverified"
    except subprocess.TimeoutExpired:
        return "unverified"
    if sra_n != cram_n:
        return "bad"
    write_integrity_record(sample, acc, sra_n, cram_n, "read_count_match",
                           cram.stat().st_size)
    return "aligned"

def worker(node: str, cfg: dict, queue: list, lock: threading.Lock,
           ledger: Path, dl_ledger: Path, attempts: dict, max_retries: int,
           stop: threading.Event) -> None:
    """One slot on one node: take the next ready sample, align it, repeat."""
    while not stop.is_set() and not _stop.is_set():
        if drain_requested(node, f"{LOG_DIR}"):
            log(f"  {node}: drain file present - leaving the pool")
            return
        with lock:
            ready = verified_downloads(dl_ledger)
            pick = None
            for i, (sample, acc) in enumerate(queue):
                if acc in ready:
                    pick = queue.pop(i)
                    break
            if pick is None:
                if not queue:
                    return
                # Nothing downloaded yet that isn't claimed; the download is
                # still running, so wait rather than exit.
                pass
        if pick is None:
            time.sleep(60)
            continue

        sample, acc = pick
        if already_aligned(acc, sample):
            log(f"  skip {sample}/{acc} (CRAM already present)")
            continue
        if not claim(acc, CLAIM_DIR):
            # Requeue rather than drop. A claim can be momentary - a dispatcher
            # that claims and instantly fails a sample releases it again - and
            # dropping on sight would leave that sample unaligned for the rest
            # of this run. Re-check later; already_aligned() drops it once done.
            if not already_aligned(acc, sample):
                with lock:
                    queue.append((sample, acc))
            time.sleep(1)
            continue

        # Wait out any alignment this node is still running from a previous
        # instance rather than stacking a second one on top of it.
        waited = 0
        while node_busy(node, cfg) and not _stop.is_set():
            if waited == 0:
                log(f"  {node}: busy with an earlier alignment, waiting")
            time.sleep(60)
            waited += 60
        if _stop.is_set():
            release_claim(acc, CLAIM_DIR)
            return

        ok = run_one(node, cfg, sample, acc, ledger)
        if ok == "node_error":
            # The node could not run the script at all (stale NAS handle,
            # missing interpreter). That says nothing about the sample, so it
            # must not count against it - and retrying on this node would only
            # repeat the error at full speed. keurig once logged 2,615 such
            # "failures" in minutes after a script redeploy left a stale handle.
            release_claim(acc, CLAIM_DIR)
            with lock:
                queue.insert(0, (sample, acc))
            log(f"  {node}: NODE ERROR running the align script - leaving the pool")
            return
        if not ok:
            release_claim(acc, CLAIM_DIR)
            attempts[acc] = attempts.get(acc, 0) + 1
            if attempts[acc] <= max_retries:
                with lock:
                    queue.append((sample, acc))
                log(f"  requeued {sample}/{acc} (attempt {attempts[acc]})")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", required=True)
    p.add_argument("--download-ledger", default=f"{LOG_DIR}/odp_download.tsv")
    # Per-host by default: a dispatcher on keurig and the ones on freya must not
    # append to the same file (O_APPEND is not atomic across SMB clients).
    p.add_argument("--ledger",
                   default=f"{LOG_DIR}/alignment_ledger.{socket.gethostname()}.tsv")
    p.add_argument("--only-passing", action="store_true",
                   help="restrict to samples at or above the coverage floor")
    p.add_argument("--nodes", default="freya,shadowfax,betsy")
    p.add_argument("--max-retries", type=int, default=2)
    p.add_argument("--limit", type=int, default=0)
    a = p.parse_args()

    with open(a.manifest, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if a.only_passing:
        rows = [r for r in rows if r["passes_coverage_floor"] == "yes"]
    # Highest coverage first: the most informative background genomes land early,
    # so an interrupted run still leaves the best data aligned.
    rows.sort(key=lambda r: -float(r["est_coverage_x"]))

    queue = [(r["individual_id"] or r["run_accession"], r["run_accession"]) for r in rows]
    queue = [q for q in queue if not already_aligned(q[1], q[0])]
    if a.limit:
        queue = queue[: a.limit]

    log(f"{len(rows)} samples in scope, {len(queue)} to align")
    chosen = [n.strip() for n in a.nodes.split(",") if n.strip() in NODES]
    log(f"nodes: {', '.join(chosen)}")

    # Probe EVERY known node, not just the ones this instance dispatches to. A
    # claim can be held by an alignment running on a node outside --nodes (a
    # hand-launched validation run, or another instance's node); probing only
    # `chosen` would call that claim stale and hand the sample to a second node,
    # which is the collision this mechanism exists to prevent.
    freed = reap_stale_claims(CLAIM_DIR, NODES)
    if freed:
        log(f"released {freed} stale claim(s) before starting")

    lock = threading.Lock()
    stop = threading.Event()
    attempts: dict[str, int] = {}
    threads = []
    for node in chosen:
        cfg = NODES[node]
        for _ in range(cfg["slots"]):
            t = threading.Thread(target=worker, args=(node, cfg, queue, lock,
                                                      Path(a.ledger), Path(a.download_ledger),
                                                      attempts, a.max_retries, stop),
                                 daemon=True)
            t.start()
            threads.append(t)
    # Treat the machine we are running on as local regardless of how the table
    # declares it. This lets a dispatcher run ON a node for itself, which is the
    # only route to a node whose inbound ssh is gated behind an interactive
    # check (keurig sits behind Tailscale SSH's re-auth prompt, so freya cannot
    # reach it, but keurig can run its own dispatcher over the shared NAS).
    _here = socket.gethostname()
    for _n, _c in NODES.items():
        _c["local"] = (_n == _here)
    log(f"running on {_here}; local nodes: "
        f"{[n for n, cc in NODES.items() if cc['local']] or 'none'}")

    signal.signal(signal.SIGTERM, _handle_signal)
    signal.signal(signal.SIGINT, _handle_signal)

    # join() with no timeout is not reliably interruptible in CPython: a signal
    # arriving while the main thread blocks there may never raise. Polling with
    # a timeout is what makes SIGTERM/SIGINT drain instead of kill.
    try:
        while any(t.is_alive() for t in threads):
            for t in threads:
                t.join(timeout=1.0)
            if _stop.is_set() and all(not t.is_alive() for t in threads):
                break
    except KeyboardInterrupt:
        _stop.set()
        log("stopping after current alignments")

    failed = {k: v for k, v in attempts.items() if v > a.max_retries}
    log(f"done. persistently failed: {len(failed)}")
    for acc, n in failed.items():
        log(f"  {acc}: {n} attempts")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
