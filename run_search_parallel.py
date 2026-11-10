#!/usr/bin/env python3
"""Run the query and best-pair search trials in parallel worker processes.

Uses the same Q# expressions, seeds, and record order as
`run_experiments.py --task search`; only the scheduling differs. Every trial
sets its own quantum and classical seed before evaluation, so the trial
outcomes do not depend on the number of workers or on completion order.
Writes trials_search.jsonl and metadata_search.json; run
`run_experiments.py --task summarize` afterwards to rebuild the tables.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

from run_experiments import (BASE_SEED, DATABASE, EPSILON, NS, QUERY, ROOT,
                             containment, file_hashes, jaccard, literal)

SEARCH_OFFSET = 1000000
PHASE_BITS, MARGIN, ROUNDS_PER_UPDATE, MAX_UPDATES = 5, 0.12, 6, 4

_qsharp = None


def source_hashes():
    hashes = file_hashes()
    hashes["run_search_parallel.py"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return hashes


def build_jobs(trials):
    jobs = []
    for pair_mode in [False, True]:
        group = "Best-pair search" if pair_mode else "Query search"
        if pair_mode:
            optimum = max(jaccard(DATABASE[i], DATABASE[j])
                          for i in range(len(DATABASE)) for j in range(i + 1, len(DATABASE)))
        else:
            optimum = max(containment(a, QUERY) for a in DATABASE)
        expression = (f"SearchDatabase({literal(DATABASE)},{literal(QUERY)},{literal(pair_mode)},"
                      f"{PHASE_BITS},{MARGIN},{ROUNDS_PER_UPDATE},{MAX_UPDATES})")
        for trial in range(trials):
            qseed = BASE_SEED + SEARCH_OFFSET + len(jobs)
            jobs.append(dict(group=group, case="dna_database", trial=trial,
                             quantum_seed=qseed, classical_seed=qseed + 10000000,
                             expression=expression, expected=optimum))
    return jobs


def init_worker():
    global _qsharp
    from qdk import qsharp
    qsharp.init(project_root=str(ROOT))
    _qsharp = qsharp


def run_trial(job):
    _qsharp.set_quantum_seed(job["quantum_seed"])
    _qsharp.set_classical_seed(job["classical_seed"])
    started = time.perf_counter()
    result = _qsharp.eval(NS + job["expression"])
    return job, tuple(result), time.perf_counter() - started


def make_row(job, result, elapsed):
    observed, reported, work = result
    expected = job["expected"]
    row = dict(job, observed=observed, encoded_candidate=observed, reported_score=reported,
               logical_work=work, simulator_seconds=elapsed,
               metric="score regret", work_unit="A/A-dagger calls")
    if job["group"] == "Query search":
        valid = 0 <= observed < len(DATABASE)
        score = containment(DATABASE[observed], QUERY) if valid else -1.0
        row["decoded_reference"] = observed
    else:
        bits = max(1, (len(DATABASE) - 1).bit_length())
        i, j = observed & ((1 << bits) - 1), observed >> bits
        valid = 0 <= i < j < len(DATABASE)
        score = jaccard(DATABASE[i], DATABASE[j]) if valid else -1.0
        row["decoded_pair"] = [i, j]
    row.update(valid=valid, true_returned_score=score, error=max(0.0, expected - score),
               success=valid and score >= expected - EPSILON - 1e-12)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search-trials", type=int, default=32, help="Trials per search task")
    parser.add_argument("--workers", type=int, default=4, help="Worker processes")
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    if args.search_trials < 1 or args.workers < 1:
        parser.error("Trial and worker counts must be positive")
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)

    hashes_started = source_hashes()
    jobs = build_jobs(args.search_trials)
    metadata = {
        "task": "search", "runner": "run_search_parallel.py",
        "utc_started": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "platform": platform.platform(),
        "qdk_version": importlib.metadata.version("qdk"),
        "simulator": "QDK local noiseless sparse-state simulator; no quantum hardware",
        "scheduling": "ProcessPoolExecutor with multiprocessing spawn; one interpreter per worker",
        "workers": args.workers,
        "base_seed": BASE_SEED,
        "quantum_seed_formula": "BASE_SEED + 1000000 + record_index",
        "classical_seed_formula": "quantum_seed + 10000000",
        "record_order": "query trials followed by best-pair trials, each in increasing trial order",
        "epsilon_evaluation": EPSILON,
        "search_phase_bits": PHASE_BITS, "search_coherent_median_repetitions": 1,
        "search_threshold_median_repetitions": 3,
        "search_margin": MARGIN, "search_rounds_per_update": ROUNDS_PER_UPDATE,
        "search_max_updates": MAX_UPDATES,
        "trials_per_search_task": args.search_trials,
        "database": DATABASE, "query": QUERY,
        "source_sha256": hashes_started,
        "limitations": [
            "Explicit multiplexed QROM and exact key-comparison membership, not physical QRAM or cuckoo hashing.",
            "Search uses one coherent QAE register and a measured median-3 verification guard; it is a finite-budget variant, not the theorem's amplified parameter setting.",
            "Experiments validate small-instance functionality, not asymptotic speedup or genomic-scale performance.",
            "Search results concern repeated measurement seeds on one constructed four-reference DNA database.",
            "Parallel wall time concerns classical simulator scheduling and is not quantum algorithm running time.",
        ],
    }

    started = time.perf_counter()
    rows = []
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=context,
                             initializer=init_worker) as pool:
        # Submit the slower best-pair trials first for better load balance.
        futures = [pool.submit(run_trial, job) for job in reversed(jobs)]
        for future in as_completed(futures):
            row = make_row(*future.result())
            rows.append(row)
            print(f"search [{len(rows)}/{len(jobs)}]: {row['group']} trial {row['trial']}: "
                  f"success={row['success']} error={row['error']:.6g}", flush=True)
    rows.sort(key=lambda r: r["quantum_seed"])

    with (out / "trials_search.jsonl").open("w") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    hashes_finished = source_hashes()
    metadata.update(
        records=len(rows), expected_records=len(jobs), complete=len(rows) == len(jobs),
        successes=sum(r["success"] for r in rows),
        simulator_wall_seconds=time.perf_counter() - started,
        summed_trial_simulator_seconds=sum(r["simulator_seconds"] for r in rows),
        utc_finished=datetime.now(timezone.utc).isoformat(),
        sources_unchanged=hashes_started == hashes_finished,
        source_sha256_finished=hashes_finished,
    )
    (out / "metadata_search.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Completed search: {len(rows)} trials, {metadata['successes']} successes, "
          f"{metadata['simulator_wall_seconds']:.2f}s wall time", flush=True)


if __name__ == "__main__":
    main()
