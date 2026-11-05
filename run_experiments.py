#!/usr/bin/env python3
"""Run actual Q# simulations and regenerate the paper table from trial records.

Python extracts sets, computes independent classical ground truth, schedules seeded
Q# runs, and summarizes outcomes. Quantum algorithms receive input sets and hash
seeds only, never precomputed intersections, similarity scores, or minima.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import platform
import random
import statistics
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
EPSILON = 0.1
BASE_SEED = 20261021
NS = "QuantumSketchExperiments."


def kmers(sequence: str, k: int = 2) -> list[int]:
    alphabet = {b: i for i, b in enumerate("ACGT")}
    result = set()
    for i in range(len(sequence) - k + 1):
        value = 0
        for base in sequence[i:i+k]:
            value = 4 * value + alphabet[base]
        result.add(value)
    return sorted(result)


SEQUENCES = ["AACCA", "ACGTA", "ACGTT", "CCGGC"]
DATABASE = [kmers(x) for x in SEQUENCES]
QUERY = DATABASE[2]


def containment(a, b):
    return len(set(a) & set(b)) / len(a)


def jaccard(a, b):
    return len(set(a) & set(b)) / len(set(a) | set(b))


def literal(x):
    return json.dumps(x, separators=(",", ":"))


def wilson(successes, trials):
    z = 1.959963984540054
    p = successes / trials
    den = 1 + z*z/trials
    center = (p + z*z/(2*trials)) / den
    radius = z * math.sqrt(p*(1-p)/trials + z*z/(4*trials*trials)) / den
    return [max(0.0, center-radius), min(1.0, center+radius)]


def build_instances():
    rng = random.Random(BASE_SEED)
    hash_cases = [{"name": "n4", "keys": DATABASE[1], "modulus": 32},
                  {"name": "n8", "keys": sorted(rng.sample(range(32), 8)), "modulus": 32},
                  {"name": "n16", "keys": sorted(rng.sample(range(32), 16)), "modulus": 32}]
    for item in hash_cases:
        item.update(multipliers=[3, 5, 7], offsets=[1, 2, 4])
    return {
        "k": 2, "encoding": "A=0,C=1,G=2,T=3; base-4 integer encoding",
        "sequences": SEQUENCES, "database": DATABASE, "query": QUERY,
        "containment_cases": [
            {"name": f"overlap_{c}_of_4", "source": [0, 1, 2, 3], "target": list(range(c))}
            for c in range(5)
        ] + [
            {"name": "non_power_of_two", "source": [0, 1, 2], "target": [1]},
            {"name": "dna_pair", "source": DATABASE[1], "target": DATABASE[2]},
        ],
        "hash_cases": hash_cases,
    }


def file_hashes():
    paths = list((ROOT / "src").glob("*.qs")) + [ROOT / "run_experiments.py"]
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", choices=["all", "estimators", "search", "minhash", "summarize"], default="all")
    parser.add_argument("--trials", type=int, default=64, help="Trials per estimation/MinHash instance")
    parser.add_argument("--search-trials", type=int, default=32, help="Trials per search task")
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    if min(args.trials, args.search_trials) < 1:
        parser.error("Trial counts must be positive")
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    instances = build_instances()
    if args.task != "summarize":
        (out / "instances.json").write_text(json.dumps(instances, indent=2) + "\n")
        from qdk import qsharp
        qsharp.init(project_root=str(ROOT))
        tasks = ["estimators", "search", "minhash"] if args.task == "all" else [args.task]
        for task in tasks:
            metadata = {
                "task": task, "utc_started": datetime.now(timezone.utc).isoformat(),
                "python": sys.version, "platform": platform.platform(),
                "qdk_version": importlib.metadata.version("qdk"),
                "simulator": "QDK local noiseless sparse-state simulator; no quantum hardware",
                "base_seed": BASE_SEED, "epsilon_evaluation": EPSILON,
                "estimation_phase_bits": 6, "estimation_median_repetitions": 3,
                "search_phase_bits": 5, "search_coherent_median_repetitions": 1,
                "search_threshold_median_repetitions": 3,
                "search_margin": 0.12, "search_rounds_per_update": 6, "search_max_updates": 4,
                "minhash_restarts_per_coordinate": 2,
                "trials_per_estimation_or_hash_case": args.trials,
                "trials_per_search_task": args.search_trials,
                "source_sha256": file_hashes(),
                "limitations": [
                    "Explicit multiplexed QROM and exact key-comparison membership, not physical QRAM or cuckoo hashing.",
                    "Search uses one coherent QAE register and a measured median-3 verification guard; it is a finite-budget variant, not the theorem's amplified parameter setting.",
                    "Experiments validate small-instance functionality, not asymptotic speedup or genomic-scale performance.",
                    "Search results concern repeated measurement seeds on one constructed four-reference DNA database.",
                ],
            }
            started = time.perf_counter()
            path = out / f"trials_{task}.jsonl"
            with path.open("w", buffering=1) as handle:
                record_number = 0
                def run(group, case, trial, expression, expected, extra=None):
                    nonlocal record_number
                    task_offset = {"estimators": 0, "search": 1000000, "minhash": 2000000}[task]
                    qseed = BASE_SEED + task_offset + record_number
                    cseed = qseed + 10000000
                    qsharp.set_quantum_seed(qseed)
                    qsharp.set_classical_seed(cseed)
                    t = time.perf_counter()
                    result = qsharp.eval(NS + expression)
                    elapsed = time.perf_counter() - t
                    observed, work = result[0], result[-1]
                    row = dict(group=group, case=case, trial=trial, quantum_seed=qseed,
                               classical_seed=cseed, expression=expression, expected=expected,
                               observed=observed, logical_work=work, simulator_seconds=elapsed)
                    if extra:
                        row.update(extra)
                    if group in ["Containment QAE", "Jaccard via QAE"]:
                        row.update(error=abs(observed-expected), success=abs(observed-expected) <= EPSILON + 1e-12,
                                   metric="absolute error", work_unit="A/A-dagger calls")
                    elif group in ["Query search", "Best-pair search"]:
                        if group == "Query search":
                            valid = 0 <= observed < len(DATABASE)
                            score = containment(DATABASE[observed], QUERY) if valid else -1.0
                        else:
                            bits = max(1, (len(DATABASE)-1).bit_length())
                            i, j = observed & ((1 << bits)-1), observed >> bits
                            valid = 0 <= i < j < len(DATABASE)
                            score = jaccard(DATABASE[i], DATABASE[j]) if valid else -1.0
                            row["decoded_pair"] = [i, j]
                        row.update(valid=valid, true_returned_score=score, reported_score=result[1],
                                   error=max(0.0, expected-score), success=valid and score >= expected-EPSILON-1e-12,
                                   metric="score regret", work_unit="A/A-dagger calls")
                    else:
                        observed = list(observed)
                        row.update(observed=observed, error=sum(x != y for x, y in zip(observed, expected)),
                                   success=observed == expected, metric="incorrect coordinates",
                                   work_unit="hash compute/uncompute calls")
                    handle.write(json.dumps(row, sort_keys=True) + "\n")
                    record_number += 1
                    if trial == 0 or (trial+1) % 8 == 0:
                        print(f"{task}: {case} {trial+1}: success={row['success']} error={row['error']:.6g}", flush=True)
                if task == "estimators":
                    for case in instances["containment_cases"]:
                        a, b = case["source"], case["target"]
                        for trial in range(args.trials):
                            run("Containment QAE", case["name"], trial,
                                f"EstimateContainment({literal(a)},{literal(b)},6,3)", containment(a, b))
                    a, b = DATABASE[1:3]
                    for trial in range(args.trials):
                        run("Jaccard via QAE", "dna_pair", trial,
                            f"EstimateJaccard({literal(a)},{literal(b)},6,3)", jaccard(a, b))
                elif task == "search":
                    for pair_mode in [False, True]:
                        group = "Best-pair search" if pair_mode else "Query search"
                        optimum = max(jaccard(DATABASE[i], DATABASE[j]) for i in range(4) for j in range(i+1,4)) if pair_mode else max(containment(a, QUERY) for a in DATABASE)
                        for trial in range(args.search_trials):
                            run(group, "dna_database", trial,
                                f"SearchDatabase({literal(DATABASE)},{literal(QUERY)},{literal(pair_mode)},5,0.12,6,4)", optimum)
                else:
                    for case in instances["hash_cases"]:
                        keys, aa, bb, modulus = [case[k] for k in ["keys", "multipliers", "offsets", "modulus"]]
                        expected = [min((a*x+b) % modulus for x in keys) for a,b in zip(aa,bb)]
                        for trial in range(args.trials):
                            run("Multi-hash MinHash", case["name"], trial,
                                f"MultiHashMinima({literal(keys)},{literal(aa)},{literal(bb)},{modulus},2)", expected)
                metadata.update(records=record_number, simulator_wall_seconds=time.perf_counter()-started,
                                utc_finished=datetime.now(timezone.utc).isoformat())
            (out / f"metadata_{task}.json").write_text(json.dumps(metadata, indent=2)+"\n")
            print(f"Completed {task}: {record_number} trials in {metadata['simulator_wall_seconds']:.2f}s", flush=True)
    if args.task in ["all", "summarize"]:
        summarize(out)


def summarize(out):
    records = []
    for path in sorted(out.glob("trials_*.jsonl")):
        records.extend(json.loads(line) for line in path.read_text().splitlines() if line.strip())
    if not records:
        raise SystemExit("No trial files found")
    groups = defaultdict(list)
    cases = defaultdict(list)
    for row in records:
        groups[row["group"]].append(row)
        cases[(row["group"], row["case"])].append(row)
    def stats(rows):
        success = sum(x["success"] for x in rows)
        return dict(trials=len(rows), successes=success, success_rate=success/len(rows),
                    wilson_95=wilson(success,len(rows)), mean_error=statistics.mean(x["error"] for x in rows),
                    max_error=max(x["error"] for x in rows),
                    mean_logical_work=statistics.mean(x["logical_work"] for x in rows),
                    work_unit=rows[0]["work_unit"], metric=rows[0]["metric"],
                    total_simulator_seconds=sum(x["simulator_seconds"] for x in rows))
    order = ["Containment QAE", "Jaccard via QAE", "Query search", "Best-pair search", "Multi-hash MinHash"]
    summary = {k:stats(groups[k]) for k in order if k in groups}
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    with (out/"summary.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=["task"]+list(next(iter(summary.values())).keys()))
        writer.writeheader()
        for key,value in summary.items(): writer.writerow(dict(task=key,**value))
    case_summary = {group+" / "+case:stats(rows) for (group,case),rows in cases.items()}
    (out/"summary_by_instance.json").write_text(json.dumps(case_summary,indent=2)+"\n")
    settings = {
        "Containment QAE": r"$|A|\in\{3,4\}$; 7 instances",
        "Jaccard via QAE": r"$|A|=|B|=4$",
        "Query search": r"$R=4$; $|A_i|=4$",
        "Best-pair search": r"$R=4$; 6 legal pairs",
        "Multi-hash MinHash": r"$n\in\{4,8,16\}$; $s=3$",
    }
    def fmt(x): return f"{x:.4f}"
    lines=[r"% Auto-generated from actual trials_*.jsonl; do not edit the numbers.",
           r"\begin{table}[t]",r"\centering",r"\small",
           r"\begin{tabular}{llllr}",r"\hline",
           r"Task & Instance size & Metric & Mean & Successes \\",r"\hline"]
    for task, st in summary.items():
        metric={"absolute error":"Absolute error", "score regret":"Score regret", "incorrect coordinates":"Wrong coordinates"}[st["metric"]]
        lines.append("{} & {} & {} & {:.4f} & {}/{} {}".format(task, settings[task], metric, st["mean_error"], st["successes"], st["trials"], chr(92)*2))
    lines += [r"\hline",r"\end{tabular}",
        r"\caption{Functional correctness in noiseless Q\# simulations. Success requires absolute error or score regret at most $0.1$, or an exactly correct complete MinHash signature. Search uses the finite-budget parameters described in the text. All trials, including failures, are included.}",
        r"\label{tab:qsharp-functional}",r"\end{table}"]
    (out/"paper_table.tex").write_text("\n".join(lines)+"\n")
    markdown=["| Task | Mean error / regret | Successes | 95% Wilson interval |", "|---|---:|---:|---:|"]
    for key, st in summary.items():
        lo,hi=st['wilson_95']
        markdown.append(f"| {key} | {st['mean_error']:.6f} | {st['successes']}/{st['trials']} | [{lo:.3f}, {hi:.3f}] |")
    (out/"paper_table.md").write_text("\n".join(markdown)+"\n")
    print("\n".join(markdown),flush=True)


if __name__ == "__main__":
    main()
