# Q# functional simulations for quantum similarity and sketch construction

Tested with Python 3.12 and Microsoft QDK 1.33.1 on the local noiseless simulator. 
All seeds are fixed (BASE_SEED = 20261021).

## Files

- `src/Core.qs`: coherent indexed key loading, exact membership, uniform index
  preparation, Grover iterate, phase-estimation QAE, containment/Jaccard.
- `src/Search.qs`: QAE threshold marking, amplitude amplification, query search,
  and best-pair search.
- `src/MinHash.qs`: reversible affine-hash arithmetic and quantum minimum finding
  for multiple prescribed hashes.
- `run_experiments.py`: builds the inputs, runs seeded Q# trials, computes
  classical ground truth, and writes the summary table.
- `run_search_parallel.py`: runs the search trials in parallel worker processes.
- `results/`: the run reported in the paper (768 trials).
  - `trials_*.jsonl`: every trial, including failures.
  - `metadata_*.json`: runtime versions, seeds, parameters, source hashes.
  - `instances.json`: all input sets and hash seeds.
  - `summary.json`, `summary.csv`, `summary_by_instance.json`, `paper_table.md`.

## Run

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

python run_experiments.py --task estimators  
python run_experiments.py --task minhash      
python run_search_parallel.py --workers 4     
python run_experiments.py --task summarize    
```


## Experimental setup

- **Containment/Jaccard:** seven set pairs (source size 3 or 4, overlap 0..4,
  empty target, non-power-of-two source) plus one Jaccard pair. QAE uses 6 phase
  qubits and the median of 3 runs.
- **Search:** distinct 2-mers of `AACCA`, `ACGTA`, `ACGTT`, `CCGGC`
  (A=0,C=1,G=2,T=3, base-4 encoding); query `ACGTT`. QAE uses 5 phase qubits
  and one coherent QAE register. Margin is 0.12, with 6 attempts per update and
  at most 4 update rounds. Search starts from a suboptimal candidate.
- **MinHash:** n ∈ {4, 8, 16}, hashes h(x) = (a·x + b) mod 32 with
  (a,b) ∈ {(3,1), (5,2), (7,4)}, and 2 minimum-finding restarts per hash.
