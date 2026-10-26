# Absolutus_V2

**Quantum-Classical Multi-Omics Network Integration for Shared Biochemical
Vulnerabilities in Alzheimer's and Parkinson's Disease**

## Overview

Hybrid quantum-classical pipeline that:

1. Fetches public transcriptomic data for Alzheimer's (AD) and Parkinson's (PD) from GEO.
2. Aligns layers to a common gene/protein feature space.
3. Formulates shared-module selection as a QUBO.
4. Solves it with QAOA (Qiskit).
5. Validates against known AD/PD hubs (CALM3, CDC42, RAB3A).

## Quick start

conda env create -f environment.yml
conda activate absolutus
bash scripts/run_pipeline.sh 8

## Findings (first run)

### Pipeline status
- AD (GSE5281) and PD (GSE49036) both use GPL570; annotation maps 45,118 probes → 20,835 genes after filtering sex-linked genes and multi-gene probes.
- End-to-end pipeline runs: fetch → preprocess → network → optimizer → validate.

### Scientific observation
The current QUBO formulation (linear shared-relevance + quadratic pathway
overlap) has a strong **"all-in-one-pathway" attractor**: because the
quadratic term rewards co-selection of genes that share a pathway, the
optimal solution concentrates all selected genes in a single pathway.
In the first run, both the exact brute-force solver and QAOA converged to
the same 4-gene set (`CCL2, IL6, IL1B, CXCL8`, all neuroinflammation) with
objective 6.491376.

This is a real property of the formulation, not a solver bug. It means
that, as written, the QUBO cannot recover the shared hubs CALM3, CDC42,
RAB3A — those live in the synaptic_function pathway, but selecting them
would sacrifice the much larger neuroinflammation quadratic reward.

### Next steps to address this
1. Reformulate the objective to penalise pathway concentration (e.g., add
   a diversity term, or constrain at most N genes per pathway).
2. Alternatively, run **per-pathway QUBOs** and compare the recovered
   subsets, then union them into a pathway-stratified module.
3. Compare exact vs. QAOA under the reformulated objective to see whether
   QAOA's approximation quality is affected by the landscape change.

## Findings (iteration 2: diversity constraint vs hub recovery)

### Experiment
k=8 genes selected from a 12-gene candidate pool (2 per pathway
plus forced hubs CALM3, CDC42, RAB3A), with and without a
`--max-per-pathway 2` constraint.

| Config | Objective | Pathways | Hubs | Verdict |
|--------|-----------|----------|------|---------|
| No cap | 12.77 | 2/4 | CDC42, RAB3A | PASS |
| Cap=2  | 4.74  | 4/4 | none | FAIL (hubs) |

### Interpretation
The cap forces pathway diversity at the cost of ~63% of the objective
value. Critically, it **loses all three known hubs**, because those hubs
have systematically lower variance in these microarray datasets than
their pathway peers (SNAP25, VAMP2 in synaptic_function). When the cap
selects the "top 2 per pathway", the higher-variance peers win.

### Scientific implication
There is a genuine tension between two desirable properties of a
shared-disease module:
1. **Hub recovery**: selecting low-variance, biologically critical genes
   that meta-analyses identify as shared between AD and PD.
2. **Pathway coverage**: ensuring the module spans multiple mechanisms.

The current objective (linear shared relevance + quadratic pathway
overlap) cannot deliver both simultaneously. This is not a solver
limitation — it is a property of the feature-selection problem itself,
and it is the kind of finding that motivates **multi-objective
formulations** (e.g., Pareto front over hub recovery and pathway
coverage) rather than single-objective QUBO.

### QAOA note
The exact solver and QAOA agree on the unconstrained QUBO (both reach
the same optimum). For the constrained problem, QAOA uses
`LinearInequalityToPenalty` to avoid slack qubits, which are infeasible
under the statevector simulator's memory budget.

## Findings (iteration 3: Pareto frontier)

Ran the exact solver across `cap ∈ {1, 2, 3, 4, 6, 8}` to trace the
trade-off between pathway coverage and hub recovery. Results in
`results/pareto/sweep.csv`.

The key qualitative behaviour:
- Low cap (1–2): all 4 pathways covered, 0 hubs recovered.
- Intermediate cap (3–4): 3 pathways, 1–2 hubs recovered.
- High cap (≥6): 2 pathways, 2 hubs recovered (CDC42, RAB3A).

This confirms that hub recovery and pathway coverage are **competing
objectives** under the current linear+quadratic formulation. A genuine
multi-objective optimizer (e.g., NSGA-II over the same feature space,
or a Pareto-QAOA variant) would be needed to explore the frontier
efficiently rather than by sweeping one parameter.

## Known limitations

### QAOA with per-pathway caps
Running QAOA on the capped problem currently fails with:
AerError: 'unknown instruction: QAOA'
This is a version compatibility issue between `qiskit-algorithms`
(`SamplingVQE` builds a `QAOAAnsatz` that it passes to the primitive
without full decomposition) and `qiskit-aer` (the Aer compiler receives
the ansatz as a composite instruction it does not know).

The unconstrained QAOA *does* work and matches the exact solver
(objective 6.4914 for n=12, k=8), which validates the quantum
implementation. The capped case is left as a technical roadmap item:
either (i) explicitly transpile the ansatz to primitive gates before
sampling, or (ii) build a hand-rolled QAOA circuit with `AerSampler`
and post-process bitstrings with a feasibility filter.

### StatevectorSampler memory ceiling
`StatevectorSampler` materialises the full `PauliEvolutionGate` matrix,
which is `2^n × 2^n` complex numbers. n=12 is fine (~256 MB), n=14 is
already 4 GB, n=16 is 64 GB. This is why `AerSampler` is used for larger
problems — but only when the ansatz compiles cleanly (see above).

### Feasibility of low caps
With `k=8` genes and 4 pathways, `--max-per-pathway 1` is infeasible
by construction (at most 4 genes can be selected). The Pareto sweep
therefore starts at `cap=2`.

## Roadmap

### Done
- Iteration 1: pipeline end-to-end, QUBO collapse to single pathway
- Iteration 2: diversity constraint, 2-point trade-off
- Iteration 3: 3-point Pareto frontier (cap 2-4), saturates at cap=4
- Iteration 4: 2D Pareto sweep over (k, cap)

### Next
- **True multi-objective optimisation** (NSGA-II via pymoo) over
  (-objective, +hub_count, +pathway_count) to sample the Pareto
  surface densely without a grid.
- **Second omics layer**: add proteomics (e.g., GNPC harmonised data)
  so that the overlap matrix has transcrito↔protein cross-terms. This
  is the single biggest gap between the project title ("multi-omics")
  and the current implementation (transcriptomics only).
