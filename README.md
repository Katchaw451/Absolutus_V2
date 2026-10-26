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

Ran the exact solver across `cap ∈ {2, 3, 4, 6, 8}` to trace the trade-off between pathway coverage and hub recovery (`cap=1` is infeasible for `k=8` with 4 pathways: at most 4 genes can be selected).
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

## Findings (iteration 4: 2D Pareto surface over k and cap)

Sweeping k ∈ {4, 6, 8, 10} × cap ∈ {2, 3, 4} with the exact solver
produces a clean 2D structure:

| k \ cap | 2              | 3              | 4              |
|---------|----------------|----------------|----------------|
| 4       | 2 pw, 0 hubs   | 2 pw, 0 hubs   | 1 pw, 0 hubs   |
| 6       | 3 pw, 0 hubs   | 2 pw, 1 hub    | 2 pw, 0 hubs   |
| 8       | 4 pw, 0 hubs   | 3 pw, 1 hub    | 2 pw, 2 hubs   |
| 10      | (infeasible)   | 4 pw, 1 hub    | 3 pw, 2 hubs   |

Three regularities emerge:

1. **k and cap act almost orthogonally.** Holding cap fixed and
   increasing k increases pathway coverage (rows). Holding k fixed and
   increasing cap increases hub recovery (columns). This is a
   consequence of the QUBO design: k controls how many features are
   drawn, cap controls how they are distributed across pathways, and
   hub recovery depends on whether the two hubs (CDC42, RAB3A) survive
   the top-cap cut within synaptic_function.

2. **k=10 strictly dominates k=8 at cap=4.** Same hub count (2), more
   pathways (3 vs 2), higher objective (13.92 vs 12.77). This is a
   genuine Pareto improvement, revealing that the k=8 Pareto frontier
   was a projection artefact of fixing k. Fixing k can hide better
   solutions; the correct object is the k-agnostic frontier.

3. **CALM3 is never recovered at any (k, cap).** The calmodulin gene
   has systematically lower variance in GSE5281/GSE49036 than its
   synaptic peers CDC42, RAB3A, SNAP25, VAMP2. Biologically this is
   expected — CALM3 is one of the most conserved eukaryotic proteins
   and has near-housekeeping expression. Variance-based feature
   selection therefore *systematically excludes* highly conserved
   disease hubs. Any method that ranks features purely by variance
   will miss this class of genes. This motivates the next iteration:
   use network connectivity (STRING, BioGRID) as an alternative or
   complementary ranking criterion.

### Infeasibility
k=10 with cap=2 is infeasible by construction: 4 pathways × cap 2 = at
most 8 genes can be selected. The sweep skips it silently.

## Findings (iteration 5: network degree vs variance)

With CALM3 now in the candidate pool (recovered via multi-gene probe
resolution), we compared the two score fields across caps:

| cap | shared_score (variance only) | priority (variance + STRING degree) |
|-----|------------------------------|--------------------------------------|
| 2   | 4 pathways, **0 hubs**       | 4 pathways, **2 hubs** (CALM3, CDC42) |
| 3   | 3 pathways, 1 hub            | 3 pathways, 2 hubs                   |
| 4   | 3 pathways, 2 hubs (CDC42, RAB3A) | 3 pathways, 2 hubs (CALM3, CDC42) |
| 6   | 2 pathways, **3 hubs**       | 2 pathways, **3 hubs**               |
| 8   | 2 pathways, 3 hubs           | 2 pathways, 3 hubs                   |

Three results:

1. **Degree is decisive under a tight cap.** At cap=2, variance-only
   recovers zero hubs; degree-enhanced recovers two (CALM3, CDC42) with
   identical pathway coverage (4/4). No trade-off — pure dominance.

2. **At cap>=6 the two scores converge** (both recover all three hubs),
   because the constraint stops binding: the 12-feature pool cannot
   violate a cap of 6 within any single pathway. The marginal value of
   the degree signal is therefore concentrated in the regime where the
   problem is genuinely constrained.

3. **CALM3 and RAB3A are opposite hub sub-types.** Variance-only prefers
   RAB3A (moderate variance, STRING degree 53); degree-enhanced prefers
   CALM3 (low variance, STRING degree 348). No single score recovers
   both simultaneously at low cap. This reframes the "CALM3 problem" as
   a finding: AD/PD hubs are not a homogeneous class, and any single
   ranking signal will bias toward one sub-type.

### Implication
The natural next experiment is a **two-parameter sweep** over
`--degree-weight` in {0, 0.1, 0.25, 0.5, 1.0, 2.0}, tracing the
Pareto frontier between CALM3 recovery (high weight) and RAB3A
recovery (low weight). This makes the hub-heterogeneity claim
quantitative rather than anecdotal.

## License

MIT License — see [LICENSE](LICENSE).

## Limitations

The current implementation is a **well-documented toy problem**, not a
scientific result. The following limitations are structural and were
identified in external review; they are listed here so that the gap
between this repository and a publishable study is explicit.

### 1. Circular validation
The candidate pool has 12 genes. The known hubs (CALM3, CDC42, RAB3A)
are forced into the pool via `--include-hubs`, and the four pathways
are hand-curated. The `PASS`/`FAIL` outcome therefore measures whether
the solver selects among 12 pre-chosen candidates, not whether the
method discovers anything. A valid design requires:
- the full gene universe (~20,000 genes), or at least several hundred;
- the hub list fixed **before** the experiment, independently of the
  pool construction;
- a permutation null (repeat the pipeline with random gene sets of
  equal size) to calibrate what "recovered" means.

### 2. The quantum component is decorative at this scale
`C(12, 8) = 495` combinations. The exact brute-force solver is
instantaneous. QAOA matching the exact solver validates the
implementation, not any quantum advantage. The constrained case (the
one that matters) fails at runtime. **There is no quantum result in
this repository — the result is purely classical.** To make QAOA
meaningful, the problem must scale beyond brute force (`n >= 25`), and
must be compared against strong classical baselines (ILP, simulated
annealing, DOMINO, DIAMOnD, heinz), not against an exact solver on a
12-variable problem.

### 3. STRING degree is study-biased
CALM3's degree of 348 reflects, in part, how promiscuous and how
well-studied calmodulin is. Recovering degree-defined hubs using a
degree-based score is partially circular. The claim that "CALM3 and
RAB3A are opposite hub sub-types" rests on two genes and is a
hypothesis, not a result. A valid test would use many more genes and
compare against degree-matched controls.

### 4. Variance is not disease relevance
The score uses raw expression variance, not case/control effect sizes.
GSE5281 (multiple brain regions) and GSE49036 (substantia nigra) are
different tissues, and the synaptic-gene signal is likely dominated by
neuronal loss and cell-composition shifts. Without deconvolution or
covariate adjustment, "shared between AD and PD" may mean only
"neurons die in both".

### 5. The QUBO formulation rewards the wrong thing
The quadratic pathway-overlap term rewards concentration (all genes in
one pathway), which is the opposite of the desired diversity. This was
discovered empirically (iteration 1) and patched with a hard cap, but
the patch does not fix the formulation. The natural formulation is
**densest-k-subgraph on a PPI network**: linear weight = evidence,
quadratic weight = network cohesion, plus a redundancy penalty. That is
where the QUBO literature on disease modules provides both a principled
objective and classical baselines.

### 6. "Multi-omics" is transcriptomics only
The repository title over-promises. Only transcriptomic data is used.
Adding proteomics (e.g., GNPC harmonised data) requires an access
request and should be started early.

### 7. Artefacts of a small pool
The "k and cap are nearly orthogonal" regularity and the
"CALM3 ↔ RAB3A trade-off" are artefacts of a 12-gene pool across 4
pathways. The trade-off is manufactured by `cap=4`: at `cap=6` all
three hubs coexist. The `--degree-weight` sweep is a step function,
not a smooth frontier, and is driven by degree-ranked substitution
(RAB3A → CALM3, HSPB1 → COX5A). Comparing objective values across
different `k` is also invalid, because the positive linear terms sum.

### 8. QAOA technical issues (separate from the formulation)
Two distinct failures were observed:
- `unknown instruction: QAOA` — the `AerSampler` receives the
  `QAOAAnsatz` as a composite instruction. Likely fixable by
  transpiling the ansatz to basis gates before sampling.
- 64 GB allocation at n=16 — the `StatevectorSampler` materialises
  `PauliEvolutionGate.to_matrix()`, which is `2^n × 2^n` complex
  numbers. This is a different problem and is not fixed by transpiling.

### What this repository is good for
- A reproducible pipeline (`environment.yml`, CI, tests, license).
- An honest record of negative results (the all-in-one-pathway
  attractor, the QAOA-with-constraints failure, the k=10/cap=2
  infeasibility).
- A documented baseline that a future `Absolutus_V3` can compare
  against, with the five structural problems above explicitly
  addressed in the new design.

## Roadmap

### Done
- Iteration 1: pipeline end-to-end, QUBO collapse to single pathway
- Iteration 2: diversity constraint, 2-point trade-off
- Iteration 3: 3-point Pareto frontier (cap 2-4), saturates at cap=4
- Iteration 4: 2D Pareto sweep over (k, cap)
- Iteration 5: STRING degree recovers CALM3 under tight cap
- Iteration 6: degree-weight sweep (step function, no sweet spot)

### Planned (Absolutus_V3)
- Full gene universe (~20,000), not a 12-gene pool
- Hub list fixed a priori, independently of pool construction
- Permutation null to calibrate hub recovery
- Case/control effect sizes with cell-type deconvolution
- Densest-k-subgraph formulation on STRING/BioGRID PPI
- Strong classical baselines (DOMINO, DIAMOnD, heinz, ILP)
- QAOA only at n >= 25, compared with SA/ILP
- Multi-omics (GNPC proteomics) as a distinct milestone

## Limitations

The current implementation is a **well-documented toy problem**, not a
scientific result. The following limitations are structural and were
identified in external review; they are listed here so that the gap
between this repository and a publishable study is explicit.

### 1. Circular validation
The candidate pool has 12 genes. The known hubs (CALM3, CDC42, RAB3A)
are forced into the pool via `--include-hubs`, and the four pathways
are hand-curated. The `PASS`/`FAIL` outcome therefore measures whether
the solver selects among 12 pre-chosen candidates, not whether the
method discovers anything. A valid design requires the full gene
universe (~20,000 genes) or at least several hundred, the hub list
fixed **before** the experiment, and a permutation null to calibrate
what "recovered" means.

### 2. The quantum component is decorative at this scale
`C(12, 8) = 495` combinations. The exact brute-force solver is
instantaneous. QAOA matching the exact solver validates the
implementation, not any quantum advantage. The constrained case (the
one that matters) fails at runtime. **There is no quantum result in
this repository — the result is purely classical.** To make QAOA
meaningful, the problem must scale beyond brute force (`n >= 25`) and
be compared against strong classical baselines (ILP, simulated
annealing, DOMINO, DIAMOnD, heinz).

### 3. STRING degree is study-biased
CALM3's degree of 348 reflects, in part, how promiscuous and how
well-studied calmodulin is. Recovering degree-defined hubs using a
degree-based score is partially circular. The claim that "CALM3 and
RAB3A are opposite hub sub-types" rests on two genes and is a
hypothesis, not a result.

### 4. Variance is not disease relevance
The score uses raw expression variance, not case/control effect sizes.
GSE5281 (multiple brain regions) and GSE49036 (substantia nigra) are
different tissues, and the synaptic-gene signal is likely dominated by
neuronal loss and cell-composition shifts. Without deconvolution or
covariate adjustment, "shared between AD and PD" may mean only
"neurons die in both".

### 5. The QUBO formulation rewards the wrong thing
The quadratic pathway-overlap term rewards concentration, which is the
opposite of the desired diversity. This was discovered empirically and
patched with a hard cap, but the patch does not fix the formulation.
The natural formulation is **densest-k-subgraph on a PPI network**:
linear weight = evidence, quadratic weight = network cohesion, plus a
redundancy penalty. That is where the QUBO literature on disease
modules provides both a principled objective and classical baselines.

### 6. "Multi-omics" is transcriptomics only
Only transcriptomic data is used. Adding proteomics (e.g., GNPC
harmonised data) requires an access request and should be started
early.

### 7. Artefacts of a small pool
The "k and cap are nearly orthogonal" regularity and the
"CALM3 <-> RAB3A trade-off" are artefacts of a 12-gene pool across 4
pathways. The trade-off is manufactured by `cap=4`: at `cap=6` all
three hubs coexist. The `--degree-weight` sweep is a step function,
not a smooth frontier, driven by degree-ranked substitution
(RAB3A -> CALM3, HSPB1 -> COX5A). Comparing objective values across
different `k` is also invalid, because the positive linear terms sum.

### 8. QAOA technical issues (separate from the formulation)
Two distinct failures were observed:
- `unknown instruction: QAOA` — the `AerSampler` receives the
  `QAOAAnsatz` as a composite instruction. Likely fixable by
  transpiling the ansatz to basis gates before sampling.
- 64 GB allocation at n=16 — the `StatevectorSampler` materialises
  `PauliEvolutionGate.to_matrix()`, which is `2^n x 2^n` complex
  numbers. This is a different problem and is not fixed by
  transpiling.

### What this repository is good for
- A reproducible pipeline (`environment.yml`, CI, tests, license).
- An honest record of negative results.
- A documented baseline that a future `Absolutus_V3` can compare
  against, with the eight structural problems above explicitly
  addressed in the new design.
