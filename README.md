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
