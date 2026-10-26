"""quantum_optimizer.py - QUBO with optional per-pathway diversity cap."""
from __future__ import annotations
import argparse, json
from itertools import combinations
from pathlib import Path
import numpy as np, pandas as pd

PATHWAYS = {
    "mitochondrial_dysfunction": {"NDUFA1","NDUFB2","SDHA","SDHB","COX5A",
                                  "ATP5A1","ATP5B","UQCRC1","UQCRC2","CYCS"},
    "neuroinflammation": {"IL6","TNF","IL1B","CXCL8","CCL2","IL10","IL4",
                          "TGFB1","NLRP3","TLR4"},
    "proteostasis": {"HSPA1A","HSPA8","PSMD1","UBB","SQSTM1","HSPB1",
                     "DNAJB1","BAG3","MAP1LC3B","ATG5"},
    "synaptic_function": {"RAB3A","SNAP25","SYN1","CALM3","CDC42",
                          "STX1A","VAMP2","SNAP29","RAB5A","DNM1"},
}


def _pathways_of(name):
    return {p for p, m in PATHWAYS.items() if name in m}


def _respects_cap(combo, names, cap):
    counts = {}
    for i in combo:
        for p in _pathways_of(names[i]):
            counts[p] = counts.get(p, 0) + 1
    return all(c <= cap for c in counts.values())


def solve_exact(rel, ov, k, names, cap=None):
    n = len(rel)
    best_obj, best_combo = -np.inf, None
    for combo in combinations(range(n), k):
        if cap is not None and not _respects_cap(combo, names, cap):
            continue
        idx = np.array(combo)
        obj = float(rel[idx].sum()) + float(ov[np.ix_(idx, idx)].sum()) / 2.0
        if obj > best_obj:
            best_obj, best_combo = obj, combo
    if best_combo is None:
        raise SystemExit(
            f"[err] no feasible solution for k={k}, cap={cap}. Relax."
        )
    return best_obj, list(best_combo)


def solve_qaoa(rel, ov, k, reps, maxiter, names, cap=None):
    from qiskit_optimization import QuadraticProgram
    from qiskit_optimization.algorithms import MinimumEigenOptimizer
    from qiskit_optimization.converters import LinearInequalityToPenalty
    from qiskit_algorithms import QAOA
    from qiskit_algorithms.optimizers import COBYLA
    from qiskit.primitives import StatevectorSampler

    n = len(rel)
    if n > 12:
        raise SystemExit(f"[err] QAOA needs n<=12 to fit statevector. Got n={n}.")

    qp = QuadraticProgram("shared_module")
    for i in range(n):
        qp.binary_var(name=f"x_{i}")
    qp.minimize(
        linear={f"x_{i}": -float(rel[i]) for i in range(n)},
        quadratic={(f"x_{i}", f"x_{j}"): -float(ov[i, j])
                   for i in range(n) for j in range(i + 1, n) if ov[i, j] != 0},
    )
    qp.linear_constraint(
        linear={f"x_{i}": 1 for i in range(n)}, sense="==", rhs=k, name="select_k"
    )
    if cap is not None:
        for pathway, members in PATHWAYS.items():
            indices = [i for i, nm in enumerate(names) if nm in members]
            if len(indices) > cap:
                qp.linear_constraint(
                    linear={f"x_{i}": 1 for i in indices},
                    sense="<=", rhs=cap, name=f"cap_{pathway}",
                )
        # Convert <= to penalty so MinimumEigenOptimizer does not add slack qubits
        qp = LinearInequalityToPenalty(penalty=100.0).convert(qp)

    qaoa = QAOA(sampler=StatevectorSampler(),
                optimizer=COBYLA(maxiter=maxiter), reps=reps)
    result = MinimumEigenOptimizer(qaoa).solve(qp)
    x = np.asarray(result.x).astype(int)
    return -float(result.fval), list(np.where(x == 1)[0])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--features", default="data/processed/features.csv")
    p.add_argument("--overlap", default="data/processed/overlap.npy")
    p.add_argument("--k", type=int, default=4)
    p.add_argument("--solver", choices=["exact", "qaoa"], default="exact")
    p.add_argument("--reps", type=int, default=1)
    p.add_argument("--maxiter", type=int, default=30)
    p.add_argument("--max-per-pathway", type=int, default=None)
    p.add_argument("--out", default="results/quantum_selection.json")
    args = p.parse_args()

    names = json.loads(Path(args.overlap).with_suffix(".names.json").read_text())
    rel = pd.read_csv(args.features, index_col=0).loc[names]["shared_score"].to_numpy()
    ov = np.load(args.overlap)

    print(f"[info] solver={args.solver}  n={len(names)}  k={args.k}  cap={args.max_per_pathway}")
    if args.solver == "exact":
        obj, idx = solve_exact(rel, ov, args.k, names, args.max_per_pathway)
    else:
        obj, idx = solve_qaoa(rel, ov, args.k, args.reps, args.maxiter,
                              names, args.max_per_pathway)

    selected = [names[i] for i in idx]
    payload = {"solver": args.solver, "k": args.k, "reps": args.reps,
               "max_per_pathway": args.max_per_pathway,
               "objective_value": float(obj), "selected_features": selected}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(payload, indent=2))
    print(f"[ok ] objective = {obj:.6f}")
    print(f"[ok ] selected {len(selected)} -> {args.out}")
    print("      " + ", ".join(selected))


if __name__ == "__main__":
    main()
