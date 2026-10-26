"""quantum_optimizer.py - QUBO solver (exact or QAOA)."""
from __future__ import annotations
import argparse, json
from itertools import combinations
from pathlib import Path
import numpy as np, pandas as pd


def solve_exact(rel, ov, k):
    n = len(rel)
    best_obj, best_combo = -np.inf, None
    for combo in combinations(range(n), k):
        idx = np.array(combo)
        obj = float(rel[idx].sum()) + float(ov[np.ix_(idx, idx)].sum()) / 2.0
        if obj > best_obj:
            best_obj, best_combo = obj, combo
    return best_obj, list(best_combo)


def solve_qaoa(rel, ov, k, reps, maxiter):
    from qiskit_optimization import QuadraticProgram
    from qiskit_optimization.algorithms import MinimumEigenOptimizer
    from qiskit_algorithms import QAOA
    from qiskit_algorithms.optimizers import COBYLA
    from qiskit.primitives import StatevectorSampler

    n = len(rel)
    if n > 12:
        raise SystemExit(f"[err] QAOA needs n<=12. Got n={n}.")

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
    qaoa = QAOA(sampler=StatevectorSampler(),
                optimizer=COBYLA(maxiter=maxiter), reps=reps)
    result = MinimumEigenOptimizer(qaoa).solve(qp)
    x = np.asarray(result.x).astype(int)
    # QAOA minimises the negated objective; flip back to maximisation value
    return -float(result.fval), list(np.where(x == 1)[0])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--features", default="data/processed/features.csv")
    p.add_argument("--overlap", default="data/processed/overlap.npy")
    p.add_argument("--k", type=int, default=4)
    p.add_argument("--solver", choices=["exact", "qaoa"], default="exact")
    p.add_argument("--reps", type=int, default=1)
    p.add_argument("--maxiter", type=int, default=30)
    p.add_argument("--out", default="results/quantum_selection.json")
    args = p.parse_args()

    names = json.loads(Path(args.overlap).with_suffix(".names.json").read_text())
    rel = pd.read_csv(args.features, index_col=0).loc[names]["shared_score"].to_numpy()
    ov = np.load(args.overlap)

    print(f"[info] solver={args.solver}  n={len(names)}  k={args.k}")
    if args.solver == "exact":
        obj, idx = solve_exact(rel, ov, args.k)
    else:
        obj, idx = solve_qaoa(rel, ov, args.k, args.reps, args.maxiter)

    selected = [names[i] for i in idx]
    payload = {"solver": args.solver, "k": args.k, "reps": args.reps,
               "objective_value": float(obj), "selected_features": selected}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(payload, indent=2))
    print(f"[ok ] objective = {obj:.6f}")
    print(f"[ok ] selected {len(selected)} -> {args.out}")
    print("      " + ", ".join(selected))


if __name__ == "__main__":
    main()
