"""Bounded search for x-symmetric quadratizations of exact k-out-of-n functions.

An x-symmetric quadratic polynomial with m auxiliary variables has the normal form

    g_y(s) = C s^2/2 + (beta0 + sum_j E_j y_j) s + A + sum_j D_j y_j + sum_{j<l} F_jl y_j y_l,

where s is the Hamming weight of x. The script decides, by mixed-integer linear
programming, whether real coefficients with absolute value at most BOUND exist such that
min_y g_y(s) equals 1 at s = k and 0 at every other weight 0..n.

Binary variables z[s, y] mark a branch y that attains the minimum at weight s:
    g_y(s) >= r(s)                      for all s, y
    g_y(s) <= r(s) + M (1 - z[s, y])    for all s, y
    sum_y z[s, y] >= 1                  for all s

An INFEASIBLE answer is a statement about the bounded box only. It is numerical evidence,
not a proof, and it is reported as such in the paper.

Requires numpy and scipy >= 1.9 (HiGHS).

Usage:  python search_normal_form.py n,k,m [n,k,m ...]
        python search_normal_form.py            (the cases n = 14, 15 quoted in the paper)
"""

import itertools
import sys

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

BOUND = 3000.0
TIME_LIMIT = 900

STATUS = {
    0: "FEASIBLE",
    1: "LIMIT REACHED (unknown)",
    2: "INFEASIBLE",
    3: "UNBOUNDED",
    4: "NUMERICAL TROUBLE",
}


def search(n, k, m, bound=BOUND, time_limit=TIME_LIMIT):
    """Return (status string, coefficient dict or None)."""
    ys = list(itertools.product((0, 1), repeat=m))
    pairs = list(itertools.combinations(range(m), 2))
    names = (
        ["C", "beta0"]
        + [f"E{j}" for j in range(m)]
        + ["A"]
        + [f"D{j}" for j in range(m)]
        + [f"F{j}{l}" for j, l in pairs]
    )
    index = {name: i for i, name in enumerate(names)}
    n_real = len(names)
    n_bin = (n + 1) * len(ys)
    n_var = n_real + n_bin

    def z(s, yi):
        return n_real + s * len(ys) + yi

    r = [int(s == k) for s in range(n + 1)]
    big_m = bound * (2 + n + n * n / 2 + m * n + m + len(pairs)) + 10

    rows, lower, upper = [], [], []
    for s in range(n + 1):
        for yi, y in enumerate(ys):
            row = np.zeros(n_var)
            row[index["C"]] = s * s / 2
            row[index["beta0"]] = s
            row[index["A"]] = 1
            for j in range(m):
                row[index[f"E{j}"]] = s * y[j]
                row[index[f"D{j}"]] = y[j]
            for j, l in pairs:
                row[index[f"F{j}{l}"]] = y[j] * y[l]
            rows.append(row.copy())
            lower.append(r[s])
            upper.append(np.inf)
            tight = row.copy()
            tight[z(s, yi)] = big_m
            rows.append(tight)
            lower.append(-np.inf)
            upper.append(r[s] + big_m)
        cover = np.zeros(n_var)
        for yi in range(len(ys)):
            cover[z(s, yi)] = 1
        rows.append(cover)
        lower.append(1)
        upper.append(np.inf)

    result = milp(
        np.zeros(n_var),
        constraints=LinearConstraint(np.array(rows), np.array(lower), np.array(upper)),
        integrality=np.concatenate([np.zeros(n_real), np.ones(n_bin)]),
        bounds=Bounds(
            np.concatenate([np.full(n_real, -bound), np.zeros(n_bin)]),
            np.concatenate([np.full(n_real, bound), np.ones(n_bin)]),
        ),
        options={"time_limit": time_limit},
    )
    status = STATUS.get(result.status, str(result.status))
    if result.status != 0:
        return status, None
    return status, {name: round(float(v), 6) for name, v in zip(names, result.x[:n_real])}


def gap_values(n):
    """Values 0 < k <= n/2 for which the known upper bound exceeds the lower bound."""
    m = (n - 1).bit_length() - 1  # ceil(log2 n) - 1
    return m, [k for k in range(1, n // 2 + 1) if n - k > 2**m]


def main(argv):
    if argv:
        cases = [tuple(map(int, arg.split(","))) for arg in argv]
    else:
        cases = []
        for n in (14, 15):
            m, ks = gap_values(n)
            cases += [(n, k, m) for k in ks]
    for n, k, m in cases:
        status, coefficients = search(n, k, m)
        print(f"n={n} k={k} m={m} |coefficients|<={BOUND:g}: {status}", flush=True)
        if coefficients:
            print(f"    {coefficients}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
