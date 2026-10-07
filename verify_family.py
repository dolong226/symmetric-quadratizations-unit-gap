"""Exact check of the subset-sum construction for exact 1- and 2-out-of-n functions.

For m >= 3 let w = (3, 5, 7, 14, 28, ..., 14 * 2^(m-4)) and N_m = 14 * 2^(m-3) - 2, and let

    g(x, y) = (s - B(y)) (s - B(y) - 1) / 2,    B(y) = a + sum_j w_j y_j,    s = |x|.

Claim (a): with a = -1, min_y g equals the exact 1-out-of-n function for every 2 <= n <= N_m.
Claim (b): with a =  0, min_y g equals the exact 2-out-of-n function for every 3 <= n <= N_m + 1.

The polynomial depends on x only through s, so it suffices to evaluate all weights 0..n.
The script also confirms that the ranges cannot be extended with these weights.
Integer arithmetic only; no dependencies. Run:  python verify_family.py
"""

from itertools import product


def weights(m):
    return [3, 5, 7][:m] + [14 * 2 ** (j - 4) for j in range(4, m + 1)]


def minimum_profile(a, w, n):
    starts = {a + sum(wj * yj for wj, yj in zip(w, y)) for y in product((0, 1), repeat=len(w))}
    return [min((s - b) * (s - b - 1) // 2 for b in starts) for s in range(n + 1)]


def largest_valid_n(a, w, k, n_first, n_last):
    """Largest n such that the construction is exact for every size n_first..n."""
    best = None
    for n in range(n_first, n_last + 1):
        if minimum_profile(a, w, n) != [int(s == k) for s in range(n + 1)]:
            break
        best = n
    return best


def main():
    for m in range(3, 9):
        w = weights(m)
        n_m = 14 * 2 ** (m - 3) - 2
        cap = 2 ** (m + 1)
        got_one = largest_valid_n(-1, w, 1, 2, cap)
        got_two = largest_valid_n(0, w, 2, 3, cap)
        assert got_one == n_m, (m, got_one, n_m)
        assert got_two == n_m + 1, (m, got_two, n_m + 1)
        print(
            f"PASS m={m} weights={w}: exact 1-of-n for 2<=n<={got_one}, "
            f"exact 2-of-n for 3<=n<={got_two}; ranges are maximal for these weights"
        )


if __name__ == "__main__":
    main()
