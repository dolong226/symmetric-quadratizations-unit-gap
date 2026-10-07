"""Exact finite checks of the equality characterization, the recognition procedure and
the classification of exact k-out-of-n functions at n = 2^p; standard library only.

These checks supplement the proofs. A bounded exhaustive search over integer
generators is an independent oracle for the reconstruction procedure; it is
not a certificate for arbitrary sizes or for nonsymmetric representations.
Run:  python verify_theory.py
"""

from fractions import Fraction
from itertools import combinations, product


def recognize(profile, m):
    """Return the canonical offset, weights and scale, or None at this budget."""
    if m < 1 or not profile or any(v < 0 for v in profile) or not any(profile):
        raise ValueError("Requires a nonzero nonnegative profile and m >= 1")
    zeros = [s for s, value in enumerate(profile) if value == 0]
    if len(zeros) != 2 ** (m + 1):
        raise ValueError("The proposed budget must saturate the zero count")
    if any(v != u + 1 for u, v in zip(zeros[::2], zeros[1::2])):
        return None
    starts = zeros[::2]
    offset = starts[0]
    target = {b - offset for b in starts}
    generated = {0}
    weights = []
    while generated != target:
        weight = min(target - generated)
        shifted = {b + weight for b in generated}
        if shifted & generated or not shifted <= target:
            return None
        weights.append(weight)
        generated |= shifted
    if len(weights) != m:
        return None
    # Deliberately compute the envelope by enumerating all starts for checking.
    envelope = [
        min((s - b) * (s - b - 1) // 2 for b in starts)
        for s in range(len(profile))
    ]
    positive = next(s for s, value in enumerate(profile) if value > 0)
    scale = Fraction(profile[positive]) / envelope[positive]
    if any(value != scale * h for value, h in zip(profile, envelope)):
        return None
    return offset, tuple(weights), scale


def brute_encodings(n, m):
    """Enumerate all positive integer generators bounded by the root domain."""
    result = {}
    for weights in combinations(range(1, n + 1), m):
        sums = sorted(
            sum(w * bit for w, bit in zip(weights, bits))
            for bits in product((0, 1), repeat=m)
        )
        if len(set(sums)) != 2**m:
            continue
        if any(v - u < 2 for u, v in zip(sums, sums[1:])):
            continue
        for offset in range(n - sums[-1]):
            zeros = tuple(s for b in sums for s in (b + offset, b + offset + 1))
            canonical = (offset, weights)
            # A second generator multiset would contradict the rigidity claim.
            assert zeros not in result or result[zeros] == canonical
            result[zeros] = canonical
    return result


def check_recognition():
    zero_sets = scaled_profiles = perturbed_profiles = 0
    for m, sizes in ((1, range(4, 11)), (2, range(8, 15)), (3, (16, 17, 18))):
        for n in sizes:
            oracle = brute_encodings(n, m)
            for zeros in combinations(range(n + 1), 2 ** (m + 1)):
                zero_sets += 1
                paired = all(v == u + 1 for u, v in zip(zeros[::2], zeros[1::2]))
                if paired:
                    starts = zeros[::2]
                    profile = [
                        min((s - b) * (s - b - 1) // 2 for b in starts)
                        for s in range(n + 1)
                    ]
                else:
                    profile = [int(s not in zeros) for s in range(n + 1)]
                answer = recognize(profile, m)
                expected = oracle.get(zeros)
                assert (answer is not None) == (expected is not None), (n, m, zeros)
                if expected is None:
                    continue
                assert answer == (*expected, Fraction(1)), (n, m, zeros, answer)
                scaled = [Fraction(3, 2) * value for value in profile]
                assert recognize(scaled, m) == (*expected, Fraction(3, 2))
                scaled_profiles += 1
                positive = [s for s, value in enumerate(profile) if value > 0]
                if len(positive) >= 2:
                    perturbed = profile.copy()
                    perturbed[positive[0]] *= 2
                    assert recognize(perturbed, m) is None
                    perturbed_profiles += 1
    return zero_sets, scaled_profiles, perturbed_profiles


def check_cardinality():
    recognized_profiles = branch_evaluations = 0
    for p in range(2, 9):
        n = 2**p
        middle = n // 2
        for k in range(n + 1):
            profile = [int(s == k) for s in range(n + 1)]
            answer = recognize(profile, p - 1)
            assert (answer is not None) == (k in (0, middle, n)), (p, k)
            recognized_profiles += 1
            if answer is not None:
                if k == middle:
                    expected = tuple(2**i for i in range(1, p - 1)) + (middle + 1,)
                    offset = 0
                else:
                    expected = tuple(2**i for i in range(1, p))
                    offset = int(k == 0)
                assert answer == (offset, expected, Fraction(1)), (p, k, answer)

    # Independent direct evaluation of the published p-auxiliary construction.
    for p in range(2, 8):
        n = 2**p
        for k in range(n + 1):
            for s in range(n + 1):
                values = []
                for z in (0, 1):
                    for total in range(0, n, 2):
                        a = s - (k - n) * z - (k + 1) * (1 - z) - total
                        values.append(a * (a - 1) // 2)
                assert min(values) == int(s == k), (p, k, s)
                branch_evaluations += len(values)
    return recognized_profiles, branch_evaluations


def main():
    zero_sets, scaled, perturbed = check_recognition()
    profiles, evaluations = check_cardinality()
    print(f"PASS: reconstruction versus exhaustive generator oracle: {zero_sets} zero sets")
    print(f"PASS: rational scaling: {scaled}; incompatible positive-value perturbations: {perturbed}")
    print(f"PASS: exact-cardinality recognition p=2..8: {profiles} profiles")
    print(f"PASS: known upper construction p=2..7: {evaluations} branch evaluations")
    print("Finite checks supplement the analytical proofs; unrestricted optimality is not tested.")


if __name__ == "__main__":
    main()
