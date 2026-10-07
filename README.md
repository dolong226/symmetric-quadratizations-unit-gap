# Verification programs for symmetric quadratizations of exact k-out-of-n functions

Companion code for the manuscript

> Long Duc Do, *Symmetric quadratizations of exact k-out-of-n functions: rigidity and the unit gap*.

The proofs in the manuscript do not depend on these programs. They are finite, independent
checks of the statements, and one bounded numerical search that the manuscript quotes as
evidence only.

| File | What it checks | Dependencies |
|------|----------------|--------------|
| `verify_theory.py` | The recognition procedure against an exhaustive oracle (13,369 zero sets); the classification at n = 2^p for p = 2..8; the known upper-bound construction for p = 2..7 | Python standard library |
| `verify_family.py` | The subset-sum construction with weights 3, 5, 7, 14, 28, ... for exact 1- and 2-out-of-n, m = 3..8, and maximality of the stated ranges for these weights | Python standard library |
| `search_normal_form.py` | Bounded mixed-integer search over the normal form of x-symmetric quadratizations; default run covers n = 14, 15 | `numpy`, `scipy >= 1.9` |

## Running

```
python verify_theory.py
python verify_family.py
pip install -r requirements.txt
python search_normal_form.py              # n = 14, 15, all k in the gap range
python search_normal_form.py 6,1,2 7,1,2  # arbitrary cases n,k,m
```

The first two programs use exact integer and rational arithmetic and finish in well under a
minute. The search can take several minutes per case.

## Scope and limits

- `verify_theory.py` and `verify_family.py` are exhaustive within the stated finite ranges.
  They do not replace the proofs for general sizes.
- `search_normal_form.py` uses floating-point mixed-integer programming with coefficients
  restricted to absolute value at most 3000. An `INFEASIBLE` answer is numerical evidence for
  that box only, not a proof of nonexistence.
- None of the programs tests quadratizations that are not symmetric in the original variables.

## License

MIT, see `LICENSE`.
