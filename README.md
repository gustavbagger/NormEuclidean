# Norm-Euclidean cyclic cubic fields
This repo contains the code used for our unconditional classification of all norm-Euclidean cyclic cubic fields. The code is split in two portions: NonResidueCutoffs and NumericalVerification. See remark 11 and section 8 in our preprint [The determination of norm-Euclidean cyclic cubic fields](https://arxiv.org/abs/2507.05862) for more details.

## Numerical verification
This algorithm enumerates primes p in the range 2x10^14 - 2x10^20 with no small cubic non-residues modulo p. This computation verifies that there are no norm-Euclidean cyclic cubic fields with conductor p in the given range.

## Computing least cubic non-residue cutoffs
Let q1 be the least (prime) cubic non-residue modulo a prime p. Then, this algorithm computes a cutoff N, non-decreasing with respect to p, such that:
If q1<= N, then the cyclic cubic field with conductor p is not Norm-Euclidean.

### Quick Start
```bash
# install using Go toolchain
go install github.com/gustavbagger/NormEuclidean

# run
NormEuclidean <numerical precision (int)>
```