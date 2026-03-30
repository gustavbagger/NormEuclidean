# Computing least cubic non-residue cutoffs
Let q1 be the least cubic non-residue modulo a prime p. Then, this algorithm computes a cutoff N, non-decreasing with respect to p, such that:
If q1<= N, then the cyclic cubic field with conductor p is not Norm-Euclidean.
This computation is a crucial step in the computatation needed to completely classify all cyclic cubic Norm-Euclidean fields. See our preprint for more details [The determination of norm-Euclidean cyclic cubic fields](https://arxiv.org/abs/2507.05862).
## Quick Start
```bash
# install using Go toolchain
go install github.com/gustavbagger/NormEuclidean

# run
NormEuclidean <numerical precision (int)>
```