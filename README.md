# RLU

Rust LU Decomposition

## About

Crate `rlu` provides sparse LU factorization with partial pivoting as
described in "Sparse Partial Pivoting in Time Proportional to Arithmetic
Operations" by John R. Gilbert and Tim Peierls.

```
@article{Gilbert1988,
  doi = {10.1137/0909058},
  url = {https://doi.org/10.1137/0909058},
  year  = {1988},
  month = {sep},
  publisher = {Society for Industrial {\&} Applied Mathematics ({SIAM})},
  volume = {9},
  number = {5},
  pages = {862--874},
  author = {John R. Gilbert and Tim Peierls},
  title = {Sparse Partial Pivoting in Time Proportional to Arithmetic Operations},
  journal = {{SIAM} Journal on Scientific and Statistical Computing}
}
```

## License

This source code is distributed, with the permission of John Gilbert
and Tim Peierls, under the BSD 3-clause license ([LICENSE](LICENSE) or
https://opensource.org/licenses/BSD-3-Clause).

This source code was translated from the original `gp` FORTRAN code into 
[Rust](https://rustlang.org) by Richard W. Lincoln. The FORTRAN source was
distributed in Sivan Toledo's work on incomplete-factorization, from PARC 
in the early 1990s, and can be found in the `ILU` package on Netlib:

http://www.netlib.org/linalg/ilu.tgz
