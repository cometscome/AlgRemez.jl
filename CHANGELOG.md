# Changelog

This project follows [Semantic Versioning](https://semver.org/).

## 0.1.0

- Add a native Julia implementation of the original AlgRemez calculation.
- Reproduce the legacy initial-point, exchange, Gaussian-elimination,
  Newton-deflation, and partial-fraction algorithms.
- Add an independently refined arbitrary-precision algorithm.
- Preserve the `calc_coefficients` and `AlgRemez_coeffs` RHMC interface.
- Add direct `AlgRemez_jll` coefficient and iteration regression tests.
- Add LatticeQCD-oriented action and molecular-dynamics tests.
- Add legacy command-line and `approx.dat`/`error.dat` output support.
