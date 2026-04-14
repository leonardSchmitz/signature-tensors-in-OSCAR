```@meta
CollapsedDocStrings=true
```


# SignatureTensors.jl

**SignatureTensors.jl** is a Julia package for computing and manipulating signature tensors of paths and membranes, built on top of the [OSCAR](https://www.oscar-system.org) computer algebra system.

> **Version:** 1.0.0.

---

## Mathematical Background

Path signatures, introduced by Terry Lyons in the framework of rough path theory
[lyons1998differential, lyons2002system, lyons2007differential](@cite),
provide a canonical and robust way to encode paths as algebraic objects via iterated integrals. This viewpoint has become fundamental in modern stochastic analysis and has been extensively developed in the literature, notably in the work of [friz2010multidimensional](@cite). Further foundational and comprehensive treatments can be found in the monographs of [lyons2007differential](@cite), and  [friz2020course](@cite). Their structure is deeply
connected to classical results on iterated integrals due to [chen1957integration, chen1977iterated](@cite). Beyond their analytical origins, signatures have found applications in a wide range
of areas, including mathematical finance [cuchiero2026](@cite), machine learning [chevyrev2026](@cite), and topological data analysis [chevyrev2018persistence](@cite).

## Overview

**SignatureTensors.jl** provides a general computational framework for working with truncated signatures of paths and membranes. It is designed to operate over both standard Julia numerical types (such as `Float64`, `Float32`, and other native floating-point formats) and exact algebraic structures in the [OSCAR](https://www.oscar-system.org) computer algebra system, including rings and polynomial rings. This dual design enables fast numerical experimentation as well as rigorous symbolic computation, leveraging OSCAR’s infrastructure for ring arithmetic, Gröbner bases, and Lie-theoretic methods.

The package supports:

- **Path signatures** — iterated-integrals signatures of piecewise linear, polynomial, spline, axis, and monomial paths over arbitrary OSCAR rings.
- **Membrane signatures** — two-parameter (id-)signatures of piecewise bilinear and polynomial membranes, extending the one-parameter theory [lotter2024signature](@cite).
- **Tensor learning** — recovery of paths from their signature tensors via polynomial systems and Gröbner bases, including the efficient algorithm from [schmitz2025efficientalgorithmtensorlearning](@cite).
- **Lie group barycenters** — computation of Fréchet means on the free nilpotent Lie group $G_{d,k}$, with several algorithm options including BCH-based and polynomial map approaches.
- **Algebraic operations** — group multiplication, inverse, logarithm, exponential, and graded projections on truncated tensor algebra elements.

---

## Installation

```julia
] activate .
] instantiate
```

Then load the package alongside OSCAR:

```julia
using Oscar
using SignatureTensors
```

---

## Quick Start

```julia
using Oscar, SignatureTensors

# Define a truncated tensor algebra: dimension d=2, truncation level k=3
d, k = 2, 3
T = TruncatedTensorAlgebra(QQ, d, k)

# Signature of the canonical axis path
C=sig(T, :axis)

# Signature of a polynomial path t ↦ (t + 2t², 3t + 4t²)
S=sig(T, :pwln, coef = QQ[1 2; 3 4])

# Path recovery from a signature tensor
recover(S,core=C)
```

---

## Documentation

```@contents
Pages = ["api.md"]
Depth = 2
```
---

## Authors

| Name | Affiliation | Email |
|------|-------------|-------|
| Gabriel Riffo | TU Berlin | [riffo@tu-berlin.de](mailto:riffo@tu-berlin.de)    |
| Leonard Schmitz | TU Berlin | [lschmitz@math.tu-berlin.de](mailto:lschmitz@math.tu-berlin.de) |

Funded by the Deutsche Forschungsgemeinschaft (DFG) — [CRC/TRR 388 *"Rough Analysis, Stochastic Dynamics and Related Fields"*](https://sites.google.com/view/trr388/) — Project A04 and B01, 516748464.

