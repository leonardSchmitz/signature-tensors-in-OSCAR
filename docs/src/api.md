```@meta
CurrentModule = SignatureTensors
```


# Documentation



---

## Types


### `TruncatedTensorAlgebra`


The core algebraic structures of the package. `TruncatedTensorAlgebra` defines the ambient space
```math
T_{d,k} := \bigoplus_{\ell=0}^{k} (\mathbb{R}^d)^{\otimes \ell}
= \mathbb{R} \oplus \mathbb{R}^d \oplus \mathbb{R}^{d\times d} \oplus \dots \oplus \mathbb{R}^{d \times \dots \times d},
```

together with a coefficient ring and a sequence type (e.g., `:iis` for iterated-integrals signatures,
`:p2id` for two-parameter id-signatures). Elements of this algebra are represented by
`TruncatedTensorAlgebraElem`, which stores the graded components as a vector of arrays.

For a sufficiently smooth path ``X:[0,1]\to \mathbb{R}^d``, its **truncated iterated-integrals signature** is
```math
\sigma^{\le k}(X) = 1 \oplus \sigma^{(1)}(X) \oplus \dots \oplus \sigma^{(k)}(X) \in T_{d,k},
```

with entries
```math
\sigma^{(\ell)}(X)_{w_1,\dots,w_\ell} := \int_{0 \le t_1 \le \dots \le t_\ell \le 1}
\dot X_{w_1}(t_1)\cdots \dot X_{w_\ell}(t_\ell) \, \mathrm{d}t_1 \cdots \mathrm{d}t_\ell,
\quad 1\le w_j \le d.
```

By Chen–Chow's theorem, the image of ``\sigma^{\le k}`` lies in the free nilpotent Lie group
``\mathcal{G}_{d,k}`` [friz2010multidimensional; Theorem 7.28](@cite).



```@docs
TruncatedTensorAlgebra

```
The constructor of TruncatedTensorAlgebra{R} is flexible with respect to the choice of
coefficients. In particular, it supports standard Julia numeric types such as Float64
and other floating-point representations.
Moreover, it integrates with the algebraic structures available in Oscar.jl, allowing
the use of coefficient rings such as rational numbers or polynomial rings.

```@docs
truncation_level(::TruncatedTensorAlgebra)
base_dimension(::TruncatedTensorAlgebra)
base_algebra(::TruncatedTensorAlgebra)
sequence_type(::TruncatedTensorAlgebra)
```

### `TruncatedTensorAlgebraElem`

Accessor functions for the components of a `TruncatedTensorAlgebraElem`.
```@docs
TruncatedTensorAlgebraElem
Base.parent(::TruncatedTensorAlgebraElem)
tensor_sequence(::TruncatedTensorAlgebraElem)
Base.getindex(x::TruncatedTensorAlgebraElem, w...)
```

Standard algebraic operations on `TruncatedTensorAlgebraElem`. All operations respect the
truncated tensor algebra structure: multiplication is the shuffle/group product (truncated
at level $k$), and `exp`/`log` map between the Lie algebra and the Lie group $G_{d,k}$.

```@docs
Base.:+
Base.:-
Base.:*
Base.:^
Base.inv
Base.exp
Base.log
Base.vec
Base.:(==)
```

---

## Signature Constructors

The primary interface for computing signatures. `sig` dispatches on the `geom_type` symbol
to select the appropriate path or membrane class, and accepts keyword arguments for
coefficients, shape, composition, regularity, and algorithm choice. See the table below for
supported geometry types.

| `geom_type` | Description | Key arguments |
|-------------|-------------|---------------|
| `:point` | Constant path (unit in $G_{d,k}$) | — |
| `:axis` | Canonical axis path | — |
| `:mono` | Moment path | — |
| `:pwln` | Piecewise linear path | `coef`, `algorithm` (`:Chen` or `:congruence`) |
| `:poly` | Polynomial path | `coef`, `algorithm` (`:congruence`, `:ARS26` or `:LS26`) |
| `:pwmon` | Piecewise monomial path | `composition`, `regularity` |
| `:spline` | Piecewise polynomial (spline) path | `coef`, `composition`, `regularity` |
| `:segment` | Single linear segment | `coef` |
| `:axis` *(membrane)* | Canonical axis membrane | `shape` |
| `:mono` *(membrane)* | Monomial membrane | `shape` |
| `:pwbln` | Piecewise bilinear membrane | `coef`, `shape`, `algorithm` (`:congruence` or `:LS26`) |
| `:poly` *(membrane)* | Polynomial membrane | `coef`, `shape`|

> For membrane types, set `sequence_type=:p2id` when constructing [`TruncatedTensorAlgebra`](@ref).

```@docs
sig
```

---

## Tensor Learning / Path Recovery

Tools for the inverse problem of recovering a path from its signature tensor. `recover`
solves the polynomial system $S = A * C$ using Gröbner bases, where $C$ is a core tensor
(axis, polynomial, spline, or membrane). The optional `algorithm=:Sch25` selects the
efficient congruence-stabilizer method that scales as $O(d^4)$ in expectation.

```@docs
recover
```

---

## Barycenters

The expected signature of a family of paths is not necessarily itself a signature. Motivated by this, [clausel2024barycenterfreenilpotentlie,amendola2025learning](@cite) introduce the notion of a Lie group barycenter on signatures, following Buser and Karcher.

Given $S_1,\dots,S_N\in\mathcal{G}_{d,k}$, the equation

```math
0 = \sum_{i=1}^N \log\!\big(M^{-1}\cdot S_i\big)
```
defines a unique element $M\in\mathcal{G}_{d,k}$, called the *Lie group barycenter*.

By [amendola2025learning; Proposition 4.9](@cite), the solution map
```math
\mathsf{bary}:\mathcal{G}_{d,k}^N \to \mathcal{G}_{d,k}, \quad (S_1,\dots,S_N)\mapsto M
```
is a (non-commutative) polynomial surjection in the graded components. Our default implementation uses this polynomial map.


| `algorithm` | Description |
|-------------|-------------|
| *(default)* | Polynomial surjection map |
| `:geodesic` | Geodesic barycenter [amendola2025learning; Prop. 4.4](@cite) |
| `:AS25trunc2` | Truncation-level-2 formula  [amendola2025learning; Thm. 4.11](@cite) |
| `:CDMSSU24aBCH` | Asymmetrized BCH series [clausel2024barycenterfreenilpotentlie](@cite) |


```@docs
bary
```

---






