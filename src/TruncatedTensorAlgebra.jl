export TruncatedTensorAlgebra,
       TruncatedTensorAlgebraElem,
       truncation_level,
       base_dimension,
       sequence_type,
       base_algebra,
       tensor_sequence,
       zero,
       one,
       sig
       #matrix_tensorAlg_congruence_TA,

       #sigAxis_ClosedForm,  # sr
       #sigAxis_p2id_ClosedForm,  # sr
       #sigAxis_p2_ClosedForm,  # sr
       #sigAxis_Chen,  # sr
       #sigAxis_p2id_Chen,  # sr
       #sigAxis_p2_Chen,  # sr

       #Moment functions, 
       #moment_path_level,  # sr
       #moment_membrane,  # sr
       #moment_membrane_p2id,  # sr
       #moment_membrane_p2,  # sr

      # Tensor algebra operations
       #applyMatrixToTTA, 
       #sig2parPoly

"""
TruncatedTensorAlgebra{R}

Structure that stores a **truncated signature** of a rough path or a
two-parameter object (membrane).

Concrete signatures are stored as elements of this algebra via
[`TruncatedTensorAlgebraElem`](@ref).

# Fields
- `base_algebra::R`  
  Base algebra over which the signature is defined.

- `base_dimension::Int`  
  Dimension `d` of the underlying path or membrane.

- `truncation_level::Int`  
  Truncation level `k` of the signature.

- `sequence_type::Symbol`  
  Type of signature:
  - `:iis`   — iterated-integrals signature of a rough path,
  - `:p2id`  — two-parameter signature for membranes

# Mathematical meaning
This object represents the truncated signature space
 S^{k}(X) = (1, S^1(X), dots, S^k(X)), 
where the precise tensor structure depends on `sequence_type`.

See also [`TruncatedTensorAlgebraElem`](@ref).
"""




struct TruncatedTensorAlgebra{R}
    base_algebra::R
    base_dimension::Int
    truncation_level::Int
    sequence_type::Symbol
end

"""
    TruncatedTensorAlgebraElem{R,E}

Concrete **truncated signature** of a rough path or membrane.

This object stores the actual tensor coefficients of the signature up to
a fixed truncation level.

# Type Parameters
- `R`: base algebra of the parent signature object.
- `E`: type of the coefficients at each tensor level.

# Fields
- `parent::TruncatedTensorAlgebra{R}`  
  Signature metadata (dimension, truncation level, type).

- `elem::Vector{Array{E}}`  
  Sequence of tensors representing the signature:
    - `elem[:] level 1  (vector),
    - `elem[:,:]`: level 2,
    - `elem[:,:,:, ..., :]`: level k,

# Mathematical meaning
This object represents a truncated signature

S^{k}(X) = (1, S^1(X), ..., S^k(X)),

where each `S^i(X)` is a tensor of order `i`.
"""

struct TruncatedTensorAlgebraElem{R,E}
    parent::TruncatedTensorAlgebra{R}
    elem::Vector{Array{E}}
end

truncation_level(F::TruncatedTensorAlgebra) = F.truncation_level
base_dimension(F::TruncatedTensorAlgebra) = F.base_dimension
base_algebra(F::TruncatedTensorAlgebra) = F.base_algebra
sequence_type(F::TruncatedTensorAlgebra) = F.sequence_type

Base.parent(a::TruncatedTensorAlgebraElem) = a.parent
tensor_sequence(a::TruncatedTensorAlgebraElem) = a.elem


"""
    TruncatedTensorAlgebra(R, d::Int, k::Int; sequence_type = :iis)

Create a container that stores **truncated signatures** of rough paths or
two-parameter objects (membranes).

# Arguments
- `R`  
  Base algebra of the signature coefficients.

- `d::Int`  
  Dimension of the underlying path or membrane.
  Must satisfy `d ≥ 0`.

- `k::Int`  
  Truncation level of the signature.
  Must satisfy `k ≥ 0`.

# Keyword Arguments
- `sequence_type::Symbol = :iis`  
  Type of signature to be stored:
  - `:iis`   — iterated-integrals signature of a rough path,
  - `:p2id`  — signature of a membrane,

# Returns
- `TruncatedTensorAlgebra{typeof(R)}`

# Errors
- Throws an error if `d < 0`.
- Throws an error if `k < 0`.
- Throws an error if `sequence_type` is not one of
  `(:iis, :p2id, :p2)`.

# Example

A = TruncatedTensorAlgebra(QQ, 2, 3; sequence_type = :iis)

This creates a container for truncated signatures of a 2-dimensional
rough path, truncated at level 3.
"""
function TruncatedTensorAlgebra(R, d::Int, k::Int; sequence_type::Symbol=:iis)

    d >= 0 || error("ambient dimension must be >= 0")
    k >= 0 || error("truncation level must be >= 0")

    if !(sequence_type in (:iis, :p2id, :p2))
        error("sequence_type must be one of :iis, :p2id, :p2")
    end

    A = TruncatedTensorAlgebra{typeof(R)}(R, d, k, sequence_type)

    return A
end


function _C0_TA(_k::Int, _order::Int, _alg, is_one::Bool)
    if _k == 0
        return fill(is_one ? one(_alg) : zero(_alg), ())
    else
        return zeros(_alg, ntuple(_ -> _order, _k)...)
    end
end

function tensor_sequence_constructor(A::TruncatedTensorAlgebra; is_one=true)
    R   = A.base_algebra
    d   = A.amb_dim
    k   = A.trunc_level
    seq = A.sequence_type

    elems = Vector{Any}(undef, k+1)

    # Level 0 (common for all Algebra)
    elems[1] = fill(is_one ? one(R) : zero(R), ())

    # IIS
    if seq == :iis
        for n in 1:k
            elems[n+1] = zeros(R, ntuple(_ -> d, n)...)
        end
        return elems
    end

    # P2ID
    if seq == :p2id
        for n in 1:k
            dims = ntuple(_ -> d, n)
            elems[n+1] = Array{Any}(undef, dims...)
        end
        return elems
    end

    # P2
    if seq == :p2
        for n in 1:k
            dims = (ntuple(_ -> d, n)..., factorial(n))
            elems[n+1] = Array{Any}(undef, dims...)
        end
        return elems
    end

    error("Unknown sequence_type = $seq")
end



function _C0_seq_TA(_trunc_level::Int, _order::Int, _alg, is_one::Bool)
    return [_C0_TA(_k, _order, _alg, is_one) for _k in 0:_trunc_level]
end


#### Base Constructor

"""
    Base.zero(T::TruncatedTensorAlgebra)

Return the **zero truncated signature** associated with the truncated tensor
algebra `T`.

This function constructs a `TruncatedTensorAlgebraElem` whose tensor
coefficients are identically zero at all levels. It represents the additive
identity in the space of truncated signatures.

# Arguments
- `T::TruncatedTensorAlgebra`  

# Returns
- `TruncatedTensorAlgebraElem`  
  An element of `T` whose tensor coefficients are zero at every level.

# Example
T = TruncatedTensorAlgebra(QQ, 3, 2; sequence_type = :iis)
z = zero(T)
"""


function Base.zero(T::TruncatedTensorAlgebra{R}) where R

    k = truncation_level(T)
    d = base_dimension(T)
    RA = base_algebra(T)   
    seq = T.sequence_type
    E = typeof(one(RA))

    elem = Vector{Array{E}}(undef, k + 1)

    # Level 0 (always)
    elem[1] = fill(zero(RA), ())   # 0-dimensional scalar

    # Levels >= 1
    if seq == :iis
        for n in 1:k
            elem[n + 1] = zeros(RA, ntuple(_ -> d, n)...)
        end

    elseif seq == :p2id
        for n in 1:k
            if n == 1
                dims = (d,)
            else
                dims = ntuple(_ -> d, n)
            end
            elem[n + 1] = zeros(RA, dims...)
        end

    elseif seq == :p2
        for n in 1:k
            if n == 1
                dims = (d,)
            else
                dims = (ntuple(_ -> d, n)..., factorial(n))
            end
            elem[n + 1] = zeros(RA, dims...)
        end

    else
        throw(ArgumentError("zero(T) not implemented for sequence_type = $seq"))
    end

    return TruncatedTensorAlgebraElem{R,E}(T, elem)
end

"""
    one(T::TruncatedTensorAlgebra)

Return the **identity truncated signature** associated with the truncated tensor
algebra `T`.


# Arguments
- `T::TruncatedTensorAlgebra`  

# Returns
- `TruncatedTensorAlgebraElem`  
  An element of `T` whose tensor coefficients are the identity signature.

# Example
T = TruncatedTensorAlgebra(QQ, 3, 2; sequence_type = :iis)
z = one(T)
"""


function Base.one(T::TruncatedTensorAlgebra{R}) where R
    if T.sequence_type == :iis

        k = truncation_level(T)
        d = base_dimension(T)
        alg = T.base_algebra
        E = typeof(one(alg))

        return TruncatedTensorAlgebraElem{R,E}(T, _C0_seq_TA(k, d, alg, true))
    else
        throw(ArgumentError("one(T) is only defined for sequence_type = :iis"))
    end
end


"""
    sig(
        T::TruncatedTensorAlgebra,
        path_type::Symbol;
        coef = [],
        composition::Vector{Int} = Int[],
        regularity::Int = 0,
        algorithm::Symbol = :default
    )

Compute the **truncated signature** associated with a specified path type
inside the truncated tensor algebra `T`.

# Arguments
- `T::TruncatedTensorAlgebra`  
  Truncated tensor algebra specifying the dimension, truncation level,
  base algebra, and signature type.

- `path_type::Symbol`  
  Type of path whose signature is to be computed. Supported values include:
  
  # Path types and algorithms

 path_type = :point
   Creates the signature of a constant path.
   This corresponds to the unit signature, where the level-0 component is 1
   and all higher levels are zero.

 path_type = :axis
   Creates the signature of an axis path, i.e. a path that moves only along
   coordinate axes.

   Additional arguments:
     algorithm = :AFS19
       Computes the signature explicitly using closed-form formulas.

     algorithm = :Chen
       Computes the signature using Chen’s identity.

  path_type = :mono
   Computes the signature of a monomial path.

 path_type = :pwln
 Computes the signature of a piecewise linear path.

   Additional arguments:
     coef
       Matrix of coefficients describing the piecewise linear path.

     algorithm = :Chen
       Computes the signature using Chen’s identity.

     algorithm = :congruence
       Computes the signature using matrix tensor congruence.

  path_type = :pwmon
"""
sig
function sig(T::TruncatedTensorAlgebra{R},
             path_type::Symbol; 
             coef=[], 
             composition::Vector{Int}=Int[],
             regularity::Int=0,
             algorithm::Symbol=:default) where R
    if path_type==:point && coef==[] && algorithm == :default
        return one(T)
    elseif path_type==:axis && coef==[] && (algorithm == :default || algorithm ==:AFS19)
        return sigAxis_TA_ClosedForm(T) 
    elseif path_type==:axis && coef==[] && algorithm == :Chen 
        return sig_axis_TA(T) 
    elseif path_type==:mono && coef==[] && algorithm == :default
        return sig_mono_TA(T) 
    elseif path_type==:pwln && algorithm == :congruence
        return sig_pwln_TA_Congruence(T,Array(coef))
    elseif path_type==:pwln && (algorithm == :Chen || algorithm == :default)
        return sig_pwln_TA_chen(T,Array(coef))
    elseif (path_type == :pwmon && algorithm == :Chen)
        return sig_pw_mono_chen(T,composition,regularity)
    elseif path_type == :pwmon && (algorithm == :ALS26 || algorithm == :default)
        return sig_pw_mono_ALS26(T,composition,regularity)
    else 
        throw(ArgumentError("sig not supported for given arguments")) 
    end 
end


function Base.show(io::IO, x::TruncatedTensorAlgebraElem)
    for (i, t) in enumerate(x.elem)
        show(io, MIME("text/plain"), t)
        if i < length(x.elem)
            println(io, "\noplus")
        end
    end
end

function Oscar.evaluate(
    p::FreeTruncSigAlgMultivElem{T},
    bs::Vector{TruncatedTensorAlgebraElem{R,E}}
) where {T,R,E}

    A = parent(bs[1])                 # algebra base
    k = truncation_level(A)           # truncation level

    # evaluate in every level 0..k
    res_seq = [evaluate_to_tensor_TA(p, bs, i) for i in 0:k]

    return TruncatedTensorAlgebraElem(A, res_seq)
end




# helpers for sig_mono 

function _CmonoTA(_k::Int, _order::Int, R0)
    if _k == 0
        return fill(one(R0), ())
    end

    E = typeof(one(R0))

    res = Array{E}(undef, ntuple(_ -> _order, _k)...)

    for idx in CartesianIndices(res)
        word = Tuple(idx)
        coeff = prod(word)
        expo  = prod(cumsum(word))
        res[idx] = QQ(coeff, expo) * one(R0)
    end

    return res
end


function _Cmono_seqTA(_trunc_level::Int, _order::Int, R0)
    return [_CmonoTA(_k, _order, R0) for _k in 0:_trunc_level]
end

function _one_hot_TA(_i::Int, _n::Int, _R)
    res = fill(zero(_R), _n)
    res[_i] = one(_R)
    return res
end


function sig_mono_TA(T::TruncatedTensorAlgebra{R}) where R
    if T.sequence_type != :iis
        throw(ArgumentError("sig_mono only defined for sequence_type = :iis"))
    end

    k = truncation_level(T)
    d = base_dimension(T)
    R0 = base_algebra(T)
    # seq must be an Vector{Array{E}} for some E
    seq = _Cmono_seqTA(k, d, R0)

    # determinate E using seq
    E = eltype(seq[1])   # type of the elements of the intern array
    return TruncatedTensorAlgebraElem{R, E}(T, seq)
end

#helpters for Pw polynomial paths

function embedding_matrix(m::Vector{Int},i::Int)
  M = sum(m)
  A = zero_matrix(QQ,M,m[i])
  A[sum(m[1:i-1])+1:sum(m[1:i]),:] = identity_matrix(QQ,m[i])
  return A
end

function nextBlock(A, r::Int, mi)
    nrows, ncols = size(A)
    vecs = zeros(QQ, ncols, r)
    for i in 1:r
        for j in 1:ncols
            vecs[j, i] = binomial(j, i)
        end
    end
    return matrix(QQ,Matrix(A) * Matrix(vecs))
end

function coreSplineTrafo(m::Vector{Int}, r::Int) 
    total_dim = sum(m)
    if r == 0
        return identity_matrix(QQ,total_dim)
    end
    zz = 1
    crb = identity_matrix(QQ,m[1])
    B = crb
    while zz < length(m)
        nr = size(B, 2)
        start_col = nr - m[zz] + 1
        end_col   = nr
        subB = B[:, start_col:end_col]
        nb = nextBlock(subB, r, zz)
        idpart = identity_matrix(QQ,m[zz+1] - r)
        B = block_diagonal_matrix([hcat(B, nb), idpart])
        zz += 1
    end
    return B
end

function sig_pw_mono_chen(TTS::TruncatedTensorAlgebra,m,r) 
  k = truncation_level(TTS)
  R = base_algebra(TTS)
  d = base_dimension(TTS)
  if d != sum(m) - r*(length(m) -1)
    error("m must be a composition of the ambient dimension")
  end
  sigs = [Array(embedding_matrix(m,i))*sig(TruncatedTensorAlgebra(R,m[i],k),:mono) for i in (1:length(m))]
  res = prod(sigs)
  if r == 0
    return res
  else
    return Array(coreSplineTrafo(m,r))*res
  end
end

function adapted_word(w::Vector{Int}, m::Vector{Int})
    ell = length(m)
    offsets = cumsum([0; m[1:end-1]])  
    blocks = [Int[] for _ in 1:ell]
    current = 1
    for x in w
        while current <= ell && x > offsets[current] + m[current]
            current += 1
        end
        if current > ell || x ≤ offsets[current]
            return false,[]
        end
        push!(blocks[current], x)
    end
    # subtract offsets to get words over [mi]
    return true,[blocks[i] .- offsets[i] for i in 1:ell]
end

function _Cpwpoly(_k::Int, m::Vector{Int}, _R)
   M = sum(m)
   res = zeros(_R, M*ones(Int,_k)...)
   for idx in CartesianIndices(res)
       if _k ==0 
         res[idx] = one(_R)
       else
         b,vs = adapted_word(collect(Tuple(idx)),m)
         if b
           res[idx] = prod([QQ(prod(v),prod(cumsum(v))) for v in vs])*one(_R)
         end
       end 
   end
   return res
end

function _Cpwpoly_seq(_trunc_level::Int, m::Vector{Int}, _R)
   Cmono = [_Cpwpoly(i,m, _R) for i in (0:_trunc_level)]
   return Cmono
end

function sig_pw_mono_ALS26(TTS::TruncatedTensorAlgebra,m::Vector{Int},r)
  d = base_dimension(TTS)
  if d != sum(m) - r*(length(m) -1)
    error("m must be a composition of the ambient dimension")
  end
  k = truncation_level(TTS) 
  R = base_algebra(TTS)
  TM = TruncatedTensorAlgebra(R,sum(m),k)
  res = TruncatedTensorAlgebraElem(TM,_Cpwpoly_seq(k,m,R))
  if r == 0
    return res
  else
    return Array(coreSplineTrafo(m,r))*res
  end
end


# Base functions

function Base.getindex(x::TruncatedTensorAlgebraElem, w...)
    k=length(w)
    b=parent(x)
    typeT=b.sequence_type
    if k < 0
        throw(BoundsError(x, k))
    end
    
    if k > length(x.elem)
        throw(BoundsError(x, k))
    end
    if typeT==:iis 
        temp= x.elem[k+1 ]
    elseif typeT==:p2
        if k==1
            temp=x.elem[k+1]
        else
            temp=x.elem[k]
        end
    elseif typeT==:p2id
        if k==1
            temp=x.elem[k+1]
        else
            temp=x.elem[k]
        end
    else 
        error("seq must be one of :iis, :p2id, :p2")
    end
    
    return temp[w...]
end


function Base.vec(x::TruncatedTensorAlgebraElem)
    return vcat([vec(xi) for xi in x.elem]...)
end


function Base.:(==)(a::TruncatedTensorAlgebraElem, b::TruncatedTensorAlgebraElem)
    if parent(a).sequence_type == :iis && parent(b).sequence_type == :iis
        return parent(a) == parent(b) && tensor_sequence(a) == tensor_sequence(b)
    else
        throw(ArgumentError("== only defined for sequence_type == :iis"))
    end
end



# =========================
# Addition and Subtraction
# =========================

"""
    Base.:+(a::TruncatedTensorAlgebraElem, b::TruncatedTensorAlgebraElem)

Element-wise addition of two truncated signatures.

- Only defined for `sequence_type == :iis`.
- Returns a new `TruncatedTensorAlgebraElem` whose tensors are the sum of
  the corresponding levels of `a` and `b`.
"""
function Base.:+(a::TruncatedTensorAlgebraElem, b::TruncatedTensorAlgebraElem)
    if parent(a).sequence_type == :iis && parent(b).sequence_type == :iis
        A = parent(a)
        res_seq = tensor_sequence(a) + tensor_sequence(b)
        return TruncatedTensorAlgebraElem(A, res_seq)
    else
        throw(ArgumentError("+ only defined for sequence_type == :iis"))
    end
end


"""
    Base.:-(a::TruncatedTensorAlgebraElem, b::TruncatedTensorAlgebraElem)

Element-wise subtraction of two truncated signatures.

- Only defined for `sequence_type == :iis`.
- Returns a new `TruncatedTensorAlgebraElem` whose tensors are the difference
  of the corresponding levels of `a` and `b`.
- Analogous to `+` but performs subtraction.
"""
function Base.:-(a::TruncatedTensorAlgebraElem, b::TruncatedTensorAlgebraElem)
    if parent(a).sequence_type == :iis && parent(b).sequence_type == :iis
        A = parent(a)
        res_seq = tensor_sequence(a) - tensor_sequence(b)
        return TruncatedTensorAlgebraElem(A, res_seq)
    else
        throw(ArgumentError("- only defined for sequence_type == :iis"))
    end
end




# =========================
# Multiplication (Chen / Congruence / Scalar)
# =========================

"""
    Base.:*(a::TruncatedTensorAlgebraElem, b::TruncatedTensorAlgebraElem)

Chen product (truncated tensor algebra multiplication) of two signatures.

- Only defined for `sequence_type == :iis`.
- Uses free truncated algebra and evaluates the product using Chen’s identity.
- Returns a new `TruncatedTensorAlgebraElem` representing the concatenation
  of `a` and `b`.
"""

function Base.:*(a::TruncatedTensorAlgebraElem{R,E}, 
                 b::TruncatedTensorAlgebraElem{R,E}) where {R,E}
    A = parent(a)
    k = truncation_level(A)

    F, s = free_trunc_sig_alg_multiv(k, 2)
    chen = prod([free_sig_from_sample(i, F) for i in 1:2])

    return evaluate(chen, [a, b])
end

"""
    Base.:*(matrix::AbstractMatrix, x::TruncatedTensorAlgebraElem)

Matrix-tensor congruence multiplication.

- Only defined for `sequence_type == :iis`.
- Left-multiplies a truncated signature by a matrix using a tensor congruence
  operation.
"""

function Base.:*(matrix::AbstractMatrix, x::TruncatedTensorAlgebraElem)
    T = parent(x)

    if T.sequence_type == :iis
            return matrix_tensorAlg_congruence_TA(matrix, x)
    else
                throw(ArgumentError("matrix * TruncatedTensorAlgebraElem ony defined for sequence_type = :iis"))
    end

end


"""
    Base.:*(a::Number, b::TruncatedTensorAlgebraElem)

Scalar multiplication of a truncated signature.

- Only defined for `sequence_type == :iis`.
- Multiplies each tensor level by the scalar `a`.
- Supports numbers (`Number`) or algebra elements (`FieldElem` / `R`).
- Returns a new `TruncatedTensorAlgebraElem`.
"""


function Base.:*(a::FieldElem, b::TruncatedTensorAlgebraElem)
    if parent(b).sequence_type == :iis
        A = parent(b)
        k = truncation_level(A)
        res_seq = tensor_sequence(b)
        res_seq[1][] = a * res_seq[1][]        # Level 0
        for j in 2:k+1
            res_seq[j] = a .* res_seq[j]
        end
        return TruncatedTensorAlgebraElem(A, res_seq)
    else
        throw(ArgumentError("scalar * only defined for sequence_type == :iis"))
    end
end

#TODO is this redundant if R == QQ? 
function Base.:*(a::R, b::TruncatedTensorAlgebraElem{R,E}) where {R,E}
    if parent(b).sequence_type == :iis
       A = parent(b)
       k = truncation_level(A)
       res_seq = tensor_sequence(b)
       res_seq[1][] = a * res_seq[1][]
       for j in 2:k+1
           res_seq[j] = a .* res_seq[j]
       end
       return TruncatedTensorAlgebraElem(A, res_seq)
    else
        throw(ArgumentError("scalar * only defined for sequence_type == :iis"))
    end
end

function Base.:*(a::TruncatedTensorAlgebraElem, b::TruncatedTensorAlgebraElem)
    if parent(a).sequence_type == :iis && parent(b).sequence_type == :iis
        A = parent(a)
        k = truncation_level(A)
        F, s = free_trunc_sig_alg_multiv(k, 2)
        chen = prod([free_sig_from_sample(i, F) for i in 1:2])
        return evaluate(chen, [a, b])
    else
        throw(ArgumentError("* only defined for sequence_type == :iis"))
    end
end


# =========================
# Inverse and Powers
# =========================

"""
    Base.:inv(a::TruncatedTensorAlgebraElem)

Multiplicative inverse (Chen inverse) of a truncated signature.

- Only defined for `sequence_type == :iis`.
- Returns `a^{-1}` in the truncated tensor algebra using free algebra inversion.
"""

function Base.:inv(a::TruncatedTensorAlgebraElem)
    if parent(a).sequence_type == :iis
        A = parent(a)
        k = truncation_level(A)
        F, s = free_trunc_sig_alg_multiv(k, 1)
        p = inv(free_sig_from_sample(1, F))
        return evaluate(p, [a])
    else
        throw(ArgumentError("inv only defined for sequence_type == :iis"))
    end
end


"""
    Base.:^(a::TruncatedTensorAlgebraElem, n::Int)

Exponentiation of a truncated signature.

- Only defined for `sequence_type == :iis`.
- For `n >= 0`, returns `a^n` (n-fold Chen product).
- For `n < 0`, uses the inverse: `a^n = (inv(a))^(-n)`.
"""

function Base.:^(a::TruncatedTensorAlgebraElem, n::Int)
    if parent(a).sequence_type == :iis
        A = parent(a)
        res = one(A)
        temp = n >= 0 ? a : inv(a)

        for i in 1:abs(n)
            res = res * temp
        end
        return res
    else
        throw(ArgumentError("^ only defined for sequence_type == :iis"))
    end
end

# =========================
# Exponential and Logarithm
# =========================

"""
    Base.:exp(a::TruncatedTensorAlgebraElem)

Exponential of a truncated signature.

- Only defined for `sequence_type == :iis`.
- Computes the exponential in the free truncated tensor algebra:
  
    exp(a) = sum_{n=0}^{infty} frac{a^n}{n!}
  
  truncated to the algebra’s level.
"""


function Base.:exp(a::TruncatedTensorAlgebraElem)
    if parent(a).sequence_type == :iis
        A = parent(a)
        k = truncation_level(A)
        F, s = free_trunc_sig_alg_multiv(k, 1)
        p = exp(free_sig_from_sample(1, F) - one(F))
        return evaluate(p, [a])
    else
        throw(ArgumentError("exp only defined for sequence_type == :iis"))
    end
end



"""
    Base.:log(a::TruncatedTensorAlgebraElem)

Logarithm of a truncated signature.

- Only defined for `sequence_type == :iis`.
- Computes the logarithm in the free truncated tensor algebra:
  
    log(a) = sum_{n=1}^{infty} (-1)^{n+1} \frac{(a - 1)^n}{n}
  
  truncated to the algebra’s level.
"""


function Base.:log(a::TruncatedTensorAlgebraElem)
    if parent(a).sequence_type == :iis
        A = parent(a)
        k = truncation_level(A)
        F, s = free_trunc_sig_alg_multiv(k, 1)
        p = log(free_sig_from_sample(1, F))
        return evaluate(p, [a])
    else
        throw(ArgumentError("log only defined for sequence_type == :iis"))
    end
end







function ideal_of_lyndon_entries(x::TruncatedTensorAlgebraElem)
  A = parent(x)
  R = base_algebra(A)
  d = ambient_dimension(A)
  k = truncation_level(A)
  lynd = generate_lyndon_words(k, Vector((1:d)))
  #res = ideal(R,zero(R))
  elem_tts = tensor_sequence(x)
  return ideal(R,[elem_tts[length(w)+1][w...] for w in lynd])
end


# common_ring: verify the rule
function common_ring(R_tensor, R_matrix)
    if R_matrix isa MPolyRing
        return R_matrix
    else
        return R_tensor
    end
end


# Help function for permutations
function permutations_ntuple(v::NTuple{N,Int}) where N
    if N == 1
        return [v]
    end
    out = NTuple{N,Int}[]
    for i in 1:N
        x = v[i]
        rest = ntuple(k -> v[k < i ? k : k+1], N-1)
        for p in permutations_ntuple(rest)
            push!(out, (x, p...))
        end
    end
    return out
end

permutations_1_to_j(j::Int) =
    permutations_ntuple(ntuple(identity, j))


function sig_segment_TA(T::TruncatedTensorAlgebra{R}, v::Vector{E}) where {R,E}
    # --- validation---
    k = truncation_level(T)
    d = base_dimension(T)

    @assert length(v) == d "dimensions do not match"

    # --- veryfy the types---
    @assert typeof(v[1]) <: E "element types do not match"
    
    # --- create a new algebra ---
    T1 = TruncatedTensorAlgebra(base_algebra(T), 1, k, :iis)

    # --- signature mono in dimension 1 ---
    C = sig_mono_TA(T1)

    # --- Convert the vector v in a column matrix d×1 ---
    v_mat = reshape(v, d, 1)

    # --- Apply the matricial congruence (column matrix* signature 1D) ---
    return matrix_tensorAlg_congruence_TA(v_mat, C)
end


function sig_segment_standard_direction_TA(T::TruncatedTensorAlgebra{R}, _i::Int) where R
    d = base_dimension(T)
    A = base_algebra(T)
    # vector one-hot in the direction _i
    v = _one_hot_TA(_i, d, A)
    return sig_segment_TA(T, v)
end




function leading_coefficient_and_zero_TA(p)
    if is_zero(p)
        return leading_coefficient(one(parent(p))) - leading_coefficient(one(parent(p))) # TODO: its zero 
    else
        return leading_coefficient(p)
    end
end


#TODO Gabriel , incorperate tensor concatenation 
function concatenate_tensors_TA(t1::TruncatedTensorAlgebraElem, t2::TruncatedTensorAlgebraElem)
    # Check same parent
    @assert parent(t1) === parent(t2) "Parents must match"

    # Concatenate each tensor level-wise
    elem1 = t1.elem
    elem2 = t2.elem

    @assert length(elem1) == length(elem2) "Truncation levels must match"

    new_elem = Vector{Any}(undef, length(elem1))

    for k in 1:length(elem1)
        tensor1 = elem1[k]
        tensor2 = elem2[k]

        reshaped_tensor1 = reshape(tensor1, size(tensor1)..., ones(Int, ndims(tensor2))...)
        reshaped_tensor2 = reshape(tensor2, ones(Int, ndims(tensor1))..., size(tensor2)...)

        new_elem[k] = reshaped_tensor1 .* reshaped_tensor2
    end

    return TruncatedTensorAlgebraElem(parent(t1), new_elem)
end


function evaluate_to_tensor_TA(
    p::FreeTruncSigAlgMultivElem{T},
    _bs::Vector{TruncatedTensorAlgebraElem{R,E}},
    index::Int
) where {T,R,E}

    bs = deepcopy(_bs)

    # Level 0
    if index == 0
        res = tensor_sequence(bs[1])[1]   # Level 0
        res[] = leading_coefficient_and_zero_TA(
                  FreeTruncSigAlgMultivElem_to_polynomial(graded_component(p, 0))
               ) * one(res[])
        return res
    end

    # data of the space
    num_vars = length(bs)
    A = parent(bs[1])
    k = truncation_level(A)

    # Component grade index of p
    fi = FreeTruncSigAlgMultivElem_to_polynomial(graded_component(p, index))

    # Sequence of tensor of livels 1..k
    seq_tensor_seq = [tensor_sequence(b)[2:k+1] for b in bs]

    # Tensor result (Right size)
    res = zero(seq_tensor_seq[1][index])

    # Evaluate Term to Term
    for t in 1:length(fi)
        fit = collect(terms(fi))[t]

        first_exp  = fit.exps[1][1]
        sample_var = div(first_exp-1, k) + 1
        tensor_var = mod(first_exp-1, k) + 1

        temp = leading_coefficient_and_zero_TA(fit) .* seq_tensor_seq[sample_var][tensor_var]

        for fak in 2:length(fit.exps[1])
            temp_exp   = fit.exps[1][fak]
            sample_var = div(temp_exp-1, k) + 1
            tensor_var = mod(temp_exp-1, k) + 1
            temp = concatenate_tensors_TA(temp, seq_tensor_seq[sample_var][tensor_var])
        end

        res += temp
    end

    return res
end



function evaluate_TA(
    p::FreeTruncSigAlgMultivElem{T},
    bs::Vector{TruncatedTensorAlgebraElem{R,E}}
) where {T,R,E}

    A = parent(bs[1])                     # algebra base
    k = truncation_level(A)               # truncation level
    
    # Evaluation for every level 0..k
    res_seq = [evaluate_to_tensor_TA(p, bs, i) for i in 0:k]

    return TruncatedTensorAlgebraElem(parent(bs[1]), res_seq)
end

function combinations_with_replacement(iter, k)
    arr = collect(iter)
    n = length(arr)

    if k == 0
        return [[]]
    end

    res = Vector{Vector{Int}}()
    comb = zeros(Int, k)

    function backtrack(pos, start)
        if pos > k
            push!(res, copy(comb))
            return
        end
        for i in start:n
            comb[pos] = arr[i]
            backtrack(pos + 1, i)
        end
    end

    backtrack(1, 1)
    return res
end



function sigAxis_TA_ClosedForm(T::TruncatedTensorAlgebra{R}) where R
    if T.sequence_type != :iis
        error("sigAxis_TA only defined for sequence_type = :iis")
    end

    d = base_dimension(T)
    k = truncation_level(T)
    A = base_algebra(T)

    E = typeof(one(A))                 # tipe of the element
    elem_out = Vector{Array{E}}(undef, k+1)

    # Level 0 (tensor 0-dimensional)
    elem_out[1] = fill(one(A), ())    # array 0-dimensional

    # Level 1
    elem_out[2] = fill(one(A), d)

    # Levels >= 2
    for j in 2:k
        dims = ntuple(_ -> d, j)
        tensor_j = zeros(A, dims...)

        for idx in combinations_with_replacement(1:d, j)
            counts = Dict{Int,Int}()

            for i in idx
                counts[i] = get(counts, i, 0) + 1
            end

            denom = prod(factorial(c) for c in values(counts))
            value = 1 // denom

            tensor_j[idx...] = try
                A(value)
            catch
                convert(E, value)
            end
        end

        elem_out[j+1] = tensor_j
    end

    return TruncatedTensorAlgebraElem(T, elem_out)
end


#TODO: iterate chen not simpultaneous
function sig_axis_TA(T::TruncatedTensorAlgebra{R}) where R
    k = truncation_level(T)
    d = base_dimension(T) 
    # 1) Create the d standard segments  ( base direction)
    sample = [sig_segment_standard_direction_TA(T, i) for i in 1:d]
    # 2) create the free truncated multivariate algebra
    F, s = free_trunc_sig_alg_multiv(k, d)
    # 3) chen product: Multiplication the generators 1..d
    chen = prod([free_sig_from_sample(i, F) for i in 1:d])
    # 4) evaluation
    return evaluate_TA(chen, sample)
end


function sig_poly_TA(T::TruncatedTensorAlgebra{R}, coeffs::AbstractMatrix) where R
    # 1) Obtain moment path (sig_mono_TA)
    mono_path = sig_mono_TA(T)

    # 2) Apply coeficients with matricial congruence
    poly_path = matrix_tensorAlg_congruence_TA(coeffs, mono_path)

    return poly_path
end


function sig_pwln_TA_Congruence(T::TruncatedTensorAlgebra{R}, coeffs::AbstractMatrix{E}) where {R,E}
    d = base_dimension(T)
    @assert size(coeffs,1) == d "Dimensions mismatch"
    k=truncation_level(T)
    A = base_algebra(T)
    m=size(coeffs, 2)
    Tm = TruncatedTensorAlgebra(A,m,k)
    return matrix_tensorAlg_congruence_TA(coeffs, sig(Tm,:axis))
end

#TODO: iterate this, not all secments simpultaneous
function sig_pwln_TA_chen(T::TruncatedTensorAlgebra{R}, P::AbstractMatrix{E}) where {R,E}
    d = base_dimension(T)
    @assert size(P,1) == d "Dimensions mismatch"
    m=size(P,2)
    #k1=truncation_level(T)
    #R1 = base_algebra(T)
    #d1=size(P,2)
    #T1 = TruncatedTensorAlgebra(R1,d1,k1)
    #seg_vecs = [P[i+1, :] .- P[i, :] for i in 1:size(P,1)-1]
    #seg_sigs = [sig_segment_TA(T1, seg_vecs[i]) for i in 1:length(seg_vecs)]
    #TODO: the differences are more an interpoation and should be another option in sig 
    seg_sigs = [sig_segment_TA(T, P[:,i]) for i in 1:m]
    return prod(seg_sigs)
end

function matrix_tensorAlg_congruence_TA(
    v::AbstractVector,
    b::TruncatedTensorAlgebraElem{R,E}
) where {R,E}
    M = reshape(v, 1, length(v))
    return matrix_tensorAlg_congruence_TA(M, b)
end


function matrix_tensorAlg_congruence_TA(
    matrix::AbstractMatrix,
    b::TruncatedTensorAlgebraElem{R,E}
) where {R,E}
    T = parent(b)
    k = truncation_level(T)
    m = base_dimension(T)
    @assert m == ncols(matrix)
    d = nrows(matrix)
    R_tensor = base_algebra(T)
    R_matrix = parent(matrix[1,1])
    Rnew = common_ring(R_tensor, R_matrix)
    Tnew = TruncatedTensorAlgebra(Rnew, d, k, :iis)
    tensorSeq = tensor_sequence(b)
    tensorSeq = [map(x -> Rnew(x), t) for t in tensorSeq]
    resSeq = [
        matrix_tensor_congruence_TA(map(x -> Rnew(x), matrix), t)
        for t in tensorSeq
    ]
    return TruncatedTensorAlgebraElem(Tnew, resSeq)
end


# -------------------------------
# Principal function (dispatch)
# -------------------------------
function sigAxis_ClosedForm(T::TruncatedTensorAlgebra{R}, m::Int, n::Int) where R
    if m * n != base_dimension(T)
        error("m * n != d")
    end

    if T.sequence_type == :p2id
        return sigAxis_p2id_ClosedForm(T, m, n)
    elseif T.sequence_type == :p2
        return sigAxis_p2_ClosedForm(T, m, n)
    else
        error("sequence_type must be :p2id or :p2")
    end
end


# -------------------------------
# seq_type == :p2id
# -------------------------------
function sigAxis_p2id_ClosedForm(T::TruncatedTensorAlgebra{R}, m::Int, n::Int) where R
    k = truncation_level(T)
    
    if m * n != d
        error("m * n != d")
    end

    Tm = TruncatedTensorAlgebra(base_algebra(T), m, k, :iis)
    Tn = TruncatedTensorAlgebra(base_algebra(T), n, k, :iis)

    sig_m = sigAxis_TA_ClosedForm(Tm)
    sig_n = sigAxis_TA_ClosedForm(Tn)
    
    elem_out = Vector{Array{eltype(sig_m.elem[2])}}(undef, k+1)

    # ─────────────────────────────
    # Level 0: copy directly
    # ─────────────────────────────
    elem_out[1] = sig_m.elem[1]

    for j in 1:k
        σ_m = sig_m.elem[j+1]
        σ_n = sig_n.elem[j+1]

        tensor_j = Array{eltype(σ_m)}(undef, ntuple(_ -> m*n, j)...)

        for idx in CartesianIndices(tensor_j)
            linear_idx = Tuple(idx)

            idx_m = ntuple(t -> (div(linear_idx[t]-1, n) % m) + 1, j)
            idx_n = ntuple(t -> (mod(linear_idx[t]-1, n) + 1), j)

            tensor_j[idx] = σ_m[idx_m...] * σ_n[idx_n...]
        end

        elem_out[j+1] = tensor_j
    end

    return TruncatedTensorAlgebraElem{R, eltype(sig_m.elem[2])}(T, elem_out)
end


# -------------------------------
# seq_type == :p2
# -------------------------------
function sigAxis_p2_ClosedForm(T::TruncatedTensorAlgebra{R}, m::Int, n::Int) where R
    k = truncation_level(T)
    d=base_dimension(T)

    if m * n != d
        error("m * n != d")
    end


    Tm = TruncatedTensorAlgebra(base_algebra(T), m, k, :iis)
    Tn = TruncatedTensorAlgebra(base_algebra(T), n, k, :iis)

    sig_m = sigAxis_TA_ClosedForm(Tm)
    sig_n = sigAxis_TA_ClosedForm(Tn)

    elem_out = Vector{Any}(undef, k+1)

    # ─────────────────────────────
    # Level 0: copy directly
    # ─────────────────────────────
    elem_out[1] = sig_m.elem[1]

    # ─────────────────────────────
    # Level ≥ 1
    # ─────────────────────────────
    for j in 1:(k)
        σ_m = sig_m.elem[j+1]
        σ_n = sig_n.elem[j+1]

        perms = permutations_1_to_j(1:j)

        tensor_j = Array{eltype(σ_m)}(
            undef,
            (ntuple(_ -> m*n, j)..., factorial(j))
        )

        total = (m*n)^j
        perm_idx = 1

        for perm in perms
            for lin in 1:total
                # multi-índice in {1,…,m*n}^j
                linear_idx = ntuple(
                    t -> ((lin - 1) ÷ (m*n)^(t-1)) % (m*n) + 1,
                    j
                )

                # separate index m / n
                idx_m = ntuple(
                    t -> (div(linear_idx[t] - 1, n) % m) + 1,
                    j
                )

                idx_n = ntuple(
                    t -> (mod(linear_idx[t] - 1, n) + 1),
                    j
                )

                idx_n_perm = ntuple(t -> idx_n[perm[t]], j)

                tensor_j[linear_idx..., perm_idx] =
                    σ_m[idx_m...] * σ_n[idx_n_perm...]
            end
            perm_idx += 1
        end

        elem_out[j+1] = tensor_j
    end
    
    return TruncatedTensorAlgebraElem{R, eltype(sig_m.elem[2])}(T, elem_out)
end



# -------------------------------
# Principal function (dispatch)
# -------------------------------
function sigAxis_Chen(T::TruncatedTensorAlgebra{R}, m::Int, n::Int) where R
    if m * n != base_dimension(T)
        error("m * n != d")
    end

    if T.sequence_type == :p2id
        return sigAxis_p2id_Chen(T, m, n)
    elseif T.sequence_type == :p2
        return sigAxis_p2_Chen(T, m, n)
    else
        error("sequence_type must be :p2id or :p2")
    end
end


# -------------------------------
# seq_type == :p2id
# -------------------------------
function sigAxis_p2id_Chen(T::TruncatedTensorAlgebra{R}, m::Int, n::Int) where R
    k = truncation_level(T)
    
    if m * n != d
        error("m * n != d")
    end

    Tm = TruncatedTensorAlgebra(base_algebra(T), m, k, :iis)
    Tn = TruncatedTensorAlgebra(base_algebra(T), n, k, :iis)

    
    sig_m = sig_axis_TA(Tm)
    sig_n = sig_axis_TA(Tn)

    elem_out = Vector{Array{eltype(sig_m.elem[2])}}(undef, k+1)

    # ─────────────────────────────
    # Level 0: copy directly
    # ─────────────────────────────
    elem_out[1] = sig_m.elem[1]

    for j in 1:k
        σ_m = sig_m.elem[j+1]
        σ_n = sig_n.elem[j+1]

        tensor_j = Array{eltype(σ_m)}(undef, ntuple(_ -> m*n, j)...)

        for idx in CartesianIndices(tensor_j)
            linear_idx = Tuple(idx)

            idx_m = ntuple(t -> (div(linear_idx[t]-1, n) % m) + 1, j)
            idx_n = ntuple(t -> (mod(linear_idx[t]-1, n) + 1), j)

            tensor_j[idx] = σ_m[idx_m...] * σ_n[idx_n...]
        end

        elem_out[j] = tensor_j
    end

    return TruncatedTensorAlgebraElem{R, eltype(sig_m.elem[2])}(T, elem_out)
end


# -------------------------------
# seq_type == :p2
# -------------------------------
function sigAxis_p2_Chen(T::TruncatedTensorAlgebra{R}, m::Int, n::Int) where R
    k = truncation_level(T)
    d=base_dimension(T)

    if m * n != d
        error("m * n != d")
    end


    Tm = TruncatedTensorAlgebra(base_algebra(T), m, k, :iis)
    Tn = TruncatedTensorAlgebra(base_algebra(T), n, k, :iis)

    sig_m = sig_axis_TA(Tm)
    sig_n = sig_axis_TA(Tn)

    elem_out = Vector{Any}(undef, k+1)

    # ─────────────────────────────
    # Level 0: copy directly
    # ─────────────────────────────
    elem_out[1] = sig_m.elem[1]

    # ─────────────────────────────
    # Level ≥ 1
    # ─────────────────────────────
    for j in 1:(k)
        σ_m = sig_m.elem[j+1]
        σ_n = sig_n.elem[j+1]

        perms = permutations_1_to_j(1:j)

        tensor_j = Array{eltype(σ_m)}(
            undef,
            (ntuple(_ -> m*n, j)..., factorial(j))
        )

        total = (m*n)^j
        perm_idx = 1

        for perm in perms
            for lin in 1:total
                # multi-índice in {1,…,m*n}^j
                linear_idx = ntuple(
                    t -> ((lin - 1) ÷ (m*n)^(t-1)) % (m*n) + 1,
                    j
                )

                # separate index m / n
                idx_m = ntuple(
                    t -> (div(linear_idx[t] - 1, n) % m) + 1,
                    j
                )

                idx_n = ntuple(
                    t -> (mod(linear_idx[t] - 1, n) + 1),
                    j
                )

                idx_n_perm = ntuple(t -> idx_n[perm[t]], j)

                tensor_j[linear_idx..., perm_idx] =
                    σ_m[idx_m...] * σ_n[idx_n_perm...]
            end
            perm_idx += 1
        end

        elem_out[j+1] = tensor_j
    end
    
    return TruncatedTensorAlgebraElem{R, eltype(sig_m.elem[2])}(T, elem_out)
end


# ==============================================================
# Moment path at level j (shared helper)
# ==============================================================
"""
    moment_path_level(R::QQMPolyRing, d::Int, j::Int)

Constructs the moment tensor of level `j` for dimension `d` over the ring `R`.
Each multi-index entry is numerator / denominator * one(R), 
with numerator = prod(indices) and denominator = prod(cumsum(indices)).
"""
function moment_path_level(R::QQMPolyRing, d::Int, j::Int)
    dims = ntuple(_ -> d, j)
    T = Array{QQMPolyRingElem}(undef, dims...)

    for idx in Iterators.product(ntuple(_ -> 1:d, j)...)
        idx_vec = collect(idx)
        numerator = prod(idx_vec)
        denominator = prod(cumsum(idx_vec))
        T[idx...] = QQ(numerator, denominator) * one(R)
    end

    return T
end

# ==============================================================
# Moment membrane for :p2id
# ==============================================================
"""
    moment_membrane_p2id(TTA::TruncatedTensorAlgebra, m::Int, n::Int)

Compute the moment tensors for a TruncatedTensorAlgebra of sequence type `:p2id`.
Returns a new TTA object with `elem[j]` filled for all levels up to `truncation_level`.
"""
function moment_membrane_p2id(TTA::TruncatedTensorAlgebra, m::Int, n::Int)
    R = TTA.base_algebra
    k = TTA.truncation_level
    d = TTA.base_dimension

    if m * n != d
        error("m * n must equal the ambient dimension of TTA")
    end

    TTA_out = TruncatedTensorAlgebra(
        R, d, k, :p2id
    )
    TTA_out.elem = Vector{Any}(undef, k)

    for j in 1:k
        σ_m = moment_path_level(R, m, j)
        σ_n = moment_path_level(R, n, j)

        dims = ntuple(_ -> m*n, j)
        tensor_j = Array{QQMPolyRingElem}(undef, dims...)

        # Fill tensor with direct product of σ_m and σ_n
        for idx in Iterators.product(ntuple(_ -> 1:(m*n), j)...)
            idx_m = [(div(i-1, n) % m) + 1 for i in idx]
            idx_n = [(mod(i-1, n) + 1) for i in idx]
            tensor_j[idx...] = σ_m[idx_m...] * σ_n[idx_n...]
        end

        TTA_out.elem[j] = tensor_j
    end

    return TTA_out
end

# ==============================================================
# Moment membrane for :p2
# ==============================================================
"""
    moment_membrane_p2(TTA::TruncatedTensorAlgebra, m::Int, n::Int)

Compute the moment tensors for a TruncatedTensorAlgebra of sequence type `:p2`.
Returns a new TTA object with `elem[j]` filled for all levels up to `truncation_level`.
"""
function moment_membrane_p2(TTA::TruncatedTensorAlgebra, m::Int, n::Int)
    R = TTA.base_algebra
    k = TTA.truncation_level
    d = TTA.ambient_dimension

    if m * n != d
        error("m * n must equal the ambient dimension of TTA")
    end

    TTA_out = TruncatedTensorAlgebra(
        R, d, k, :p2
    )
    TTA_out.elem = Vector{Any}(undef, k)

    for j in 1:k
        σ_m = moment_path_level(R, m, j)
        σ_n = moment_path_level(R, n, j)

        perms = collect(permutations_1_to_j(j))
        tensor_j = Array{QQMPolyRingElem}(undef, (ntuple(_ -> m*n, j)..., factorial(j)))

        # Loop over permutations
        for (perm_idx, perm) in enumerate(perms)
            σ_n_perm = Array{QQMPolyRingElem}(undef, size(σ_n)...)
            for idx in Iterators.product(ntuple(_ -> 1:n, j)...)
                idx_perm = idx[perm]
                σ_n_perm[idx...] = σ_n[idx_perm...]
            end

            for idx in Iterators.product(ntuple(_ -> 1:(m*n), j)...)
                idx_m = [(div(i-1, n) % m) + 1 for i in idx]
                idx_n = [(mod(i-1, n) + 1) for i in idx]
                tensor_j[idx..., perm_idx] = σ_m[idx_m...] * σ_n_perm[idx_n...]
            end
        end

        TTA_out.elem[j] = tensor_j
    end

    return TTA_out
end

# ==============================================================
# Optional dispatch by seq_type
# ==============================================================
function moment_membrane(TTA::TruncatedTensorAlgebra, m::Int, n::Int)
    if TTA.sequence_type == :p2id
        return moment_membrane_p2id(TTA, m, n)
    elseif TTA.sequence_type == :p2
        return moment_membrane_p2(TTA, m, n)
    else
        error("sequence_type must be :p2id or :p2")
    end
end



"""
    applyMatrixToTTA(A::AbstractMatrix, X::TruncatedTensorAlgebra)

Applies a linear transformation `A` to all levels of a truncated tensor algebra `X`
by mode-product along all modes. Handles both `:p2id` and `:p2` signatures.

- `A` : matrix of size (d_new × d), linear transformation
- `X` : input truncated tensor algebra
Returns a **new** TruncatedTensorAlgebra with updated ambient dimension `d_new`.
"""
function applyMatrixToTTA(A::AbstractMatrix, X::TruncatedTensorAlgebra)
    d_new, d = size(A)
    d == X.base_dimension || error("size(A,2) must match ambient dimension")

    R = X.base_algebra
    k = X.truncation_level
    seq_type = X.sequence_type

    # Create new TTA with updated ambient dimension
    X_new = TruncatedTensorAlgebra(R, d_new, k, sequence_type=seq_type)

    for j in 1:k
        T = X.tensor_sequence[j]

        # Level 0 (scalar one) remains the same
        if T === one(R)
            X_new.tensor_sequence[j] = T
            continue
        end

        if seq_type == :p2id
            # For p2id, apply mode-product along all modes
            for mode in 1:j
                T = mode_product(T, A, mode, R)
            end
        elseif seq_type == :p2
            # For p2, we must account for permutations
            perms = permutations_1_to_j(j)
            T_perm = Array{typeof(one(R))}(undef, (ntuple(_ -> d_new, j)..., factorial(j)))

            for (perm_idx, perm) in enumerate(perms)
                T_temp = similar(T)
                # permute indices
                for idx in Iterators.product(ntuple(_ -> 1:size(T,1), j)...)
                    idx_perm = idx[perm]
                    T_temp[idx...] = T[idx_perm...]
                end
                # apply mode-product along all modes
                for mode in 1:j
                    T_temp = mode_product(T_temp, A, mode, R)
                end
                T_perm[:, :, :, perm_idx] = T_temp  # last axis stores permutations
            end
            T = T_perm
        else
            error("sequence_type must be :p2 or :p2id")
        end

        X_new.tensor_sequence[j] = T
    end

    return X_new
end


"""
    sig2parPoly(T::TruncatedTensorAlgebra, a::AbstractMatrix)

Given a truncated tensor algebra `T` (seq_type `:p2id` or `:p2`), and a matrix of coefficients `a`,
constructs the **moment membrane tensor** in the algebra of `T`, and then applies the matrix `a`
via matrix-tensor congruence. Returns a **new** truncated tensor algebra with updated ambient dimension.

- `T` : input truncated tensor algebra
- `a` : coefficient matrix (d_new × d)
"""
function sig2parPoly(T::TruncatedTensorAlgebra, a::AbstractMatrix)
    d_new, d = size(a)
    d == T.base_dimension || error("size(a,2) must match T.ambient_dimension")

    R = T.base_algebra
    k = T.truncation_level
    seq_type = T.sequence_type

    # Step 1: Build moment membrane tensor in same algebra as T
    T_moment = momentMembraneTensor(T)   # <-- returns TruncatedTensorAlgebra with same algebra

    # Step 2: Apply matrix-tensor congruence using a
    T_new = applyMatrixToTTA(a, T_moment)

    return T_new
end
