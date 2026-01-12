export 
  trunc_tensor_seq,
  matrix_tensorSeq_congruence, 
  sig_mono, 
  sig_axis, 
  sig_segment_standard_direction, 
  sig_segment, 
  bary_2samples, 
  bary_2nd_trunc_closedform, 
  bary_2nd_trunc, 
  bary, 
  bary_Nminus1_samples_fixed_inverse, 
  tensor_sequence, 
  ideal_of_entries, 
  embedding_matrix, # should not be exported
  coreSplineTrafo, # should not be exported 
  sig_pw_mono, 
  generate_lyndon_words, # should not be here
  dim_free_nil_lie_alg, # should not be here
  ideal_of_lyndon_entries
  
"""
    TruncTensorSeq(base_ring::QQMPolyRing, trunc_level::Int, amb_dim::Int)

Represents a truncated tensor sequence over a polynomial ring `QQMPolyRing`, with truncation
level `trunc_level` and ambient dimension `amb_dim`. Used to compute signatures of paths.

# Arguments
- `base_ring::QQMPolyRing`: the polynomial ring (e.g., `QQ[a11,a21,a12,a22]`)
- `trunc_level::Int`: truncation level (maximum tensor order)
- `amb_dim::Int`: ambient dimension

# Example
  >julia
   d, m, k = 2, 2, 2
   R, a = polynomial_ring_sig_transform(d, m)
   TTS = trunc_tensor_seq(R, k, m)
"""
TruncTensorSeq

struct TruncTensorSeq
  base_ring::QQMPolyRing
  trunc_level::Int
  amb_dim::Int
end

truncation_level(F::TruncTensorSeq) = F.trunc_level
ambient_dimension(F::TruncTensorSeq) = F.amb_dim
base_ring(F::TruncTensorSeq) = F.base_ring

"""
    TruncTensorSeqElem(parent::TruncTensorSeq, elem::Vector{Array{QQMPolyRingElem}})

Represents a single element of a truncated tensor sequence, corresponding to the
signature of a specific path.

# Arguments
- `parent::TruncTensorSeq`: the parent truncated tensor sequence
- `elem::Vector{Array{QQMPolyRingElem}}`: vector of arrays holding the tensor coefficients

# Example
    d, m, k = 2, 2, 2
    R, a = polynomial_ring_sig_transform(d, m)
    TTS = trunc_tensor_seq(R, k, m)
    seq = sig_axis(TTS)
    TTS_elem = TruncTensorSeqElem(TTS, seq)

# Methods
- `tensor_sequence(TTS_elem)`: returns the tensor data (vector of arrays)
- `Base.parent(TTS_elem)`: returns the parent `TruncTensorSeq`
"""
TruncTensorSeqElem

struct TruncTensorSeqElem
  parent::TruncTensorSeq
  elem::Vector{Array{QQMPolyRingElem}}
end

Base.parent(a::TruncTensorSeqElem) = a.parent
tensor_sequence(a::TruncTensorSeqElem) = a.elem

"""
    trunc_tensor_seq(base_ring::QQMPolyRing, trunc_level::Int, amb_dim::Int)

Convenience constructor for `TruncTensorSeq`. Creates a truncated tensor sequence
with the given base polynomial ring, truncation level, and ambient dimension.

# Arguments
- `base_ring::QQMPolyRing`: the polynomial ring
- `trunc_level::Int`: truncation level
- `amb_dim::Int`: ambient dimension

# Example
    R, a = polynomial_ring_sig_transform(2,2)
    TTS = trunc_tensor_seq(R, 2, 2)
"""
trunc_tensor_seq

function trunc_tensor_seq(base_ring::QQMPolyRing,trunc_level::Int,amb_dim::Int)
  TTS = TruncTensorSeq(base_ring, trunc_level, amb_dim)
  return TTS
end


"""
    trunc_tensor_seq_elem(TTS::TruncTensorSeq, seq::Vector{Array{QQMPolyRingElem}})

Constructs a `TruncTensorSeqElem` from a truncated tensor sequence `TTS` and a sequence of polynomial arrays `seq`.
"""
trunc_tensor_seq_elem

function trunc_tensor_seq_elem(TTS::TruncTensorSeq,seq::Vector{Array{QQMPolyRingElem}})
  return TruncTensorSeqElem(TTS, seq)
end

##################
# Constructor helpers
##################

function _C0(_k::Int, _order::Int, _R::QQMPolyRing, is_one::Bool)
   res = Array{QQMPolyRingElem}(undef, _order*ones(Int,_k)...)
   for idx in CartesianIndices(res)  
       if is_one && _k ==0 
         res[idx] = one(_R)
       else
         res[idx] = zero(_R)
       end 
   end
   return res
end

function _C0_seq(_trunc_level::Int, _order::Int, _R::QQMPolyRing, is_one::Bool)
   C0 = [_C0(_k,_order, _R, is_one) for _k in (0:_trunc_level)]
   return C0
end

function _Cmono(_k::Int, _order::Int, _R::QQMPolyRing)
   res = Array{QQMPolyRingElem}(undef, _order*ones(Int,_k)...)
   for idx in CartesianIndices(res)
       if _k ==0 
         res[idx] = one(_R)
       else
         res[idx] = QQ(prod(Tuple(idx)),prod(cumsum(Tuple(idx))))*one(_R)
       end 
   end
   return res
end

function _Cmono_seq(_trunc_level::Int, _order::Int, _R::QQMPolyRing)
   Cmono = [_Cmono(_k,_order, _R) for _k in (0:_trunc_level)]
   return Cmono
end

function _one_hot(_i::Int,_n::Int,_R::QQMPolyRing)
  res = zeros(_R,_n)
  res[_i] = one(_R)
  return res
end

#function sig_segment_standard_direction(_i::Int,_order::Int,_trunc_level::Int,_R::QQMPolyRing)
#   return sig_segment(_one_hot(_i,_order,_R),_trunc_level)
#end
#
#function sig_segment(v::Vector{QQMPolyRingElem},_trunc_level::Int)
#   n = length(v)
#   C = Cmono_seq(_trunc_level,1,parent(v[1]))
#   return matrix_tensorSeq_congruence(v,C)
#end

function Caxis_seq(_trunc_level::Int, _order::Int, _R::QQMPolyRing)
   sample = [sig_segment_standard_direction(i,_order,_trunc_level,_R) for i in (1:_order)]
   F,s = free_trunc_sig_alg_multiv(_trunc_level,_order)
   chen = trunc_multivariate(prod([free_sig_elem_from_sample(i,s) for i in (1:_order)]),_trunc_level)
   return eval_to_tensorSeq_multivariate(chen,sample)
end

##################
# Constructors
##################

function Base.zero(TTS::TruncTensorSeq)
  k = truncation_level(TTS) 
  d = ambient_dimension(TTS)
  R = base_ring(TTS)
  return TruncTensorSeqElem(TTS,_C0_seq(k,d,R,false))
end

function Base.one(TTS::TruncTensorSeq)
  k = truncation_level(TTS) 
  d = ambient_dimension(TTS)
  R = base_ring(TTS)
  return TruncTensorSeqElem(TTS,_C0_seq(k,d,R,true))  ### this is not working 
end

"""
    sig_mono(TTS::TruncTensorSeq)

Constructs a `TruncTensorSeqElem` representing the **monomial signature sequence**
for a given truncated tensor sequence `TTS`.  

This function builds the canonical **monomial sequence** of tensors used in the 
computation of signatures for paths in the algebraic setting.

# Arguments
- `TTS::TruncTensorSeq`: the parent truncated tensor sequence.

# Returns
- `TruncTensorSeqElem`: a single element containing the monomial tensor sequence.

# Example
d, m, k = 2, 2, 2
R, a = polynomial_ring_sig_transform(d, m)
TTS = trunc_tensor_seq(R, k, m)
mono_seq = sig_mono(TTS)
"""
sig_mono

function sig_mono(TTS::TruncTensorSeq)
  k = truncation_level(TTS) 
  d = ambient_dimension(TTS)
  R = base_ring(TTS)
  return TruncTensorSeqElem(TTS,_Cmono_seq(k,d,R))
end

function sig_segment(TTS::TruncTensorSeq,v::Vector{QQMPolyRingElem})
  k = truncation_level(TTS) 
  d = ambient_dimension(TTS)
  R = base_ring(TTS)
  n = length(v)
  @assert n == d "dimensions do not match"
  @assert R == parent(v[1]) "rings do not match"
  TTS1 = trunc_tensor_seq(R,k,1)
  C = sig_mono(TTS1)
  #C = sig_mono(_trunc_level,1,parent(v[1]))
  return matrix_tensorSeq_congruence(v,C)
end

"""
    sig_segment_standard_direction(TTS::TruncTensorSeq, j::Int)

Returns the truncated signature of the standard segment `Ej` according to Theorem 6.7.

# Example
    R, a = polynomial_ring_sig_transform(2,2)
    TTS = trunc_tensor_seq(R, 2, 2)
    sigEj = sig_segment_standard_direction(TTS, 1)
"""
sig_segment_standard_direction


function sig_segment_standard_direction(TTS::TruncTensorSeq,_i::Int)
  d = ambient_dimension(TTS)
  R = base_ring(TTS)
  return sig_segment(TTS,_one_hot(_i,d,R))
end

"""
    sig_axis(TTS::TruncTensorSeq)

Returns the tensor sequence axes of a truncated tensor sequence `TTS`.
Used as input for `matrix_tensorSeq_congruence`.

# Example
    R, a = polynomial_ring_sig_transform(2,2)
    TTS = trunc_tensor_seq(R, 2, 2)
    axes = sig_axis(TTS)
"""
sig_axis

function sig_axis(TTS::TruncTensorSeq)
  k = truncation_level(TTS) 
  d = ambient_dimension(TTS)
  R = base_ring(TTS)
  sample = [sig_segment_standard_direction(TTS,i) for i in (1:d)]
  F,s = free_trunc_sig_alg_multiv(k,d)
  chen = prod([free_sig_from_sample(i,F) for i in (1:d)])
  return evaluate(chen,sample) #TODO: this should be the same as multiplication of the signatuers directly, without fromal chen and eval 
end

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


"""
    sig_pw_mono(TTS::TruncTensorSeq, m)

Constructs the **piecewise monomial signature** of a truncated tensor sequence `TTS` 
according to a composition `m` of the ambient dimension.

This function computes the tensor sequence signature by splitting the ambient dimension 
into blocks specified by `m` and applying `matrix_tensorSeq_congruence` and `sig_mono` 
to each block. The results are then multiplied together to obtain the piecewise signature.

# Arguments
- `TTS::TruncTensorSeq`: the parent truncated tensor sequence.
- `m::Vector{Int}`: a composition of the ambient dimension (vector of positive integers summing to `ambient_dimension(TTS)`).

# Returns
- `TruncTensorSeqElem`: the piecewise monomial signature element corresponding to the composition `m`.

# Example
d, k = 4, 2
R, vars = polynomial_ring_sig_transform(d,2)
TTS = trunc_tensor_seq(R, k, d)
m = [2,2]  # composition of the ambient dimension
pw_mono_sig = sig_pw_mono(TTS, m)
"""
sig_pw_mono

function sig_pw_mono(TTS::TruncTensorSeq,m,r=0)
  k = truncation_level(TTS) 
  R = base_ring(TTS)
  d = ambient_dimension(TTS)
  if d != sum(m) - r*(length(m) -1)
    error("m must be a composition of the ambient dimension") 
  end 
  sigs = [matrix_tensorSeq_congruence(Array(embedding_matrix(m,i)),sig_mono(trunc_tensor_seq(R,k,m[i]))) for i in (1:length(m))]
  res = prod(sigs)
  if r == 0 
    return res
  else 
    return matrix_tensorSeq_congruence(Array(coreSplineTrafo(m,r)),res)
  end
end

function Base.Matrix(a::TruncTensorSeqElem)
  return tensor_sequence(a)[3]
end 

function Oscar.QQMatrix(a::TruncTensorSeqElem) 
  R = base_ring(parent(a))
  M = tensor_sequence(a)[3]
  a
  z = size(M)
  res = zero_matrix(QQ, z[1], z[2])
  for i in (1:z[1])
    for j in (1:z[2])
      if M[i,j] != zero(QQ)
        res[i,j] = leading_coefficient(M[i,j])
      end
    end 
  end 
  return res
end 

#################
# ideal constructors 
################# 

"""
    ideal_of_entries(mat::AbstractMatrix)

Returns the ideal in the polynomial ring generated by the entries of `mat`.

# Example
    I = ideal_of_entries(s_pwln - sY)
"""
ideal_of_entries

function ideal_of_entries(a::TruncTensorSeqElem)
  A = parent(a)
  R = base_ring(A)
  k = truncation_level(A)
  temp = (tensor_sequence(a))[2:k+1]
  gens4learn = vcat([vec(temp[i]) for i in (1:length(temp))]...)
  return  ideal(R,gens4learn)
end

function nextword(k::Int, w::Vector{Int}, alphabet::Vector{Int})
    # repeat w enough times and truncate to length n
    reps = (k ÷ length(w)) + 1
    x = repeat(w, reps)[1:k]
    # remove trailing maximal letters
    while !isempty(x) && x[end] == alphabet[end]
        pop!(x)
    end
    if !isempty(x)
        last_char = x[end]
        next_char_index = findfirst(==(last_char), alphabet) + 1
        x[end] = alphabet[next_char_index]
    end
    return x
end


function generate_lyndon_words(k::Int, alphabet::Vector{Int})
    lwords = Vector{Vector{Int}}()
    w = [alphabet[1]]
    while length(w) <= k
        push!(lwords, copy(w))
        w = nextword(k, w, alphabet)
        isempty(w) && break
    end
    return lwords
end

function dim_free_nil_lie_alg(d::Int,k::Int)
  return sum([div(moebius_mu(a)*d^(div(ell,a)),ell) for ell in (1:k) for a in (1:k) if divides(ell,a)[1]])
end

function ideal_of_lyndon_entries(tts::TruncTensorSeqElem)
  A = parent(tts)
  R = base_ring(A)
  d = ambient_dimension(A)
  k = truncation_level(A)
  lynd = generate_lyndon_words(k, Vector((1:d)))
  res = ideal(R,zero(R))
  elem_tts = tensor_sequence(tts)
  return ideal(R,[elem_tts[length(w)+1][w...] for w in lynd])
end



##################
# printing
##################

function Base.show(io::IO, ::MIME"text/plain",f::TruncTensorSeqElem)
  A = parent(f)
  seq = f.elem
  trunc_level = truncation_level(A)
  for i in (1:trunc_level)
    show(io, "text/plain", seq[i])
    println()
    println("⊕")
  end
  show(io, "text/plain", seq[trunc_level+1])
end

##################
# evaluation 
##################

function Oscar.evaluate(p::FreeTruncSigAlgMultivElem, bs::Vector{TruncTensorSeqElem})
  A = parent(bs[1])
  k = truncation_level(A)
  res_seq = [evaluate_to_tensor(p,bs,i) for i in (0:k)]
  return trunc_tensor_seq_elem(A,res_seq)
end

function leading_coefficient_and_zero(p)
  if is_zero(p)
    return zero(QQ)
  else 
    return leading_coefficient(p)
  end 
end 

function evaluate_to_tensor(p::FreeTruncSigAlgMultivElem,  _bs::Vector{TruncTensorSeqElem}, index::Int)
   bs = deepcopy(_bs)
   if index == 0
     res = tensor_sequence(bs[1])[1]
     res[] = leading_coefficient_and_zero(FreeTruncSigAlgMultivElem_to_polynomial(graded_component(p,0)))*one(res[])
     return res
   end 
   num_vars = length(bs)
   A = parent(bs[1])
   k = truncation_level(A)
   fi = FreeTruncSigAlgMultivElem_to_polynomial(graded_component(p,index))
   seq_tensor_seq = [tensor_sequence(b)[2:k+1] for b in bs]
   res = 0*seq_tensor_seq[1][index] # zero tensor of correct size
   for t in (1:length(fi))
       fit = collect(terms(fi))[t]
       first_exp = fit.exps[1][1]
       sample_var = div(first_exp-1,k)+1
       tensor_var = mod(first_exp-1,k)+1
       temp = leading_coefficient_and_zero(fit).*seq_tensor_seq[sample_var][tensor_var]
       for fak in (2:length(fit.exps[1]))
          temp_exp = fit.exps[1][fak]
          sample_var = div(temp_exp-1,k)+1
          tensor_var = mod(temp_exp-1,k)+1
          temp = concatenate_tensors(temp,seq_tensor_seq[sample_var][tensor_var])
       end
       res += temp
   end
   return res
end


##################
# arithmetic tensor sequences
##################

function Base.:(==)(a::TruncTensorSeqElem, b::TruncTensorSeqElem)
    return parent(a) == parent(b) && tensor_sequence(a) == tensor_sequence(b)
end

function Base.:+(a::TruncTensorSeqElem,
                 b::TruncTensorSeqElem)
  A = parent(a)
  res_seq = tensor_sequence(a) + tensor_sequence(b)
  return trunc_tensor_seq_elem(A,res_seq)
end


function Base.:-(a::TruncTensorSeqElem,
                 b::TruncTensorSeqElem)
  A = parent(a)
  res_seq = tensor_sequence(a) - tensor_sequence(b)
  return trunc_tensor_seq_elem(A,res_seq)
end

function Base.:*(a::FieldElem,
                 b::TruncTensorSeqElem)
  A = parent(b)
  k = truncation_level(A)
  res_seq = tensor_sequence(b)
  res_seq[1][] = a*res_seq[1][]
  for j in (2:k+1) 
    res_seq[j] = a.*res_seq[j]
  end 
  return trunc_tensor_seq_elem(A,res_seq)
end

function Base.:*(a::TruncTensorSeqElem, 
                 b::TruncTensorSeqElem)
  A = parent(a)
  k = truncation_level(A)
  F,s = free_trunc_sig_alg_multiv(k,2)
  chen = prod([free_sig_from_sample(i,F) for i in (1:2)]) 
  return evaluate(chen,[a,b])
end

function Base.:inv(a::TruncTensorSeqElem)
  A = parent(a)
  k = truncation_level(A)
  F,s = free_trunc_sig_alg_multiv(k,1)
  p = inv(free_sig_from_sample(1,F))
  return evaluate(p,[a])
end

function Base.:^(a::TruncTensorSeqElem, 
                 n::Int)
  A = parent(a)
  #k = truncation_level(A)
  res = one(A)
  if n >= 0
    temp = a
  else 
    temp = inv(a)
  end 
  for i in (1:abs(n))  
    res = res * temp
  end 
  return res
end

function Base.:exp(a::TruncTensorSeqElem)
  A = parent(a)
  k = truncation_level(A)
  F,s = free_trunc_sig_alg_multiv(k,1)
  p = exp(free_sig_from_sample(1,F)-one(F))
  return evaluate(p,[a])
end

function Base.:log(a::TruncTensorSeqElem)
  A = parent(a)
  k = truncation_level(A)
  F,s = free_trunc_sig_alg_multiv(k,1)
  p = log(free_sig_from_sample(1,F))
  return evaluate(p,[a])
end

"""
    matrix_tensorSeq_congruence(matrix::Array, b::TruncTensorSeqElem)

Returns a `TruncTensorSeqElem` representing the tensor congruence of the vector/matrix `a`
with respect to the tensor axes.

# Example
    R, a = polynomial_ring_sig_transform(2,2)
    TTS = trunc_tensor_seq(R, 2, 2)
    sigX = matrix_tensorSeq_congruence([1, 2], sig_axis(TTS))
"""
matrix_tensorSeq_congruence 

function matrix_tensorSeq_congruence(matrix::Array, b::TruncTensorSeqElem)
  TTSm = parent(b) 
  k = truncation_level(TTSm)
  m = ambient_dimension(TTSm) 
  R = base_ring(TTSm)
  @assert m == size(matrix, 2) "Matrix and tensor dimensions do not align!"
  d = size(matrix, 1)
  TTSd = trunc_tensor_seq(R,k,d)
  tensorSeq = tensor_sequence(b)
  resSeq = [matrix_tensor_congruence(matrix, tensor) for tensor in tensorSeq]
  return trunc_tensor_seq_elem(TTSd,resSeq)
end

"""
    bary_2samples(a::TruncTensorSeqElem, b::TruncTensorSeqElem)

Computes the **barycenter (geometric mean)** of two truncated tensor sequence elements `a` and `b` 
in the free truncated signature algebra.  

This is done using the **log–exp formula** in the algebra, which corresponds to the 2-point 
barycenter in the Lie group of truncated signatures:


    bary(a,b) = a cdot exp Big(frac{1}{2} log(a^{-1} cdot b)Big)


# Arguments
- `a::TruncTensorSeqElem`: first truncated tensor sequence element
- `b::TruncTensorSeqElem`: second truncated tensor sequence element

# Returns
- `TruncTensorSeqElem`: the barycenter element of `a` and `b`.

# Example

d, m, k = 2, 2, 2
R, vars = polynomial_ring_sig_transform(d, m)
TTS = trunc_tensor_seq(R, k, m)
x = sig_axis(TTS)[1]  # example first path
y = sig_axis(TTS)[2]  # example second path
bary = bary_2samples(x, y)
"""
bary_2samples

function bary_2samples(a::TruncTensorSeqElem, b::TruncTensorSeqElem)
  return a*exp(QQ(1,2)*log(inv(a)*b))
end


"""
    bary(signatures::Vector{TruncTensorSeqElem})

Computes the barycenter of a list of truncated tensor sequence elements.

# Example
    R, a = polynomial_ring_sig_transform(2,1)
    TTS = trunc_tensor_seq(R, 2, 1)
    sigX1 = matrix_tensorSeq_congruence([1//2, 1], sig_axis(TTS))
    sigX2 = matrix_tensorSeq_congruence([1, 1//2], sig_axis(TTS))
    sY = bary([sigX1, sigX2])
"""
bary

function bary(bs::Vector{TruncTensorSeqElem})
  TTSm = parent(bs[1])
  res = one(TTSm)
  k = truncation_level(TTSm)
  N = length(bs)
  p = bary_defining_polynomial_system(k,N) 
  s = gens_in_shape(parent(p))
  y = s[:,N+1]
  for j in (1:k)
    pj = graded_component(p,j)
    res_seq = tensor_sequence(res)
    bs_new = [bs;res]
    res_seq[j+1] = tensor_sequence(evaluate(pj+y[j],bs_new))[j+1]
    res = trunc_tensor_seq_elem(TTSm,res_seq) 
  end 
  return res 
end


"""
    bary_Nminus1_samples_fixed_inverse(bs::Vector{TruncTensorSeqElem}, y::TruncTensorSeqElem)

Computes the barycenter of `N` truncated tensor sequence elements where the **last element `y` is fixed**, 
by solving the polynomial system defining the barycenter in the free truncated signature algebra.  

This function iteratively computes each **graded component** of the barycenter sequence 
using the polynomials from `bary_defining_polynomial_system`, adjusting the components to satisfy 
the barycenter condition with the fixed inverse of the last element.

# Arguments
- `bs::Vector{TruncTensorSeqElem}`: vector of the first `N-1` elements in the barycenter computation.
- `y::TruncTensorSeqElem`: the `N`-th element which remains fixed during the computation.

# Returns
- `TruncTensorSeqElem`: the barycenter element of `[bs..., y]` satisfying the fixed-inverse constraint.

# Example
R, vars = polynomial_ring_sig_transform(2,2)
TTS = trunc_tensor_seq(R, 3, 2)
x1 = sig_axis(TTS)[1]
x2 = sig_axis(TTS)[2]
x3 = sig_axis(TTS)[3]
bary_fixed = bary_Nminus1_samples_fixed_inverse([x1, x2], x3)
"""
bary_Nminus1_samples_fixed_inverse

function bary_Nminus1_samples_fixed_inverse(bs::Vector{TruncTensorSeqElem},y::TruncTensorSeqElem)
  TTSm = parent(bs[1])
  res = one(TTSm)
  k = truncation_level(TTSm)
  N = length(bs) + 1
  p = bary_defining_polynomial_system(k,N) 
  s = gens_in_shape(parent(p))
  xN = s[:,N]
  for j in (1:k)
    pj = graded_component(p,j)
    res_seq = tensor_sequence(res)
    bs_new = [[bs[1:N-1];res];y]
    res_seq[j+1] = tensor_sequence(evaluate(-QQ(N,1)*pj+xN[j],bs_new))[j+1]
    res = trunc_tensor_seq_elem(TTSm,res_seq) 
  end 
  return res 
end


"""
    bary_2nd_trunc(bs::Vector{TruncTensorSeqElem})

Computes the **second-order barycenter** of a collection of truncated tensor sequence elements 
`bs = [b_1, ..., b_N]` in the free truncated signature algebra.  

The barycenter is computed in the **log–exp Lie algebra sense**, generalizing the 2-sample 
barycenter to `N` elements at truncation level 2:


bary_2(bs) = exp Big( frac{1}{N} sum_{i=1}^N log(b_i) Big)


# Arguments
- `bs::Vector{TruncTensorSeqElem}`: vector of truncated tensor sequence elements to average.

# Returns
- `TruncTensorSeqElem`: the second-order barycenter element.

# Example
R, vars = polynomial_ring_sig_transform(2,2)
TTS = trunc_tensor_seq(R, 2, 2)
x = sig_axis(TTS)[1]
y = sig_axis(TTS)[2]
z = sig_axis(TTS)[3]
bary = bary_2nd_trunc([x, y, z])
"""
bary_2nd_trunc

function bary_2nd_trunc(bs::Vector{TruncTensorSeqElem})
  N = length(bs)
  return exp(QQ(1,N)*sum(log.(bs))) 
end

"""
    bary_2nd_trunc_closedform(bs::Vector{TruncTensorSeqElem})

Computes the **second-order barycenter in closed form** for a collection of truncated tensor 
sequence elements `bs` at truncation level 2.  

This method avoids the log exp computation by explicitly computing the **first, second, and 
second-order tensor components**:


\text{bary}_2^\text{closed}(bs) = \bigg[ \text{level 1}, \text{level 2 vector}, \text{level 2 matrix} \bigg]


where the matrix is corrected to account for cross terms between first-order components.

# Arguments
- `bs::Vector{TruncTensorSeqElem}`: vector of truncated tensor sequence elements.

# Returns
- `TruncTensorSeqElem`: barycenter element containing the truncated sequence `[level1, level2_vector, level2_matrix]`.

# Example
R, vars = polynomial_ring_sig_transform(2,2)
TTS = trunc_tensor_seq(R, 2, 2)
x = sig_axis(TTS)[1]
y = sig_axis(TTS)[2]
z = sig_axis(TTS)[3]
bary_closed = bary_2nd_trunc_closedform([x, y, z])
"""
bary_2nd_trunc_closedform

function bary_2nd_trunc_closedform(bs::Vector{TruncTensorSeqElem})
  TTSm = parent(bs[1])
  N = length(bs)
  ls = tensor_sequence.(bs)
  vec = QQ(1,N).*sum(ls[i][2] for i in (1:N))
  mat = QQ(1,N).*sum(ls[i][3] for i in (1:N)) - QQ(1,2*N).*sum(ls[i][2]*transpose(ls[i][2]) for i in (1:N)) + QQ(1,2*N^2).*sum(ls[i1][2]*transpose(ls[i2][2]) for i1 in (1:N) for i2 in (1:N))
  return trunc_tensor_seq_elem(TTSm,[ls[1][1],vec,mat])
end



