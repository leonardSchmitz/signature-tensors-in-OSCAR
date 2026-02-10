export TruncatedTensorAlgebra,
       TruncatedTensorAlgebraElem,
       truncation_level,
       ambient_dimension,
       base_dimension,
       base_ring,
       base_algebra,
       tensor_sequence,
       zero,
       one,
       sig_mono_TA,     # soon removed from the export 
       sigAxis_TA_ClosedForm, # sr
       sig_axis_TA, # sr
       sig_poly_TA, # sr
       sig_pwln_TA_Congruence, # sr
       sig_pwl_TA_chen, # sr
       sig_segment_TA, # sr
       sig_segment_standard_direction_TA, # sr
       matrix_tensorAlg_congruence_TA,

       # Signature axis functions
       sigAxis_ClosedForm,  # sr
       sigAxis_p2id_ClosedForm,  # sr
       sigAxis_p2_ClosedForm,  # sr
       sigAxis_Chen,  # sr
       sigAxis_p2id_Chen,  # sr
       sigAxis_p2_Chen,  # sr

       # Moment functions
       moment_path_level,  # sr
       moment_membrane,  # sr
       moment_membrane_p2id,  # sr
       moment_membrane_p2,  # sr

      # Tensor algebra operations
       mode_product,      # soon be removed to tensor operations?
       applyMatrixToTTA, 
       sig2parPoly

    



struct TruncatedTensorAlgebra{R}
    base_algebra::R
    base_dimension::Int
    truncation_level::Int
    sequence_type::Symbol
end

struct TruncatedTensorAlgebraElem{R,E}
    parent::TruncatedTensorAlgebra{R}
    elem::Vector{Array{E}}
end

truncation_level(F::TruncatedTensorAlgebra) = F.truncation_level
ambient_dimension(F::TruncatedTensorAlgebra) = F.base_dimension
base_dimension(F::TruncatedTensorAlgebra) = F.base_dimension
base_ring(F::TruncatedTensorAlgebra) = F.base_algebra
base_algebra(F::TruncatedTensorAlgebra) = F.base_algebra

Base.parent(a::TruncatedTensorAlgebraElem) = a.parent
tensor_sequence(a::TruncatedTensorAlgebraElem) = a.elem


function TruncatedTensorAlgebra(R, d::Int, k::Int, seq::Symbol=:iis)

    d ≥ 0 || error("ambient dimension must be ≥ 0")
    k ≥ 0 || error("truncation level must be ≥ 0")

    if seq ∉ (:iis, :p2id, :p2)
        error("seq must be one of :iis, :p2id, :p2")
    end

    A = TruncatedTensorAlgebra{typeof(R)}(R, d, k, seq)

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
    R   = A.base_ring
    d   = A.amb_dim
    k   = A.trunc_level
    seq = A.sequence_type

    elems = Vector{Any}(undef, k+1)

    # =========================
    # Level 0 (common for all Algebra)
    # =========================
    elems[1] = fill(is_one ? one(R) : zero(R), ())

    # =========================
    # IIS
    # =========================
    if seq == :iis
        for n in 1:k
            elems[n+1] = zeros(R, ntuple(_ -> d, n)...)
        end
        return elems
    end

    # =========================
    # P2ID
    # =========================
    if seq == :p2id
        for n in 1:k
            dims = ntuple(_ -> d, n)
            elems[n+1] = Array{Any}(undef, dims...)
        end
        return elems
    end

    # =========================
    # P2
    # =========================
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

function Base.zero(T::TruncatedTensorAlgebra)

    k = truncation_level(T)
    d = base_dimension(T)
    R = base_algebra(T)   # or base_ring(T)
    seq = T.sequence_type

    elem = Vector{Any}(undef, k + 1)

    ####################
    # Level 0 (always)
    ####################
    elem[1] = fill(zero(R), ())   # 0-dimensional scalar

    ####################
    # Levels ≥ 1
    ####################
    if seq == :iis
        for n in 1:k
            elem[n + 1] = zeros(R, ntuple(_ -> d, n)...)
        end

    elseif seq == :p2id
        for n in 1:k
            if n == 1
                dims = (d,)
            else
                dims = ntuple(_ -> d, n)
            end
            elem[n + 1] = zeros(R, dims...)
        end

    elseif seq == :p2
        for n in 1:k
            if n == 1
                dims = (d,)
            else
                dims = (ntuple(_ -> d, n)..., factorial(n))
            end
            elem[n + 1] = zeros(R, dims...)
        end

    else
        throw(ArgumentError("zero(T) not implemented for sequence_type = $seq"))
    end

    return TruncatedTensorAlgebraElem(T, elem)
end




function Base.one(T::TruncatedTensorAlgebra{R}) where R
    if T.sequence_type == :iis

        k = truncation_level(T)
        d = ambient_dimension(T)
        alg = T.base_algebra
        E = typeof(one(alg))

        return TruncatedTensorAlgebraElem{R,E}(T, _C0_seq_TA(k, d, alg, true))
    else
        throw(ArgumentError("one(T) is only defined for sequence_type = :iis"))
    end
end

function Base.show(io::IO, x::TruncatedTensorAlgebraElem)
    for (i, t) in enumerate(x.elem)
        show(io, MIME("text/plain"), t)
        if i < length(x.elem)
            println(io, "\n⊕")
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



# 1) return 1 of the correct type
function one_elem(R)
    if R === QQField
        return QQ(1)          # element of field QQ
    elseif R <: Type
        return one(R)         # if R is of the element type
    else
        return one(R)         # if R is of the ring type
    end
end

# 2) return 0
function zero_elem(R)
    if R === QQField
        return QQ(0)
    elseif R <: Type
        return zero(R)
    else
        return zero(R)
    end
end


##################
#.  sig_mono 
##################

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
    res = fill(zero_elem(_R), _n)
    res[_i] = one_elem(_R)
    return res
end


function sig_mono_TA(T::TruncatedTensorAlgebra{R}) where R
    if T.sequence_type != :iis
        throw(ArgumentError("sig_mono only defined for sequence_type = :iis"))
    end

    k = truncation_level(T)
    d = ambient_dimension(T)

    # base algebra
    R0 = base_algebra(T)

    # seq must be an Vector{Array{E}} for some E
    seq = _Cmono_seqTA(k, d, R0)

    # determinate E using seq
    E = eltype(seq[1])   # type of the elements of the intern array
    TruncatedTensorAlgebraElem{R, E}(T, seq)


    return TruncatedTensorAlgebraElem{R, E}(T, seq)
end



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

#function Base.getindex(x::TruncatedTensorAlgebraElem, r::UnitRange{Int})
#    return x.elem[r .+ 1]
#end





# ----------------------------
# congruence for matrices
# ----------------------------
# ------------------------------------------------------------
# helper: obtain the father ring of an element 
# ------------------------------------------------------------
ring_of(::Type{T}) where T = parent(one(T))
ring_of(x) = ring_of(typeof(x))

# ------------------------------------------------------------
# common_ring: verify the rule
# ------------------------------------------------------------
function common_ring(R_tensor, R_matrix)
    if R_matrix isa MPolyRing
        return R_matrix
    else
        return R_tensor
    end
end



function matrix_tensorAlg_congruence_TA(
    matrix::AbstractMatrix,
    b::TruncatedTensorAlgebraElem{R,E}
) where {R,E}

    T = parent(b)
    k = truncation_level(T)
    m = ambient_dimension(T)
    @assert m == size(matrix, 2)
    d = size(matrix, 1)

    # Ring
    R_tensor = base_ring(T)

    R_matrix = parent(matrix[1,1])

    # Simple Rule 
    Rnew = common_ring(R_tensor, R_matrix)

    # New Algebra 
    Tnew = TruncatedTensorAlgebra(Rnew, d, k, :iis)

    # Convert tensor in the new Ring
    tensorSeq = tensor_sequence(b)
    tensorSeq = [map(x -> Rnew(x), t) for t in tensorSeq]

    # Apply congruence
    resSeq = [
        matrix_tensor_congruence_TA(map(x -> Rnew(x), matrix), t)
        for t in tensorSeq
    ]

    return TruncatedTensorAlgebraElem(Tnew, resSeq)
end


# ----------------------------
# congruence for vectors (treat as 1×d matrices)
# ----------------------------
function matrix_tensorAlg_congruence_TA(
    v::AbstractVector,
    b::TruncatedTensorAlgebraElem{R,E}
) where {R,E}
    M = reshape(v, 1, length(v))
    return matrix_tensorAlg_congruence_TA(M, b)
end

function sig_segment_TA(T::TruncatedTensorAlgebra{R}, v::Vector{E}) where {R,E}
    # --- validation---
    k = truncation_level(T)
    d = ambient_dimension(T)

    @assert length(v) == d "dimensions do not match"

    # --- veryfy the types---
    @assert typeof(v[1]) <: E "element types do not match"

    # --- create a new algebra ---
    T1 = TruncatedTensorAlgebra(base_ring(T), 1, k, :iis)

    # --- signature mono in dimension 1 ---
    C = sig_mono_TA(T1)

    # --- Convert the vector v in a column matrix d×1 ---
    v_mat = reshape(v, d, 1)

    # --- Apply the matricial congruence (column matrix* signature 1D) ---
    return matrix_tensorAlg_congruence_TA(v_mat, C)
end


function sig_segment_standard_direction_TA(T::TruncatedTensorAlgebra{R}, _i::Int) where {R}
    d = ambient_dimension(T)
    # vector one-hot in the direction _i
    v = _one_hot_TA(_i, d, R)
    return sig_segment_TA(T, v)
end



function Base.:*(a::R, b::TruncatedTensorAlgebraElem{R,E}) where {R,E}
    A = parent(b)
    k = truncation_level(A)

    res_seq = tensor_sequence(b)

    # first level 
    res_seq[1][] = a * res_seq[1][]

    # all other levels
    for j in 2:k+1
        res_seq[j] = a .* res_seq[j]
    end

    return TruncatedTensorAlgebraElem(A, res_seq)
end



function Base.:*(a::TruncatedTensorAlgebraElem{R,E}, 
                 b::TruncatedTensorAlgebraElem{R,E}) where {R,E}
    A = parent(a)
    k = truncation_level(A)

    F, s = free_trunc_sig_alg_multiv(k, 2)
    chen = prod([free_sig_from_sample(i, F) for i in 1:2])

    return evaluate(chen, [a, b])
end

function Base.:*(matrix::AbstractMatrix, x::TruncatedTensorAlgebraElem)
    T = parent(x)

    if T.sequence_type == :iis
            return matrix_tensorAlg_congruence_TA(matrix, x)
    else
                throw(ArgumentError("matrix * TruncatedTensorAlgebraElem ony defined for sequence_type = :iis"))
    end

end



function leading_coefficient_and_zero_TA(p)
    if is_zero(p)
        return zero(typeof(p))
    else
        return leading_coefficient(p)
    end
end



function concatenate_tensors_TA(t1::AbstractArray, t2::AbstractArray)
    # Concatenate tensors by outer product
    reshaped_tensor1 = reshape(t1, size(t1)..., ones(Int, ndims(t2))...)
    reshaped_tensor2 = reshape(t2, ones(Int, ndims(t1))..., size(t2)...)
    return reshaped_tensor1 .* reshaped_tensor2
end


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
    res = []
    comb = zeros(Int, k)

    function backtrack(pos, start)
        if pos > k
            push!(res, copy(comb))
            return
        end
        for i in start:n
            comb[pos] = arr[i]
            backtrack(pos+1, i)
        end
    end

    return res
end



function sigAxis_TA_ClosedForm(T::TruncatedTensorAlgebra)
    if T.sequence_type != :iis
        error("sigAxis_TA only defined for sequence_type = :iis")
    end

    R = base_algebra(T)
    d = base_dimension(T)
    k = truncation_level(T)

    E = typeof(one(R))                 # tipe of the element
    elem_out = Vector{Array{E}}(undef, k+1)

    # ======================
    # Level 0 (tensor 0-dimensional)
    # ======================
    elem_out[1] = fill(one(R), ())    # array 0-dimensional

    # ======================
    # Level 1
    # ======================
    elem_out[2] = fill(one(R), d)

    # ======================
    # Levels >= 2
    # ======================
    for j in 2:k
        dims = ntuple(_ -> d, j)
        tensor_j = zeros(E, dims...)

        for idx in combinations_with_replacement(1:d, j)
            counts = Dict{Int,Int}()

            for i in idx
                counts[i] = get(counts, i, 0) + 1
            end

            denom = prod(factorial(c) for c in values(counts))
            value = 1 // denom

            tensor_j[idx...] = try
                R(value)
            catch
                convert(E, value)
            end
        end

        elem_out[j+1] = tensor_j
    end

    return TruncatedTensorAlgebraElem(T, elem_out)
end




#With chen
function sig_axis_TA(T::TruncatedTensorAlgebra{R}) where R
    k = truncation_level(T)
    d = ambient_dimension(T) 

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


function sig_pwln_TA_Congruence(T::TruncatedTensorAlgebra{R}, coeffs::AbstractMatrix) where R
    # 1) Obtain moment path (sig_mono_TA)
    mono_path = sig_axis_TA(T)

    # 2) Apply coeficients with matricial congruence
    poly_path = matrix_tensorAlg_congruence_TA(coeffs, mono_path)

    return poly_path
end

function sig_pwl_TA_chen(T::TruncatedTensorAlgebra{R}, P::AbstractMatrix{E}) where {R,E}
    d = ambient_dimension(T)
    @assert size(P,2) == d "Dimensions mismatch"

    seg_vecs = [P[i+1, :] .- P[i, :] for i in 1:size(P,1)-1]
    seg_sigs = [sig_segment_TA(T, seg_vecs[i]) for i in 1:length(seg_vecs)]
    return prod(seg_sigs)
end



function matrix_tensorAlg_congruence_TA(matrix::AbstractMatrix, x::TruncatedTensorAlgebraElem)
    y = deepcopy(x)

    for lvl in 2:length(y.elem)
        y.elem[lvl] = matrix_tensor_congruence_TA(matrix, y.elem[lvl])
    end

    return y
end


##################
# arithmetic tensor sequences (only for :iis)
##################

function Base.:(==)(a::TruncatedTensorAlgebraElem, b::TruncatedTensorAlgebraElem)
    if parent(a).sequence_type == :iis && parent(b).sequence_type == :iis
        return parent(a) == parent(b) && tensor_sequence(a) == tensor_sequence(b)
    else
        throw(ArgumentError("== only defined for sequence_type == :iis"))
    end
end

function Base.:+(a::TruncatedTensorAlgebraElem, b::TruncatedTensorAlgebraElem)
    if parent(a).sequence_type == :iis && parent(b).sequence_type == :iis
        A = parent(a)
        res_seq = tensor_sequence(a) + tensor_sequence(b)
        return TruncatedTensorAlgebraElem(A, res_seq)
    else
        throw(ArgumentError("+ only defined for sequence_type == :iis"))
    end
end

function Base.:-(a::TruncatedTensorAlgebraElem, b::TruncatedTensorAlgebraElem)
    if parent(a).sequence_type == :iis && parent(b).sequence_type == :iis
        A = parent(a)
        res_seq = tensor_sequence(a) - tensor_sequence(b)
        return TruncatedTensorAlgebraElem(A, res_seq)
    else
        throw(ArgumentError("- only defined for sequence_type == :iis"))
    end
end

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






# -------------------------------
# Principal function (dispatch)
# -------------------------------
function sigAxis_ClosedForm(T::TruncatedTensorAlgebra{R}, m::Int, n::Int) where R
    if m * n != ambient_dimension(T)
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

    elem_out = Vector{Array{eltype(sig_m.elem[2])}}(undef, k)

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
function sigAxis_p2id_ClosedForm(T::TruncatedTensorAlgebra{R}, m::Int, n::Int) where R
    k = truncation_level(T)
    d=ambient_dimension(T)

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

        perms = Combinatorics.permutations(1:j)

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
    if m * n != ambient_dimension(T)
        error("m * n != d")
    end

    if T.sequence_type == :p2id
        return sigAxis_p2id_Chen(T, m, n)
    elseif T.sequence_type == :p2
        return sigAxis_p2id_Chen(T, m, n)
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

    elem_out = Vector{Array{eltype(sig_m.elem[2])}}(undef, k)

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
function sigAxis_p2id_Chen(T::TruncatedTensorAlgebra{R}, m::Int, n::Int) where R
    k = truncation_level(T)
    d=ambient_dimension(T)

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

        perms = Combinatorics.permutations(1:j)

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
    R = TTA.base_ring
    k = TTA.truncation_level
    d = TTA.ambient_dimension

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
    R = TTA.base_ring
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

        perms = collect(Combinatorics.permutations(1:j))
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
    if TTA.seq_type == :p2id
        return moment_membrane_p2id(TTA, m, n)
    elseif TTA.seq_type == :p2
        return moment_membrane_p2(TTA, m, n)
    else
        error("sequence_type must be :p2id or :p2")
    end
end

"""
    mode_product(T::AbstractArray, A::AbstractMatrix, mode::Int, R)

Compute the mode-n product of tensor `T` with matrix `A` along the specified `mode`.

- `T` : input tensor of any order (dimensions can be for p2id or p2)
- `A` : matrix of size (d_new × d) to multiply along the `mode` dimension
- `mode` : the axis along which to multiply
- `R` : base ring for elements (e.g., QQMPolyRing)

Returns a new tensor with dimension `d_new` along the `mode` axis.
- If `T` is a level-1 tensor from TruncatedTensorAlgebra, multiplies a vector of ones.
- For higher levels, performs a linear transformation along the selected mode.
"""
function mode_product(T::AbstractArray, A::AbstractMatrix, mode::Int, R)
    d_new, d = size(A)

    # Check that mode dimension matches
    dims = size(T)
    dims[mode] == d || error("Mode $mode has size $(dims[mode]), expected $d")

    # Convert A to base ring elements
    A_ring = reshape([one_elem(R) * a for a in A], size(A))

    # Move the mode-th axis to the first dimension
    perm = (mode, setdiff(1:ndims(T), mode)...)
    T_perm = permutedims(T, perm)

    # Matricize tensor: mode dimension becomes rows
    T_mat = reshape(T_perm, d, :)

    # Multiply matrix by tensor matricization
    T_out_mat = A_ring * T_mat   # Linear transform along mode

    # Reshape back to tensor with new mode dimension
    new_dims = (d_new, dims[1:mode-1]..., dims[mode+1:end]...)
    T_out = reshape(T_out_mat, new_dims)

    # Permute axes back to original order
    invp = invperm(perm)
    return permutedims(T_out, invp)
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
    d == X.ambient_dimension || error("size(A,2) must match ambient dimension")

    R = X.base_ring
    k = X.truncation_level
    seq_type = X.seq_type

    # Create new TTA with updated ambient dimension
    X_new = TruncatedTensorAlgebra(R, d_new, k, seq_type=seq_type)

    for j in 1:k
        T = X.tensor_sequence[j]

        # Level 0 (scalar one) remains the same
        if T === one_elem(R)
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
            perms = collect(Combinatorics.permutations(1:j))
            T_perm = Array{typeof(one_elem(R))}(undef, (ntuple(_ -> d_new, j)..., factorial(j)))

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
            error("seq_type must be :p2 or :p2id")
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
    d == T.ambient_dimension || error("size(a,2) must match T.ambient_dimension")

    R = T.base_ring
    k = T.truncation_level
    seq_type = T.seq_type

    # Step 1: Build moment membrane tensor in same algebra as T
    T_moment = momentMembraneTensor(T)   # <-- returns TruncatedTensorAlgebra with same algebra

    # Step 2: Apply matrix-tensor congruence using a
    T_new = applyMatrixToTTA(a, T_moment)

    return T_new
end



