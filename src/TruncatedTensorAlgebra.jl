export TruncatedTensorAlgebra,
       TruncatedTensorAlgebraElem,
       truncation_level,
       ambient_dimension,
       base_dimension,
       base_ring,
       base_algebra,
       tensor_sequence,
       TruncatedTensorAlgebra,
       zero,
       one,
       sig_mono_TA,
       sigAxis_TA_ClosedForm,
       sig_axis_TA,
       sig_poly_TA,
       sig_pwln_TA_Congruence,
       sig_pwl_TA_chen,
       sig_segment_TA,
       sig_segment_standard_direction_TA,
       concatenate_tensors_TA,
       matrix_tensorAlg_congruence_TA,
       matrix_tensor_congruence_TA,
       matrix_tensor_multiply_TA


    



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
    # Nivel 0 (común a todos)
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
    R = base_algebra(T)   # o base_ring(T)
    seq = T.sequence_type

    elem = Vector{Any}(undef, k + 1)

    ####################
    # nivel 0 (siempre)
    ####################
    elem[1] = fill(zero(R), ())   # 0-dimensional scalar

    ####################
    # niveles ≥ 1
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

    # evaluamos en cada nivel 0..k
    res_seq = [evaluate_to_tensor_TA(p, bs, i) for i in 0:k]

    return TruncatedTensorAlgebraElem(A, res_seq)
end



# 1) devuelve 1 del tipo correcto
function one_elem(R)
    if R === QQField
        return QQ(1)          # elemento del campo QQ
    elseif R <: Type
        return one(R)         # si R es tipo de elemento
    else
        return one(R)         # si R es objeto tipo anillo
    end
end

# 2) devuelve 0 del tipo correcto
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

    # seq debe ser un Vector{Array{E}} para algún E
    seq = _Cmono_seqTA(k, d, R0)

    # determinar E a partir de seq
    E = eltype(seq[1])   # tipo de los elementos del array interno
    TruncatedTensorAlgebraElem{R, E}(T, seq)


    return TruncatedTensorAlgebraElem{R, E}(T, seq)
end


function Base.getindex(x::TruncatedTensorAlgebraElem, w...)
    k=length(w)
    if k < 0
        throw(BoundsError(x, k))
    end
    if k > length(x.elem)
        throw(BoundsError(x, k))
    end
    temp= x.elem[k+1 ]
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
# helper: obtener el anillo padre de un elemento (tipo Oscar)
# ------------------------------------------------------------
ring_of(::Type{T}) where T = parent(one(T))
ring_of(x) = ring_of(typeof(x))

# ------------------------------------------------------------
# common_ring: cumple exactamente tu regla
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

    # Anillos
    R_tensor = base_ring(T)

    # Aquí está la corrección:
    # el ring real del que provienen los elementos de la matriz
    R_matrix = parent(matrix[1,1])

    # Regla simple
    Rnew = common_ring(R_tensor, R_matrix)

    # Nuevo álgebra
    Tnew = TruncatedTensorAlgebra(Rnew, d, k, :iis)

    # Convertir tensores al nuevo anillo
    tensorSeq = tensor_sequence(b)
    tensorSeq = [map(x -> Rnew(x), t) for t in tensorSeq]

    # Aplicar congruencia
    resSeq = [
        matrix_tensor_congruence_TA(map(x -> Rnew(x), matrix), t)
        for t in tensorSeq
    ]

    return TruncatedTensorAlgebraElem(Tnew, resSeq)
end


# ----------------------------
# congruencia para vectores
# ----------------------------
function matrix_tensorAlg_congruence_TA(
    v::AbstractVector,
    b::TruncatedTensorAlgebraElem{R,E}
) where {R,E}
    M = reshape(v, 1, length(v))
    return matrix_tensorAlg_congruence_TA(M, b)
end

function sig_segment_TA(T::TruncatedTensorAlgebra{R}, v::Vector{E}) where {R,E}
    # --- validaciones ---
    k = truncation_level(T)
    d = ambient_dimension(T)

    @assert length(v) == d "dimensions do not match"

    # --- aseguramos que el tipo del elemento concuerde con el anillo ---
    @assert typeof(v[1]) <: E "element types do not match"

    # --- creamos un algebra nuevo de dimensión 1 ---
    T1 = TruncatedTensorAlgebra(base_ring(T), 1, k, :iis)

    # --- signature mono en dimensión 1 ---
    C = sig_mono_TA(T1)

    # --- convierte el vector v en matriz columna d×1 ---
    v_mat = reshape(v, d, 1)

    # --- aplica la congruencia matricial (matriz columna * signature 1D) ---
    return matrix_tensorAlg_congruence_TA(v_mat, C)
end


function sig_segment_standard_direction_TA(T::TruncatedTensorAlgebra{R}, _i::Int) where {R}
    d = ambient_dimension(T)
    # vector one-hot en la dirección _i
    v = _one_hot_TA(_i, d, R)
    return sig_segment_TA(T, v)
end



function Base.:*(a::R, b::TruncatedTensorAlgebraElem{R,E}) where {R,E}
    A = parent(b)
    k = truncation_level(A)

    res_seq = tensor_sequence(b)

    # escala el primer nivel
    res_seq[1][] = a * res_seq[1][]

    # escala todos los niveles
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
                throw(ArgumentError("matrix * TruncatedTensorAlgebraElem solo definido para sequence_type = :iis"))
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
    # Concatenate tensors by outer product (como antes)
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

    # caso base: nivel 0
    if index == 0
        res = tensor_sequence(bs[1])[1]   # nivel 0
        res[] = leading_coefficient_and_zero_TA(
                  FreeTruncSigAlgMultivElem_to_polynomial(graded_component(p, 0))
               ) * one(res[])
        return res
    end

    # datos del espacio
    num_vars = length(bs)
    A = parent(bs[1])
    k = truncation_level(A)

    # componente grado index de p
    fi = FreeTruncSigAlgMultivElem_to_polynomial(graded_component(p, index))

    # secuencia de tensores niveles 1..k
    seq_tensor_seq = [tensor_sequence(b)[2:k+1] for b in bs]

    # tensor resultado (tamaño correcto)
    res = zero(seq_tensor_seq[1][index])

    # evaluar término a término
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
    
    # Evaluamos en cada nivel 0..k
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

    backtrack(1, 1)
    return res
end



function sigAxis_TA_ClosedForm(T::TruncatedTensorAlgebra)
    if T.sequence_type != :iis
        error("sigAxis_TA only defined for sequence_type = :iis")
    end

    R = base_algebra(T)
    d = base_dimension(T)
    k = truncation_level(T)

    E = typeof(one(R))                 # tipo del elemento
    elem_out = Vector{Array{E}}(undef, k+1)

    # ======================
    # Nivel 0 (tensor 0-dimensional)
    # ======================
    elem_out[1] = fill(one(R), ())    # array 0-dimensional

    # ======================
    # Nivel 1
    # ======================
    elem_out[2] = fill(one(R), d)

    # ======================
    # Niveles >= 2
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

    # 1) creamos los d segmentos estándar (direcciones base)
    sample = [sig_segment_standard_direction_TA(T, i) for i in 1:d]

    # 2) creamos el free truncated multivariate algebra
    F, s = free_trunc_sig_alg_multiv(k, d)

    # 3) chen product: multiplicación de los generadores 1..d
    chen = prod([free_sig_from_sample(i, F) for i in 1:d])

    # 4) evaluación
    return evaluate_TA(chen, sample)
end


function sig_poly_TA(T::TruncatedTensorAlgebra{R}, coeffs::AbstractMatrix) where R
    # 1) Obtener moment path (sig_mono_TA)
    mono_path = sig_mono_TA(T)

    # 2) Aplicar coeficientes con congruencia matricial
    poly_path = matrix_tensorAlg_congruence_TA(coeffs, mono_path)

    return poly_path
end


function sig_pwln_TA_Congruence(T::TruncatedTensorAlgebra{R}, coeffs::AbstractMatrix) where R
    # 1) Obtener moment path (sig_mono_TA)
    mono_path = sig_axis_TA(T)

    # 2) Aplicar coeficientes con congruencia matricial
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

function matrix_tensor_multiply_TA(matrix::AbstractMatrix, tensor::AbstractArray, index::Int)
    @assert 1 ≤ index ≤ ndims(tensor) "Index must be a valid axis of the tensor."
    @assert size(tensor, index) == size(matrix, 2) "Matrix and tensor dimensions do not align!"

    # 1) Asegurar que matrix esté en el mismo anillo que el tensor
    matrix = convert.(eltype(tensor), matrix)

    perm = [index; 1:index-1; index+1:ndims(tensor)]
    tensor_perm = permutedims(tensor, perm)

    reshaped_tensor = reshape(tensor_perm, size(tensor, index), :)
    result = matrix * reshaped_tensor

    result_shape = (size(matrix, 1), size(tensor)[[1:index-1; index+1:end]]...)
    result = reshape(result, result_shape)

    inv_perm = invperm(perm)
    return permutedims(result, inv_perm)
end

function matrix_tensor_congruence_TA(matrix::AbstractMatrix, tensor::AbstractArray)
    k = ndims(tensor)
    if k == 0
        return tensor
    end

    res = tensor
    for i in 1:k
        res = matrix_tensor_multiply_TA(matrix, res, i)
    end
    return res
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

        res_seq[1][] = a * res_seq[1][]        # nivel 0
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




