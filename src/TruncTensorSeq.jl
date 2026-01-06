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
  ideal_of_entries
  


struct TruncTensorSeq
  base_ring::QQMPolyRing
  trunc_level::Int
  amb_dim::Int
end

truncation_level(F::TruncTensorSeq) = F.trunc_level
ambient_dimension(F::TruncTensorSeq) = F.amb_dim
base_ring(F::TruncTensorSeq) = F.base_ring

struct TruncTensorSeqElem
  parent::TruncTensorSeq
  elem::Vector{Array{QQMPolyRingElem}}
end

Base.parent(a::TruncTensorSeqElem) = a.parent
tensor_sequence(a::TruncTensorSeqElem) = a.elem

function trunc_tensor_seq(base_ring::QQMPolyRing,trunc_level::Int,amb_dim::Int)
  TTS = TruncTensorSeq(base_ring, trunc_level, amb_dim)
  return TTS
end

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

function sig_segment_standard_direction(TTS::TruncTensorSeq,_i::Int)
  d = ambient_dimension(TTS)
  R = base_ring(TTS)
  return sig_segment(TTS,_one_hot(_i,d,R))
end

function sig_axis(TTS::TruncTensorSeq)
  k = truncation_level(TTS) 
  d = ambient_dimension(TTS)
  R = base_ring(TTS)
  sample = [sig_segment_standard_direction(TTS,i) for i in (1:d)]
  F,s = free_trunc_sig_alg_multiv(k,d)
  chen = prod([free_sig_from_sample(i,F) for i in (1:d)])
  return evaluate(chen,sample)
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

function ideal_of_entries(a::TruncTensorSeqElem)
  A = parent(a)
  R = base_ring(A)
  k = truncation_level(A)
  temp = (tensor_sequence(a))[2:k+1]
  gens4learn = vcat([vec(temp[i]) for i in (1:length(temp))]...)
  return  ideal(R,gens4learn)
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

function bary_2samples(a::TruncTensorSeqElem, b::TruncTensorSeqElem)
  return a*exp(QQ(1,2)*log(inv(a)*b))
end

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

function bary_2nd_trunc(bs::Vector{TruncTensorSeqElem})
  N = length(bs)
  return exp(QQ(1,N)*sum(log.(bs))) 
end

function bary_2nd_trunc_closedform(bs::Vector{TruncTensorSeqElem})
  TTSm = parent(bs[1])
  N = length(bs)
  ls = tensor_sequence.(bs)
  vec = QQ(1,N).*sum(ls[i][2] for i in (1:N))
  mat = QQ(1,N).*sum(ls[i][3] for i in (1:N)) - QQ(1,2*N).*sum(ls[i][2]*transpose(ls[i][2]) for i in (1:N)) + QQ(1,2*N^2).*sum(ls[i1][2]*transpose(ls[i2][2]) for i1 in (1:N) for i2 in (1:N))
  return trunc_tensor_seq_elem(TTSm,[ls[1][1],vec,mat])
end



