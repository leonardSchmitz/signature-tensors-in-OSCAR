export  
  #_index_number, # internal, to be removed
  #_index2number, 
  #_polyRingSigTransfElem2dimOrder, 
  #polynomial_ring_sig_transform,  # tensor operations
  #Cmono_seq, 
  #Caxis_seq, 
  concatenate_tensors, 
  matrix_tensor_multiply, 
  matrix_tensor_congruence
  #matrix_tensorSeq_congruence, 
  #sig_segment_standard_direction, 
  #sig_segment, 
  #univariate_free_signature_algebra,  # univariate free signature algebra
  #free_sig_elem, 
  #trunc_univariate, 
  #graded_component_univariate, 
  #exp_trunc_univariate, 
  #log_trunc_univariate, 
  #inverse_trunc_univariate, 
  #multivariate_free_signature_algebra, # multivariate free signature algebra
  #free_sig_elem_from_sample, 
  #trunc_multivariate, 
  #graded_component_multivariate, 
  #exp_trunc_multivariate, 
  #log_trunc_multivariate, 
  #inverse_trunc_multivariate,
  #prod_trunc_multivariate, 
  #power_trunc_multivariate,
  #eval_to_tensorSeq_univariate, # evaluation from free algebra to tensor sequence (univariate) 
  #eval_to_tensorSeq_multivariate, # evaluation from free algebra to tensor sequence (multivariate)
  #free_barycenter, # free barycenter map   
  #poly_for_log_bary,
  #free_barycenter_2samples, 
  #poly_for_log_bary_inverse
  
  

#---------------------
# nice printing
#---------------------

#function _index_number(i::Int)
#  dgs = reverse(digits(i))
#  return join(["₀₁₂₃₄₅₆₇₈₉"[3*d+1] for d in dgs], "")
#end
#
#function _index2number(i::String)
#    subscript_map = Dict('₀' => 0, '₁' => 1, '₂' => 2, '₃' => 3, '₄' => 4, 
#                         '₅' => 5, '₆' => 6, '₇' => 7, '₈' => 8, '₉' => 9)
#    digits = [subscript_map[c] for c in i]
#    result = 0
#    for d in digits
#        result = result * 10 + d
#    end
#    return result
#end
#
#function _superscript_number(i::Int)
#    dgs = reverse(digits(i))
#    superscript_digits = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']
#    return '⁽'*join([superscript_digits[d + 1] for d in dgs])*'⁾'
#end

## for function polynomial_ring_sig_transform(_dim::Int,_order::Int)
#function _magic_symbols_mat_a(_dim::Int, _order::Int)
#  a = Array{String}(undef, _dim, _order)
#  for i in 1:_dim
#    for j in 1:_order
#      if i > 9 || j > 9
#        a[i, j] = "a$(_index_number(i)),$(_index_number(j))"
#      else
#        a[i, j] = "a$(_index_number(i))$(_index_number(j))"
#      end
#    end
#  end
#  return a
#end
#
## for function univariate_free_signature_algebra(_trunc_level::Int)
#function _magic_symbols_sig(trunc_level::Int=4)
#  s = Vector{String}(undef, trunc_level)
#  for i in 1:trunc_level
#    s[i] = "s$(_superscript_number(i))"
#  end
#  return s
#end

#### very old, should be omitted for ever:
#### for function multivariate_free_signature_algebra(_trunc_level::Int,_num_samples::Int)
###function _magic_symbols_sig_multivariate(trunc_level::Int, num_samples::Int)
###  s = Array{String}(undef, trunc_level,num_samples)
###  for i in 1:trunc_level
###     for j in 1:num_samples
###      if i > 9 || j > 9
###        #s[i, j] = "s$(_index_number(i))₋$(_index_number(j))"
###        s[i, j] = "s$(_index_number(j))₋$(_superscript_number(i))"
###      else
###        #s[i, j] = "s$(_index_number(i))$(_index_number(j))"
###        s[i, j] = "s$(_index_number(j))$(_superscript_number(i))"
###      end
###    end
###  end
###  return s
###end

# for function multivariate_free_signature_algebra(_trunc_level::Int,_num_samples::Int)
#function _magic_symbols_sig_multivariate(trunc_level::Int, num_samples::Int, invert_bary::Bool = false)
#  s = Array{String}(undef, trunc_level,num_samples)
#  for i in 1:trunc_level
#     for j in 1:num_samples
#      if i > 9 || j > 9
#        #s[i, j] = "s$(_index_number(i))₋$(_index_number(j))"
#        s[i, j] = "s$(_index_number(j))₋$(_superscript_number(i))"
#      else
#        #s[i, j] = "s$(_index_number(i))$(_index_number(j))"
#        s[i, j] = "s$(_index_number(j))$(_superscript_number(i))"
#      end
#    end
#  end
#  if invert_bary 
#     for i in 1:trunc_level
#        s[i, num_samples] = "y$(_superscript_number(i))"
#     end 
#  end 
#  return s
#end



#---------------------
# tensor operations 
#---------------------

#function polynomial_ring_sig_transform(_dim::Int,_order::Int)
#   return  polynomial_ring(QQ, _magic_symbols_mat_a(_dim,_order))
##   return  polynomial_ring(QQ, :a => (1:_dim,1:_order))
#end

#function _polyRingSigTransfElem2dimOrder(_f::QQMPolyRingElem)
#   A = parent(_f)
#   num_gens = number_of_generators(A)
#   str = string(gens(A)[num_gens])[2:end]
#   if length(str) == 2
#     return (_index2number(str[1:1]),_index2number(str[end:end]))
#   else 
#     return 0
#   end 
#end

#function _polyRingSigTransfElem2dimOrder(_f::QQMPolyRingElem)
#   A = parent(_f)
#   num_gens = number_of_generators(A)
#   str = string(gens(A)[num_gens])[2:end]
#   return  parse.(Int, match(r"\[(\d+),\s*(\d+)\]", str).captures)
#end

#function _Cmono(_k::Int, _order::Int, _R::QQMPolyRing)
#   #res = Array{QQMPolyRingElem}(undef,m)
#   res = Array{QQMPolyRingElem}(undef, _order*ones(Int,_k)...)
#   for idx in CartesianIndices(res)
#       res[idx] = QQ(prod(Tuple(idx)),prod(cumsum(Tuple(idx))))*one(_R)
#   end
#   return res
#end
#
#function Cmono_seq(_trunc_level::Int, _order::Int, _R::QQMPolyRing)
#   Cmono = [_Cmono(_k,_order, _R) for _k in (1:_trunc_level)]
#   return Cmono
#end

## TODO  (use Chen and linear path here)
#function _Caxis(_k::Int, _order::Int, _R::QQMPolyRing)
#function Caxis_seq(_k::Int, _order::Int, _R::QQMPolyRing)


function concatenate_tensors(tensor1::Array, tensor2::Array)
    reshaped_tensor1 = reshape(tensor1, size(tensor1)..., ones(Int, ndims(tensor2))...)
    reshaped_tensor2 = reshape(tensor2, ones(Int, ndims(tensor1))..., size(tensor2)...)
    return reshaped_tensor1 .* reshaped_tensor2
end

function matrix_tensor_multiply(matrix::Array, tensor::Array, index::Int)
    """
    Perform matrix-tensor multiplication along the specified axis of a k-tensor.
    Parameters:
        matrix::Array: The matrix to multiply (2D array).
        tensor::Array: The input tensor (k-dimensional array).
        index::Int: The axis of the tensor to perform the multiplication.
    Returns:
        A new tensor resulting from the multiplication.
    """
    # Ensure the dimensions match for contraction
    @assert 1 ≤ index ≤ ndims(tensor) "Index must be a valid axis of the tensor."
    @assert size(tensor, index) == size(matrix, 2) "Matrix and tensor dimensions do not align!"
    # Move the target axis to the first dimension
    perm = [index; 1:index-1; index+1:ndims(tensor)]
    tensor_perm = permutedims(tensor, perm)
    # Reshape the tensor to (size along index, rest)
    reshaped_tensor = reshape(tensor_perm, size(tensor, index), :)
    # Perform matrix multiplication
    result = matrix * reshaped_tensor  # Result shape: (matrix rows, remaining dimensions)
    # Reshape back to the original dimensions
    result_shape = (size(matrix, 1), size(tensor)[[1:index-1; index+1:end]]...)
    result = reshape(result, result_shape)
    # Undo the permutation to restore original axes order
    inv_perm = invperm(perm)
    return permutedims(result, inv_perm)
end

function matrix_tensor_congruence(matrix::Array, tensor::Array)
    k = length(size(tensor))
    if k == 0 
      return tensor
    end 
    res = tensor
    #res = Array{eltype(tensor)}(undef, size(matrix, 1)*ones(Int,k)...)
    for i in (1:k)
        @assert size(tensor, i) == size(matrix, 2) "Matrix and tensor dimensions do not align!"
        res = matrix_tensor_multiply(matrix, res, i)
    end
    return res
end

## this function can be removed: (its implemented in TuncTensorSeq)
#function matrix_tensorSeq_congruence(matrix::Array, tensorSeq::Array)
#    return [matrix_tensor_congruence(matrix, tensor) for tensor in tensorSeq]
#end


##---------------------
## univariate free signature algebra
##---------------------
#
#function univariate_free_signature_algebra(_trunc_level::Int)
#  return free_associative_algebra(QQ, _magic_symbols_sig(_trunc_level))
#end
#
#function free_sig_elem(_s::Vector{Generic.FreeAssociativeAlgebraElem{QQFieldElem}})
#  trunc_level = size(_s)[1]
#  return sum(_s[i] for i in (1:trunc_level)) + one(_s[1])
#end
#
#function trunc_univariate(f::FreeAssociativeAlgebraElem,_trunc_level::Int)
#  A = parent(f)
#  res = zero(A)
#  if iszero(f)
#    return res
#  end
#  for i in (1:length(f))
#    if sum((f.exps)[i]) <= _trunc_level
#      res += collect(terms(f))[i] ## this does not work ? in my provious testing it worked ... 
#    end
#  end
#  return res
#end
#
#function graded_component_univariate(f::FreeAssociativeAlgebraElem,ell::Int)
#  A = parent(f)
#  res = zero(A)
#  if iszero(f)
#    return res
#  end
#  for i in (1:length(f))
#    if sum((f.exps)[i]) == ell
#      res += collect(terms(f))[i]
#    end
#  end
#  return res
#end
#
#function exp_trunc_univariate(f::FreeAssociativeAlgebraElem,_trunc_level::Int)
#  A = parent(f)
#  res = trunc_univariate(f,_trunc_level)
#  for ell in (2:_trunc_level)
#    res += trunc_univariate(QQ(1,factorial(ell))*f^ell,_trunc_level)
#  end
#  return trunc_univariate(res+one(A),_trunc_level)
#end
#
#
#function log_trunc_univariate(f::FreeAssociativeAlgebraElem,_trunc_level::Int)
#  A = parent(f)
#  res = zero(A)
#  for ell in (1:_trunc_level)
#    res += trunc_univariate(QQ((-1)^(ell+1),ell)*(f-one(A))^ell,_trunc_level)
#  end
#  return trunc_univariate(res,_trunc_level)
#end
#
#function inverse_trunc_univariate(_f::FreeAssociativeAlgebraElem,_trunc_level::Int)
#  log_f = log_trunc_univariate(_f,_trunc_level)
#  return exp_trunc_univariate(-1*log_f,_trunc_level)
#end
#
##---------------------
## multivariate free signature algebra
##---------------------
#
#
#function multivariate_free_signature_algebra(_trunc_level::Int,_num_samples::Int, invert_bary::Bool = false)
#   return free_associative_algebra(QQ, _magic_symbols_sig_multivariate(_trunc_level,_num_samples, invert_bary))
#end 
#
#function free_sig_elem_from_sample(_sample_index::Int, _s::Matrix{Generic.FreeAssociativeAlgebraElem{QQFieldElem}})
#  trunc_level = size(_s)[1]
#  return sum(_s[i,_sample_index] for i in (1:trunc_level)) + one(_s[1,1])
#end
#
### TODO can we get the truncation level from f (from the shape of s)
## symbols(parent(sig_1)) works, why not vars? 
#function trunc_multivariate(f::FreeAssociativeAlgebraElem,_trunc_level::Int)
#  A = parent(f)
#  res = zero(A)
#  if iszero(f)
#    return res
#  end
#  for i in (1:length(f))
#    exp_vec = (f.exps)[i].-1
#    mod_exp_vec = exp_vec.%_trunc_level
#    if sum(mod_exp_vec.+1) <= _trunc_level
#      res += collect(terms(f))[i]
#    end
#  end
#  return res
#end
#
#function graded_component_multivariate(f::FreeAssociativeAlgebraElem,ell::Int,_trunc_level::Int)
#  A = parent(f)
#  res = zero(A)
#  if iszero(f)
#    return res
#  end
#  for i in (1:length(f))
#    exp_vec = (f.exps)[i].-1
#    mod_exp_vec = exp_vec.%_trunc_level
#    if sum(mod_exp_vec.+1) == ell
#      res += collect(terms(f))[i]
#    end
#  end
#  return res
#end
#
#function prod_trunc_multivariate(_f::FreeAssociativeAlgebraElem,_g::FreeAssociativeAlgebraElem,_trunc_level::Int)
#  return trunc_multivariate(_f*_g,_trunc_level::Int)  # could be optimized, truncate on each sumand 
#end 
#
#function power_trunc_multivariate(_f::FreeAssociativeAlgebraElem,_a::Int,_trunc_level::Int)
#  # could also be optimized using fast exponentiation
#  res = _f
#  for i in (2:_a)
#    res = prod_trunc_multivariate(res,_f,_trunc_level)
#  end 
#  return res 
#end
#
### TODO can we get the truncation level from f (from the shape of s)
## symbols(parent(sig_1)) works, why not vars? 
#function exp_trunc_multivariate(_f::FreeAssociativeAlgebraElem,_trunc_level::Int)
#  A = parent(_f)
#  f = trunc_multivariate(_f,_trunc_level)
#  res = f
#  for ell in (2:_trunc_level)
#    res += QQ(1,factorial(ell))*power_trunc_multivariate(f,ell,_trunc_level)
#  end
#  return res+one(A)
#end
#
#
### TODO can we get the truncation level from f (from the shape of s)
## symbols(parent(sig_1)) works, why not vars? 
#function log_trunc_multivariate(_f::FreeAssociativeAlgebraElem,_trunc_level::Int)
#  A = parent(_f)
#  res = zero(A)
#  f = trunc_multivariate(_f,_trunc_level)-one(A)
#  for ell in (1:_trunc_level)
#    res += QQ((-1)^(ell+1),ell)*power_trunc_multivariate(f,ell,_trunc_level)
#  end
#  return res
#end
#
#function inverse_trunc_multivariate(_f::FreeAssociativeAlgebraElem,_trunc_level::Int)
#  log_f = log_trunc_multivariate(_f,_trunc_level)
#  return exp_trunc_multivariate(-1*log_f,_trunc_level)
#end
#
#
##---------------------
## evaluation from free algebra to tensor sequence (univariate) 
##---------------------
#
#function eval_to_tensorSeq_univariate(poly::FreeAssociativeAlgebraElem, tensorSeq::Array)
#   return [eval_to_tensor_univariate(poly,tensorSeq,i) for i in (1:length(tensorSeq))]
#end
#
#function eval_to_tensor_univariate(poly::FreeAssociativeAlgebraElem, tensorSeq::Array, index::Int)
#   fi = graded_component_univariate(poly,index)
#   res = 0*tensorSeq[index]
#   for t in (1:length(fi))
#       fit = collect(terms(fi))[t]
#       temp = leading_coefficient(fit).*tensorSeq[fit.exps[1][1]] #exps provides vector
#       for fak in (2:length(fit.exps[1]))
#          temp = concatenate_tensors(temp,tensorSeq[fit.exps[1][fak]]) # exps provides vector
#       end
#       res += temp
#   end
#   return res
#end
#
##---------------------
## evaluation from free algebra to tensor sequence (multivariate) 
##---------------------
#
#function eval_to_tensorSeq_multivariate(poly::FreeAssociativeAlgebraElem, seq_tensor_seq::Array)
#   trunc_level = length(seq_tensor_seq[1])
#   return [eval_to_tensor_multivariate(poly,seq_tensor_seq,i) for i in (1:trunc_level)]
#end
#
#function eval_to_tensor_multivariate(poly::FreeAssociativeAlgebraElem, seq_tensor_seq::Array, index::Int)
#   num_samples = length(seq_tensor_seq)
#   trunc_level = length(seq_tensor_seq[1])
#   fi = graded_component_multivariate(poly,index,trunc_level)
#   res = 0*seq_tensor_seq[1][index] # zero tensor of correct size
#   for t in (1:length(fi))
#       fit = collect(terms(fi))[t]
#       first_exp = fit.exps[1][1]
#       sample_var = div(first_exp-1,trunc_level)+1
#       tensor_var = mod(first_exp-1,trunc_level)+1
#       temp = leading_coefficient(fit).*seq_tensor_seq[sample_var][tensor_var] 
#       for fak in (2:length(fit.exps[1]))
#          temp_exp = fit.exps[1][fak]
#          sample_var = div(temp_exp-1,trunc_level)+1
#          tensor_var = mod(temp_exp-1,trunc_level)+1
#          temp = concatenate_tensors(temp,seq_tensor_seq[sample_var][tensor_var])
#       end
#       res += temp
#   end
#   return res
#end
#
#
###---------------------
### signature of axis path 
###---------------------
##
##function _one_hot(_i::Int,_n::Int,_R::QQMPolyRing)
##  res = zeros(_R,_n)
##  res[_i] = one(_R)
##  return res 
##end 
##
##function sig_segment_standard_direction(_i::Int,_order::Int,_trunc_level::Int,_R::QQMPolyRing)
##   return sig_segment(_one_hot(_i,_order,_R),_trunc_level)
##end
##
##function sig_segment(v::Vector{QQMPolyRingElem},_trunc_level::Int)
##   n = length(v)
##   C = Cmono_seq(_trunc_level,1,parent(v[1]))
##   return matrix_tensorSeq_congruence(v,C)
##end 
##
##
##function Caxis_seq(_trunc_level::Int, _order::Int, _R::QQMPolyRing)
##   sample = [sig_segment_standard_direction(i,_order,_trunc_level,_R) for i in (1:_order)]
##   F,s = multivariate_free_signature_algebra(_trunc_level,_order)
##   chen = trunc_multivariate(prod([free_sig_elem_from_sample(i,s) for i in (1:_order)]),_trunc_level)
##   return eval_to_tensorSeq_multivariate(chen,sample)
##end 
#
##---------------------
## barycenter (truncated)
##---------------------
#
#
### TODO: for N=2 more efficient formula which works for all truncation levels
### incorperate other weights w. 
##
##function _poly_for_grad_log_bary(c::Array,_grad_comp::Int)
##  num_samples = size(c)[2]
##  if _grad_comp == 1
##    return QQ(1,num_samples)*sum(c[1,i] for i in (1:num_samples))
##  elseif _grad_comp == 2
##    return QQ(1,num_samples)*sum(c[2,i] for i in (1:num_samples))
##  elseif _grad_comp == 3
##    temp = (c[3,:]
##            -QQ(1,12).*c[1,:].*c[1,:].*_poly_for_grad_log_bary(c,1)
##            +QQ(2,12).*c[1,:].*_poly_for_grad_log_bary(c,1).*c[1,:]
##            -QQ(1,12).*_poly_for_grad_log_bary(c,1).*c[1,:].*c[1,:])
##    return QQ(1,num_samples)*sum(temp)
##  else
##    error("_poly_for_grad: Graded component >= 4 not supported yet.") 
##  end
##end
##
##function poly_for_log_bary(_trunc_level::Int,_num_samples::Int)
##   F,s = multivariate_free_signature_algebra(_trunc_level,_num_samples)
##   log_sample = Array{FreeAssociativeAlgebraElem}(undef, _trunc_level,_num_samples);
##   for i in (1:_trunc_level)
##     for j in (1:_num_samples)
##       log_sample[i,j] = graded_component_multivariate(free_sig_elem_from_sample(j,s),i,_trunc_level)
##     end
##   end
##   res = _poly_for_grad_log_bary(log_sample,1)
##   for i in (2:_trunc_level)
##     res += _poly_for_grad_log_bary(log_sample,i)
##   end
##   return trunc_multivariate(res,_trunc_level)
##end
##
##function free_barycenter(_trunc_level::Int,_num_samples::Int)
##   F,s = multivariate_free_signature_algebra(_trunc_level,_num_samples)
##   log_sample = Array{FreeAssociativeAlgebraElem}(undef, _trunc_level,_num_samples);
##   for i in (1:_trunc_level)
##     for j in (1:_num_samples)
##       log_sample[i,j] = graded_component_multivariate(log_trunc_multivariate(free_sig_elem_from_sample(j,s),_trunc_level),i,_trunc_level)
##     end
##   end
##   res = _poly_for_grad_log_bary(log_sample,1)
##   for i in (2:_trunc_level)
##     res += _poly_for_grad_log_bary(log_sample,i)
##   end 
##   return exp_trunc_multivariate(res,_trunc_level)
##end 
##
##
##function free_barycenter_2samples(_trunc_level::Int)
##   F,s = multivariate_free_signature_algebra(_trunc_level,2)
##   x = free_sig_elem_from_sample(1,s)
##   x_inv = inverse_trunc_multivariate(x,_trunc_level)
##   y = free_sig_elem_from_sample(2,s)
##   chen = trunc_multivariate(x_inv*y,_trunc_level)
##   temp = exp_trunc_multivariate(1//2*log_trunc_multivariate(chen,_trunc_level),_trunc_level)
##   return trunc_multivariate(x*temp,_trunc_level)
##end
##
##
##function _poly_for_grad_log_bary_inverse(c::Array,_grad_comp::Int)
##  num_samples = size(c)[2]
##  if _grad_comp == 1
##    return num_samples*(c[1,num_samples] - QQ(1,num_samples)*sum(c[1,i] for i in (1:num_samples-1)))
##  elseif _grad_comp == 2
##    return num_samples*(c[2,num_samples] - QQ(1,num_samples)*sum(c[2,i] for i in (1:num_samples-1)))
##  elseif _grad_comp == 3
##    temp = (c[3,:]  ## todo 
##            -QQ(1,12).*c[1,:].*c[1,:].*_poly_for_grad_log_bary(c,1)
##            +QQ(2,12).*c[1,:].*_poly_for_grad_log_bary(c,1).*c[1,:]
##            -QQ(1,12).*_poly_for_grad_log_bary(c,1).*c[1,:].*c[1,:])
##    return QQ(1,num_samples)*sum(temp)
##  else
##    error("_poly_for_grad: Graded component >= 4 not supported yet.") 
##  end
##end
##
##function poly_for_log_bary_inverse(_trunc_level::Int,_num_samples::Int)
##   F,s = multivariate_free_signature_algebra(_trunc_level,_num_samples, true)
##   log_sample = Array{FreeAssociativeAlgebraElem}(undef, _trunc_level,_num_samples);
##   for i in (1:_trunc_level)
##     for j in (1:_num_samples)
##       log_sample[i,j] = graded_component_multivariate(free_sig_elem_from_sample(j,s),i,_trunc_level)
##     end
##   end
##   res = _poly_for_grad_log_bary_inverse(log_sample,1)
##   for i in (2:_trunc_level)
##     res += _poly_for_grad_log_bary_inverse(log_sample,i)
##   end
##   return trunc_multivariate(res,_trunc_level)
##end
##
##
##function poly_free_sig(_trunc_level::Int,_num_sampes::Int)
##   F,s = multivariate_free_signature_algebra(_trunc_level,_num_samples+1, true)
##   samples = [free_sig_elem_from_sample(i,s) for i in (1:_num_samples)]
##   bary = free_sig_elem_from_sample(_num_samples+1,s)
##   bk = sum([prod_trunc_multivariate(bary,log_trunc_multivariate(xi,k),k) for xi in samples])
##   return 0
##end
#



