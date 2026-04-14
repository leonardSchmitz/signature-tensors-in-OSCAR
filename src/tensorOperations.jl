export  
  concatenate_tensors, 
  matrix_tensor_multiply, 
  matrix_tensor_congruence, 
  matrix_tensor_congruence_TA,
  concatenate_tensors_TA,
  mode_product, 
  axis_core_tensor_normalized_3, 
  add_slice_all_3!, 
  multiply_slice_all_3!, 
  swap_slice_all_3!, 
  matrix_tensor_multiply_3!, 
  matrix_tensor_multiply_3, 
  matrix_tensor_congruence_3!, 
  matrix_tensor_congruence_3, 
  one_hot



function concatenate_tensors(tensor1::Array, tensor2::Array)
    reshaped_tensor1 = reshape(tensor1, size(tensor1)..., ones(Int, ndims(tensor2))...)
    reshaped_tensor2 = reshape(tensor2, ones(Int, ndims(tensor1))..., size(tensor2)...)
    return reshaped_tensor1 .* reshaped_tensor2
end

"""
    Perform matrix-tensor multiplication along the specified axis of a k-tensor.
    Parameters:
        matrix::Array: The matrix to multiply (2D array).
        tensor::Array: The input tensor (k-dimensional array).
        index::Int: The axis of the tensor to perform the multiplication.
    Returns:
        A new tensor resulting from the multiplication.
"""

function matrix_tensor_multiply(matrix::Array, tensor::Array, index::Int)

    @assert 1 <= index <= ndims(tensor) "Index must be a valid axis of the tensor."
    @assert size(tensor, index) == size(matrix, 2) "Matrix and tensor dimensions do not align!"

    perm = [index; 1:index-1; index+1:ndims(tensor)]
    tensor_perm = permutedims(tensor, perm)

    reshaped_tensor = reshape(tensor_perm, size(tensor, index), :)

    result = matrix * reshaped_tensor 

    result_shape = (size(matrix, 1), size(tensor)[[1:index-1; index+1:end]]...)
    result = reshape(result, result_shape)

    inv_perm = invperm(perm)
    return permutedims(result, inv_perm)
end

function matrix_tensor_congruence(matrix::Array, tensor::Array)
    k = length(size(tensor))
    if k == 0 
      return tensor
    end 
    res = tensor

    for i in (1:k)
        @assert size(tensor, i) == size(matrix, 2) "Matrix and tensor dimensions do not align!"
        res = matrix_tensor_multiply(matrix, res, i)
    end
    return res
end


function concatenate_tensors_TA(t1::AbstractArray, t2::AbstractArray)
    # Concatenate tensors by outer product
    reshaped_tensor1 = reshape(t1, size(t1)..., ones(Int, ndims(t2))...)
    reshaped_tensor2 = reshape(t2, ones(Int, ndims(t1))..., size(t2)...)
    return reshaped_tensor1 .* reshaped_tensor2
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


function matrix_tensor_multiply_TA(matrix::AbstractMatrix, tensor::AbstractArray, index::Int)
    @assert 1 <= index <= ndims(tensor) "Index must be a valid axis of the tensor."
    @assert size(tensor, index) == size(matrix, 2) "Matrix and tensor dimensions do not align!"

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


"""
    mode_product(T::AbstractArray, A::AbstractMatrix, mode::Int, R)

Compute the **mode-n product** of a tensor `T` with a matrix `A` along a specified axis (`mode`).

# Arguments
- `T::AbstractArray` : Input tensor of any order (can represent vectors, matrices, or higher-order tensors).  
- `A::AbstractMatrix` : Matrix of size `(d_new × d)` to multiply along the `mode`-th dimension of `T`.  
- `mode::Int` : The axis of `T` along which the multiplication is applied (1-based indexing).  
- `R` : Base ring or type for elements (e.g., `QQMPolyRing`), used to lift matrix entries into the tensor's algebra.

# Returns
A new tensor of the same order as `T` but with the `mode`-th dimension replaced by `d_new`.

# Details
- For **level-1 tensors** (vectors) from a `TruncatedTensorAlgebra`, this multiplies a vector of ones.  
- For **higher-order tensors**, performs a linear transformation along the specified mode by:
  1. Moving the `mode` axis to the first dimension.
  2. Matricizing the tensor along that axis.
  3. Multiplying by `A` (converted to the base ring `R`).
  4. Reshaping the result back to the tensor form.
  5. Permuting the axes back to the original order.
# Example
```julia
using Oscar
T = QQ.(rand(-40:40,3,4,2))          # 3×4×2 tensor
A = QQ.(rand(-20:20,5,3))            # 5×3 matrix to multiply along mode 1
R = QQ               # base ring
T_new = mode_product(T, A, 1, R)
size(T_new)               # returns (5,4,2)
```
Errors
Throws an error if the mode-th dimension of T does not match the number of columns d of A.
"""
function mode_product(T::AbstractArray, A::AbstractMatrix, mode::Int, R)
    d_new, d = size(A)

    dims = size(T)
    dims[mode] == d || error("Mode $mode has size $(dims[mode]), expected $d")

    A_ring = reshape([one(R) * a for a in A], size(A))

    perm = (mode, setdiff(1:ndims(T), mode)...)
    T_perm = permutedims(T, perm)

    T_mat = reshape(T_perm, d, :)

    T_out_mat = A_ring * T_mat   

    new_dims = (d_new, dims[1:mode-1]..., dims[mode+1:end]...)
    T_out = reshape(T_out_mat, new_dims)

    invp = invperm(perm)
    return permutedims(T_out, invp)
end

# functions from Schmitz 2026 - An efficient algorithm for tensor learing, use for test

function axis_core_tensor_normalized_3(_d)
  C = fill(zero(QQ), _d, _d, _d)
  for al in (1:_d)
    for be in (1:_d)
      for ga in (1:_d)
        if al == be && be == ga
          C[al,be,ga] = QQ(1,1)
        end
        if (al < be && be == ga)||(al == be && be < ga)
          C[al,be,ga] = QQ(3,1)
        end
        if al < be && be < ga
          C[al,be,ga] = QQ(6,1)
        end
      end
    end
  end 
  return C
end

function add_slice_all_3!(G,a,i::Int,j::Int)
  G[j,:,:] += a.*G[i,:,:]
  G[:,j,:] += a.*G[:,i,:]
  G[:,:,j] += a.*G[:,:,i]
  return G
end

function multiply_slice_all_3!(G,a,i::Int)
  G[i,:,:] = a.*G[i,:,:]
  G[:,i,:] = a.*G[:,i,:]
  G[:,:,i] = a.*G[:,:,i]
  return G
end

function swap_slice_all_3!(G,i::Int,j::Int)
  temp = G[i,:,:]
  G[i,:,:] = G[j,:,:]
  G[j,:,:] = temp
  temp = G[:,i,:]
  G[:,i,:] = G[:,j,:]
  G[:,j,:] = temp
  temp = G[:,:,i]
  G[:,:,i] = G[:,:,j]
  G[:,:,j] = temp
  return G
end

function matrix_tensor_multiply_3!(matrix::Array, tensor::Array{T,3}, index::Int) where T
    @assert 1 <=index <=3 "Index must be 1, 2, or 3."
    @assert size(tensor, index) == size(matrix, 2) "Dimension mismatch."
    if index == 1

        resh = reshape(tensor, size(tensor,1), :)
        result = matrix * resh
        copyto!(tensor, reshape(result, size(matrix,1), size(tensor,2), size(tensor,3)))
    elseif index == 2
        tensor_perm = permutedims(tensor, (2,1,3))
        resh = reshape(tensor_perm, size(tensor,2), :)
        result = matrix * resh
        result = reshape(result, size(matrix,1), size(tensor,1), size(tensor,3))
        copyto!(tensor, permutedims(result, (2,1,3)))
    else  
        tensor_perm = permutedims(tensor, (3,1,2))
        resh = reshape(tensor_perm, size(tensor,3), :)
        result = matrix * resh
        result = reshape(result, size(matrix,1), size(tensor,1), size(tensor,2))
        copyto!(tensor, permutedims(result, (2,3,1)))
    end
    return tensor
end

function matrix_tensor_congruence_3!(matrix::Array, tensor::Array{T,3}) where T
    @assert size(tensor,1) == size(matrix,2)
    @assert size(tensor,2) == size(matrix,2)
    @assert size(tensor,3) == size(matrix,2)
    matrix_tensor_multiply_3!(matrix, tensor, 1)
    matrix_tensor_multiply_3!(matrix, tensor, 2)
    matrix_tensor_multiply_3!(matrix, tensor, 3)
    return tensor
end

function matrix_tensor_multiply_3(matrix::Array, tensor::Array, index::Int)
    """
    Perform matrix-tensor multiplication along the specified axis of a k-tensor.
    Parameters:
        matrix::Array: The matrix to multiply (2D array).
        tensor::Array: The input tensor (k-dimensional array).
        index::Int: The axis of the tensor to perform the multiplication.
    Returns:
        A new tensor resulting from the multiplication.
    """

    @assert 1 <=index <=ndims(tensor) "Index must be a valid axis of the tensor."
    @assert size(tensor, index) == size(matrix, 2) "Matrix and tensor dimensions do not align!"

    perm = [index; 1:index-1; index+1:ndims(tensor)]
    tensor_perm = permutedims(tensor, perm)

    reshaped_tensor = reshape(tensor_perm, size(tensor, index), :)

    result = matrix * reshaped_tensor 

    result_shape = (size(matrix, 1), size(tensor)[[1:index-1; index+1:end]]...)
    result = reshape(result, result_shape)

    inv_perm = invperm(perm)
    return permutedims(result, inv_perm)
end

function matrix_tensor_congruence_3(matrix::Array, tensor::Array)
    k = length(size(tensor))
    if k == 0 
      return tensor
    end 
    res = tensor

    for i in (1:k)
        @assert size(tensor, i) == size(matrix, 2) "Matrix and tensor dimensions do not align!"
        res = matrix_tensor_multiply_3(matrix, res, i)
    end
    return res
end

function one_hot(_i::Int, _n::Int, _R)
    res = fill(zero(_R), _n)
    res[_i] = one(_R)
    return res
end