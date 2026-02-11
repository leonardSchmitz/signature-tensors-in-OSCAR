export  
  concatenate_tensors, 
  matrix_tensor_multiply, 
  matrix_tensor_congruence, 
  matrix_tensor_congruence_TA,
#  concatenate_tensors_TA,
  matrix_tensor_congruence_TA


#TODO: remove ≤ 

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
    @assert 1 <= index <= ndims(tensor) "Index must be a valid axis of the tensor."
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

    # 1) Matric and tensor type 
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


