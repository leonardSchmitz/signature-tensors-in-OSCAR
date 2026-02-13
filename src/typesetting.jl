export 
  _index_number, 
  _superscript_number, 
  polynomial_ring_sig_transform, 
  matrix_to_latex,
  tensor_to_latex_folded 


function _index_number(i::Int)
  dgs = reverse(digits(i))
  return join(["₀₁₂₃₄₅₆₇₈₉"[3*d+1] for d in dgs], "")
end

function _index2number(i::String)
    subscript_map = Dict('₀' => 0, '₁' => 1, '₂' => 2, '₃' => 3, '₄' => 4,
                         '₅' => 5, '₆' => 6, '₇' => 7, '₈' => 8, '₉' => 9)
    digits = [subscript_map[c] for c in i]
    result = 0
    for d in digits
        result = result * 10 + d
    end
    return result
end

function _superscript_number(i::Int)
    dgs = reverse(digits(i))
    superscript_digits = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']
    return '⁽'*join([superscript_digits[d + 1] for d in dgs])*'⁾'
end

"""
    polynomial_ring_sig_transform(d::Int, m::Int)

Constructs a polynomial ring `QQ[a_ij]` with `d` rows and `m` columns and returns 
a tuple `(R, a)` where `R` is the ring and `a` is a matrix of generators.

# Example
    R, a = polynomial_ring_sig_transform(2, 3)
"""
polynomial_ring_sig_transform

function polynomial_ring_sig_transform(_dim::Int,_order::Int)
   return  polynomial_ring(QQ, _magic_symbols_mat_a(_dim,_order))
#   return  polynomial_ring(QQ, :a => (1:_dim,1:_order))
end

function _magic_symbols_mat_a(_dim::Int, _order::Int)
  a = Array{String}(undef, _dim, _order)
  for i in 1:_dim
    for j in 1:_order
      if i > 9 || j > 9
        a[i, j] = "a$(_index_number(i)),$(_index_number(j))"
      else
        a[i, j] = "a$(_index_number(i))$(_index_number(j))"
      end
    end
  end
  return a
end


function matrix_to_latex(M)
    rows = [join(row, " & ") for row in eachrow(M)]
    body = join(rows, " \\")
    return "\begin{bmatrix}" * body * "\end{bmatrix}"
end

function tensor_to_latex_folded(G::Array{T,3}; name="G") where T
    d1, d2, d3 = size(G)
    @assert d1 == d2 == d3 "Tensor must be cubic d×d×d"
    rows_latex = String[]
    for i in 1:d1
        # Collect all entries in the i-th “row of slices”
        row_entries = []
        for k in 1:d3
            for j in 1:d2
                push!(row_entries, string(G[i,j,k]))
            end
            if k < d3
                push!(row_entries, "\vrule")  # vertical separator between slices
            end
        end
        # Join entries of this row with &
        push!(rows_latex, join(row_entries, " & "))
    end
    body = join(rows_latex, " \\")
    col_format = join(["c" for _ in 1:(d2*d3 + (d3-1))], "")
    return "\begin{array}{" * col_format * "}" *
           body *
           "\end{array}"
end