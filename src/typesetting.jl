export 
  _index_number, 
  _superscript_number, 
  polynomial_ring_sig_transform


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
