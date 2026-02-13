export
  tensor_learning_3

#function generic_transform(_d,_max_int=10)
#   A = QQ.(rand(-1*_max_int:_max_int, _d, _d))
#   if is_invertible( matrix(QQ,A))
#     return A
#   else
#     return generic_transform(_d,_max_int)
#   end
#end

function toTriangularTrafo_last_2x2_up(G)
  d = size(G)[1]
  Q = Matrix(identity_matrix(QQ,d))
  if G[d,d-1,d-1] == G[d-1,d,d-1] 
    return Q
  end
  if G[d,d-1,d]!=G[d-1,d,d]
    lambda23 = (G[d,d-1,d-1]- G[d-1,d,d-1])*inv(-G[d, d-1, d]+G[d-1, d, d])
    Q1 = identity_matrix(QQ,d) + lambda23.*QQMatrixE(d,d-1,d)
    add_slice_all_3!(G,lambda23,d,d-1)
    return Matrix(Q1)
  else
    swap_slice_all_3!(G,d-1,d)
    Q0 = Matrix(block_diagonal_matrix([identity_matrix(QQ,d-2),permutation_matrix(QQ,symmetric_group(2)[1])]))
    return Q0
  end
  println("PROBLEM s = d-1")
  return Q0
end

function generic_last_3x3(d)
  return Array(block_diagonal_matrix([identity_matrix(QQ,d-3),matrix(QQ,generic_transform_GL(3,2))]))
end

function toTriangularTrafo_last_3x3_up(G)
  d = size(G)[1]
  if G[d, d-1, d]==G[d-1, d, d]
    Qrand = generic_last_3x3(d)
    matrix_tensor_congruence_3!(Qrand,G)
    return Matrix(toTriangularTrafo_last_3x3_up(G))*Matrix(Qrand)
  end
  lambda23 = (G[d,d-1,d-1]- G[d-1,d,d-1])*inv(-G[d, d-1, d]+G[d-1, d, d])
  Q1 = identity_matrix(QQ,d) + lambda23.*QQMatrixE(d,d-1,d)
  Q = Matrix(Q1)
  add_slice_all_3!(G,lambda23,d,d-1)
  lambda13 = (G[d,d-1,d-2]- G[d-1,d,d-2])*inv(-G[d, d-1, d]+G[d-1, d, d])
  Q2 = identity_matrix(QQ,d) + lambda13.*QQMatrixE(d,d-2,d)
  Q = Matrix(Q2)*Q
  add_slice_all_3!(G,lambda13,d,d-2)
  Q3 = identity_matrix(QQ,d)
  if G[d-1, d-2,d-1]==G[d-2, d-1, d-1]
    Qrand = generic_last_3x3(d)
    matrix_tensor_congruence_3!(Qrand,G)
    return Matrix(toTriangularTrafo_last_3x3_up(G))*Matrix(Qrand)*Q
  end
  lambda12 = (G[d-1,d-2,d-2]- G[d-2,d-1,d-2])*inv(-G[d-1, d-2,d-1]+G[d-2, d-1, d-1])
  Q3 = identity_matrix(QQ,d) + lambda12.*QQMatrixE(d,d-2,d-1)
  add_slice_all_3!(G,lambda12,d-1,d-2)
  Q = Matrix(Q3)*Q
  if G[d-2,d-2,d-2]==zero(QQ)
    Qrand = generic_last_3x3(d)
    matrix_tensor_congruence_3!(Qrand,G)
    return Matrix(toTriangularTrafo_last_3x3_up(G))*Matrix(Qrand)*Q
  else
    return Q
  end
end

function linrel_mat_vec(G,s)
  d = size(G)[1]
  nreq = binomial(d-1,2)
  nrun = d - s
  M = zeros(QQ,nreq,nrun)
  B = zeros(QQ,nreq)
  i = 1
  for _a in (s+1:d)
    for _b in (_a+1:d)
      M[i,:] = [(G[_a,_b,s+gamma]-G[_b,_a,s+gamma]) for gamma in (1:nrun)]
      B[i] = G[_b,_a,s] - G[_a,_b,s]
      i = i + 1
    end
  end
  return M,B
end

function toTriangularTrafo_last_sxs_up_RAND(s,G)
  d = size(G)[1]
  T = identity_matrix(QQ,d)
  M,B = linrel_mat_vec(G,s)
  for att in (1:100) # Et hätt noch emmer jot jejange
    issolv,x = can_solve_with_solution(matrix(QQ,M),B,side=:right)
    if issolv && rank(matrix(QQ,M))==d-s
      Qup = Matrix(identity_matrix(QQ,d))
      Qup[s,s+1:d] = x
      for j in (s+1:d)
        add_slice_all_3!(G,Qup[s,j],j,s)
      end
      return Matrix(Qup)*Matrix(T)
    end
    T_temp = Array(block_diagonal_matrix([identity_matrix(QQ,s-1),matrix(QQ,generic_transform_GL(d-s+1,3))]))
    matrix_tensor_congruence_3!(Array(T_temp),G)
    T = Matrix(T_temp)*Matrix(T)
    M,B = linrel_mat_vec(G,s)
  end
  print("ERROR: we tried 100 random coordinate changes and it didn't suffice at s=")
  println(s)
  return identity_matrix(QQ,d)
end

#function QQMatrixE(m::Int,i::Int,j::Int)
#  @assert i<=m
#  @assert j<=m
#  M = zero_matrix(QQ,m,m)
#  M[i,j] = one(QQ)
#  return M
#end

function diagTraf(s::Int,G1)
  d = size(G1)[1]
  P2 = (inv(root(G1[s,s,s],3))-1)*QQMatrixE(d,s,s)
  return P2 + identity_matrix(QQ,d)
end

function lowTraf(s::Int,G2)
  d = size(G2)[1]
  P3 = identity_matrix(QQ,d)
  for j in (s+1:d)
    P3[j,s] = -G2[j,s,s]
  end
  return P3
end

function tensor_learning_3(G_)
  G = copy(G_)
  d = size(G)[1]
  Q = Matrix(identity_matrix(QQ,d))
  for s in (1:d)
    if s == d
       Qup = identity_matrix(QQ,d)
    elseif s == d-2
       Qup = toTriangularTrafo_last_3x3_up(G)
    elseif s == d-1
       Qup = toTriangularTrafo_last_2x2_up(G)
    else
       Qup = toTriangularTrafo_last_sxs_up_RAND(s,G)
    end
    Qdiag = diagTraf(s,G)
    multiply_slice_all_3!(G,inv(root(G[s,s,s],3)),s)
    Qlow = lowTraf(s,G)
    for j in (s+1:d)
      add_slice_all_3!(G,Qlow[j,s],s,j)
    end
    Q = Matrix(Qlow)*Matrix(Qdiag)*Matrix(Qup)*Matrix(Q)
  end
  return Q
end