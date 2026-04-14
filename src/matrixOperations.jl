export 
  QQMatrixU, # elemantary matrix constructors
  QQMatrixE, 
  QQMatrixQ, 
  QQMatrixSigAxis, # signature matrices
  QQMatrixBaryAxis1, 
  QQMatrixBaryAxisAxis, 
  QQMatrixBaryAxisAxisAxis, 
  QQMatrixBaryNAxis, 
  QQMatrixAxisNFtrafo, # transformation matrices
  QQMatrixUdmUdTNFtrafo,
  QQMatrixWNFP1m1oddHelpertrafo, 
  QQMatrixSigAxisNF, # normal forms
  QQMatrixSigNF, 
  QQMatrixUdmUdTNF, 
  generic_transform, 
  generic_transform_GL
  #QQMatrixUmUT_MPinverse

######################################
# elemantary matrix constructors     # 
######################################

# upper triangular matrix
function QQMatrixU(d::Int)
  if d == 1 
    return QQ[0;]
  else  
    return Ud = strictly_upper_triangular_matrix([one(QQ) for i in (1:div((d-1)*d,2))]); 
  end 
end


function QQMatrixE(m::Int,i::Int,j::Int)
  @assert i<=m
  @assert j<=m
  M = zero_matrix(QQ,m,m)
  M[i,j] = one(QQ)
  return M
end

function QQMatrixQ(d::Int)
  return diagonal_matrix([(-one(QQ))^(i+1) for i in (1:d)])
end

function generic_transform_GL(_d,_max_int=10)
   A = QQ.(rand(-1*_max_int:_max_int, _d, _d))
   if is_invertible( matrix(QQ,A))
     return A
   else
     return generic_transform_GL(_d,_max_int)
   end
end

function generic_transform(_d,_m,_max_int=10)
   return QQ.(rand(-1*_max_int:_max_int, _d, _m))
end



############################################
# signature matrices                       #
############################################


function QQMatrixSigAxis(d::Int)
  Ud = strictly_upper_triangular_matrix([one(QQ) for i in (1:div((d-1)*d,2))]);
  return QQ(1,2)*identity_matrix(QQ,d)+Ud
end

function QQMatrixBaryAxis1(d::Int)
  Ud = strictly_upper_triangular_matrix([one(QQ) for i in (1:div((d-1)*d,2))]);
  return QQ(1,4)*(Ud - transpose(Ud)) + QQ(1,8)*ones_matrix(QQ,d,d)
end

function QQMatrixBaryAxisAxis(m1::Int,m2::Int)
  d = m1 + m2
  tm1 = QQMatrixU(m1) - transpose(QQMatrixU(m1))
  tm2 = QQMatrixU(m2) - transpose(QQMatrixU(m2))
  return QQ(1,4)*block_diagonal_matrix([tm1,tm2]) + QQ(1,8)*ones_matrix(QQ,d,d)
end

function QQMatrixBaryAxisAxisAxis(m1::Int,m2::Int,m3::Int)
  d = m1 + m2 + m3
  N = 3
  tm1 = QQMatrixU(m1) - transpose(QQMatrixU(m1))
  tm2 = QQMatrixU(m2) - transpose(QQMatrixU(m2))
  tm3 = QQMatrixU(m3) - transpose(QQMatrixU(m3))
  return QQ(1,2*N)*block_diagonal_matrix([tm1,tm2,tm3]) + QQ(1,2*N*N)*ones_matrix(QQ,d,d)
end

function QQMatrixBaryNAxis(m::Vector{Int})
  d = sum(m)
  N = length(m)
  ts = [QQMatrixU(mi) - transpose(QQMatrixU(mi)) for mi in m]
  return QQ(1,2*N)*block_diagonal_matrix(ts) + QQ(1,2*N*N)*ones_matrix(QQ,d,d)
end

######################################
# hardcoded normal forms             #
######################################

function QQMatrixSigAxisNF(d::Int)
  H2m1 = QQ[0 1;-1 0]
  if is_odd(d)
    res = QQ[1;]
  else
    res = QQ[0 -1 ; 1 1]
  end
  d2 = Int(floor((d-1)/2))
  for i in (1:d2)
    res = block_diagonal_matrix([res,H2m1])
  end
  return res
end

function QQMatrixUdmUdTNF(d::Int)
  H2m1 = QQ[0 1;-1 0]
  if d == 2 
    return H2m1 
  end 
  if is_odd(d)
    res = QQ[0;]
  else
    res = H2m1
  end
  d2 = Int(floor((d-1)/2))
  for i in (1:d2)
    res = block_diagonal_matrix([H2m1,res])
  end
  return res
end

function multiply_row_column(a,s,i)
   return multiply_row(multiply_column(a,s,i),s,i)
end

function add_row_column(a,s,i,j)
   return add_row(add_column(a,s,i,j),s,i,j)
end

function QQMatrixSigNF(S)
   d = size(S)[1]
   for i in (1:d)
     S = multiply_row_column(S,1/sqrt(S[i,i]),i)
   end 
   Q = identity_matrix(QQ,d) - QQMatrixE(d,2,1)
   S = transpose(Q)*S*Q
   S = multiply_row_column(S,1/S[2,1],1)
   for i in (3:d)
     S = add_row_column(S,-S[i,1],2,i)
   end 
   S = multiply_row_column(S,1/sqrt(S[3,3]),3)
   S = multiply_row_column(S,1/sqrt(S[4,4]),4)
   S = add_row_column(S,-1,4,3)
   S = add_row_column(S,-S[2,3],1,3)
   S = add_row_column(S,-1,4,2)
   for i in (3:d)
     S = add_row_column(S,S[i,2],1,i)
   end 
   return S
end


###########################
# transformation matrices #
###########################

function QQMatrixBaryAxis1NFtrafo(d::Int)
  dp12 = Int(floor((d+1)/2))
  t1 = identity_matrix(QQ,d)
  t2 = sum([-1*QQMatrixE(d,j,2*i) for i in (1:dp12) for j in (2*i+1:d)])
  t3 = sum([QQMatrixE(d,j,2*i-1) for i in (1:dp12) for j in (2*i-1:d)])
  return t1+t2+t3
end

function QQMatrixAxisNFtrafo(d::Int)
  if d == 2
    return  QQMatrixE(d,2,1) +  QQMatrixE(d,1,2) -  QQMatrixE(d,1,1)
  end 
  d2 = Int(floor((d-1)/2))
  t1 = identity_matrix(QQ,d)
  t2 = sum([(-1)*QQMatrixE(d,i,1) for i in (2:d)])
  t3 = sum([QQMatrixE(d,j,d-2*(i-1))-QQMatrixE(d,j,d-2*(i-1)-1) for i in (1:d2) for j in (1:d-2*i)])
  t4 = t1+t2+t3
  if is_even(d)
    return swap_rows(t4,1,2)
  end
  return t4
end

function QQMatrixUdmUdTNFtrafo(d::Int)
  if d == 2
    return identity_matrix(QQ,d)
  end 
  d2 = Int(floor((d-1)/2))
  t1 = identity_matrix(QQ,d)
  t2 = sum([QQMatrixE(d,j,2*(i-1)+1)-QQMatrixE(d,j,2*i) for i in (1:d2) for j in (2*i+1:d)])
  return t1+t2
end

function QQMatrixWNFP1m1oddHelpertrafo(m1::Int,m2::Int)
  d = m1 + m2
  return prod([prod([add_row(identity_matrix(QQ,m1+m2),(-1)^(m1-j)*one(QQ),m1-j,m1+m2-i) for i in (0:m2-1)]) for j in (0:m1-1)])
end


