# Tests for 
# - Lemma 7.3
# - Corollary 7.5
# - Lemma 7.6
# - Lemma 7.7
# - Lemma 7.8
# - Corollary 7.10

#####################
# - Lemma 7.3       #
#####################

## test normal forms and transformations
for d in (2:15)
  Pd = QQMatrixUdmUdTNFtrafo(d)
  NF = QQMatrixUdmUdTNF(d)
  C = QQMatrixU(d) - transpose(QQMatrixU(d))
  @assert Pd*C*transpose(Pd)==NF
end

#####################
# - Corollary 7.5   #
#####################

for m in (2:2:15)
  C = QQMatrixU(m) - transpose(QQMatrixU(m))
  Q = QQMatrixQ(m)
  @assert inv(C)==Q*C*Q
end

#####################
# - Lemma 7.6       #
#####################

## test normal forms and transformations
for d in (2:15)
  Pd = QQMatrixAxisNFtrafo(d)
  NF = QQMatrixSigAxisNF(d)
  C = QQMatrixSigAxis(d)
  @assert Pd*C*transpose(Pd)==QQ(1,2)*NF # we obmit sqrt(2) in the transformation
end

#####################
# - Lemma 7.7       #
#####################

# test definition W_m=bary(sig(Ax^{m,1},...,sig(Ax^{m,N})) 
for a in [[2,3,4],[1,2,3],[4,4,4],[3,1,3],[3,3,5],[3,4,4],[1,1,1]]
  @assert QQMatrixBaryAxisAxisAxis(a[1],a[2],a[3]) == QQMatrixBaryNAxis(a) 
  N = length(a)
  B1 = QQMatrixU(a[1]) - transpose(QQMatrixU(a[1]))
  B2 = QQMatrixU(a[2]) - transpose(QQMatrixU(a[2]))
  B3 = QQMatrixU(a[3]) - transpose(QQMatrixU(a[3]))
  B = QQ(1,2*N)*block_diagonal_matrix([B1,B2,B3]) + QQ(1,2*N*N)*ones_matrix(QQ,sum(a),sum(a))
  @assert QQMatrixBaryNAxis(a) == B
end 
for a in [[2,3],[1,2],[4,4],[3,1],[3,3],[3,4],[1,1]]
  @assert QQMatrixBaryAxisAxis(a[1],a[2]) == QQMatrixBaryNAxis(a) 
  N = length(a)
  B1 = QQMatrixU(a[1]) - transpose(QQMatrixU(a[1]))
  B2 = QQMatrixU(a[2]) - transpose(QQMatrixU(a[2]))
  B = QQ(1,2*N)*block_diagonal_matrix([B1,B2]) + QQ(1,2*N*N)*ones_matrix(QQ,sum(a),sum(a))
  @assert QQMatrixBaryNAxis(a) == B
end 
for d in (2:8) 
  R, a = polynomial_ring_sig_transform(d,1);      
  TTSd = trunc_tensor_seq(R,2,d);
  A = sig_axis(TTSd)
  @assert QQMatrix(A) == QQMatrixSigAxis(d)
  @assert QQMatrix(bary([A,one(TTSd)])) ==  QQMatrixBaryAxis1(d)
end 

#####################
# - Lemma 7.9 and Corollary 7.11
#####################

for a in [[8,2,4,2],[2,2],[8,4,2],[8,4,2],[6,6,10,2],[2,8,8,2],[2,2,2,2,2],[8],
          [8,2,1,2],[1,2],[3,4,2],[7,3,5],[6,6,10,2],[2,8,7,2],[2,2,2,2,1],[9],
          [8,2,1,2,12],[1,1,1],[3,12,2],[1,1,1,3,5],[6,6,1,2],[2,1,7,2],[1,1,1,1,1],[1]]
  print('*')
  m = sum(a)
  N = length(a)
  W = QQMatrixBaryNAxis(a);
  @assert transpose(W) == -W + QQ(1,N^2)*ones_matrix(QQ,m,m)
  W_sym = QQ(1,2)*(W + transpose(W))
  @assert W_sym == QQ(1,2*N^2)*ones_matrix(QQ,m,m)
  W_skew = QQ(1,2)*(W - transpose(W))
  @assert W_skew == QQ(1,2*N)*block_diagonal_matrix([QQMatrixU(ai)-transpose(QQMatrixU(ai)) for ai in a])
  odds = length([ai for ai in a if isodd(ai)])
  if odds == 0
    Q = QQMatrixQ(m)
    @assert inv(W) == 4*N^2*Q*W*Q
  else
    @assert rank(W) == m - odds + 1
  end 
end

## 2 segments (older tests)
for d in (2:14)
 for m1 in (1:d-1)
   m2 = d-m1; k = 2;
   R, a = polynomial_ring_sig_transform(d,m1+m2);
   TTSd = trunc_tensor_seq(R,k,d);
   Sis = [sig_segment_standard_direction(TTSd,i) for i in (1:d)];
   S1 = prod([Sis[i] for i in (1:m1)]);
   S2 = prod([Sis[i] for i in (m1+1:d)]);
   M = QQMatrix(bary([S1,S2]))
   if is_odd(m1) && is_odd(m2)
     @assert rank(M)==d-1
   else
     @assert rank(M)==d
   end
 end
end

## test 3 segments
m = 12; d =2;
m1 = 3;
m2 = 3;
m3 = 3;
k = 3;
R, a = polynomial_ring_sig_transform(d,m);
TTSm = trunc_tensor_seq(R,k,m);
Sis = [sig_segment_standard_direction(TTSm,i) for i in (1:m)];
S1 = prod([Sis[i] for i in (1:m1)]);
S2 = prod([Sis[i] for i in (m1+1:m1+m2)]);
S3 = prod([Sis[i] for i in (m1+m2+1:m1+m2+m3)]);
BS1S2S3 = bary([S1,S2,S3]);
@assert rank(QQMatrix(BS1S2S3))== 7


