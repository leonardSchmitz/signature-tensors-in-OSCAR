d = 5;  # path dimension
k = 2;  # truncation level
N = 6; # sample size
#alpha = rand((1:10),N);
alpha = [1,1,1,1,1,1];
R, a = polynomial_ring_sig_transform(d,1);                          # d = 2, m = 1
TTSm = trunc_tensor_seq(R,k,1); 
TTSd = trunc_tensor_seq(R,k,d);
A = [R.(rand(0:100, d, alpha[i])) for i in (1:N)];
axP = [sig_axis(trunc_tensor_seq(R,k,alpha[i])) for i in (1:N)];
sX = [matrix_tensorSeq_congruence(A[i],axP[i]) for i in (1:N)];
sY = bary(sX);
s_pwln = matrix_tensorSeq_congruence(a,sig_axis(TTSm));             # Ansatz
I = ideal_of_entries(s_pwln - sY);                                  # relations 
@assert dim(I) != -1

for alpha in [[2,2],[3,1,2],[2,2,1,2],[2,1,1,2], 
          [1,2,1,2],[1,1,1],[1,1,1,3],[1,1,1,1,1]]
  print('*')
  N = length(alpha)
  sum_alpha = sum(alpha)
  d = sum_alpha
  odds = length([ai for ai in alpha if isodd(ai)])
  if odds == sum_alpha
    Bd2 = sum_alpha
  else 
    Bd2 = sum_alpha - odds + 1 
  end 
  R, a = polynomial_ring_sig_transform(d,Bd2);                          # d = 2, m = Bd2
  TTSBd2 = trunc_tensor_seq(R,k,Bd2);
  TTSd = trunc_tensor_seq(R,k,d);
  A = [R.(rand(0:100, d, alpha[i])) for i in (1:N)];
  axP = [sig_axis(trunc_tensor_seq(R,k,alpha[i])) for i in (1:N)];
  sX = [matrix_tensorSeq_congruence(A[i],axP[i]) for i in (1:N)];
  sY = bary(sX);
  s_pwln = matrix_tensorSeq_congruence(a,sig_axis(TTSBd2));             # Ansatz
  I = ideal_of_entries(s_pwln - sY);                                    # relations 
  @assert is_one(I) != true
end

# see also
# tests/W_alpha_tests.jl 
# and 
# tests/signature_matrix_tests.jl
