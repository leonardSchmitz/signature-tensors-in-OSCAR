R, a = polynomial_ring_sig_transform(2,1)                          # d = 2, m = 1
TTSm = trunc_tensor_seq(R,2,1)                                     # k = 2
sigX1 = matrix_tensorSeq_congruence([QQ(1,2);1],sig_axis(TTSm))
sigX2 = matrix_tensorSeq_congruence([1;QQ(1,2)],sig_axis(TTSm))
sY = bary([sigX1,sigX2])
s_pwln = matrix_tensorSeq_congruence(a,sig_axis(TTSm))             # Ansatz
I = ideal_of_entries(s_pwln - sY)                                  # relations 

groebner_basis(I)                                           # solve 
# Gröbner basis with elements
#   1: 4*a₂₁ - 3
#   2: 4*a₁₁ - 3
# with respect to the ordering
#   degrevlex([a₁₁, a₂₁])
