# Example 6.4

###################
# B_22(1,1) = 1   #
###################
R, a = polynomial_ring_sig_transform(2,1)                          # d = 2, m = 1
TTSm = trunc_tensor_seq(R,2,1)                                     # k = 2
sigX1 = matrix_tensorSeq_congruence([QQ(1,2);1],sig_axis(TTSm))
sigX2 = matrix_tensorSeq_congruence([1;QQ(1,2)],sig_axis(TTSm))
sY = bary([sigX1,sigX2])
s_pwln = matrix_tensorSeq_congruence(a,sig_axis(TTSm))             # Ansatz
I = ideal_of_entries(s_pwln - sY)                                  # relations 
groebner_basis(I)                                                  # solve 
# Gröbner basis with elements
#   1: 4*a₂₁ - 3
#   2: 4*a₁₁ - 3
# with respect to the ordering
#   degrevlex([a₁₁, a₂₁])

###################
# B_22(1,1,1) = 1 #
###################
R, a = polynomial_ring_sig_transform(2,1)                          # d = 2, m = 1
TTSm = trunc_tensor_seq(R,2,1)                                     # k = 2
sigX1 = matrix_tensorSeq_congruence([QQ(1,2);1],sig_axis(TTSm))
sigX2 = matrix_tensorSeq_congruence([1;QQ(1,2)],sig_axis(TTSm))
sigX3 = matrix_tensorSeq_congruence([QQ(3,4);0],sig_axis(TTSm))
sY = bary([sigX1,sigX2,sigX3])
s_pwln = matrix_tensorSeq_congruence(a,sig_axis(TTSm))             # Ansatz
I = ideal_of_entries(s_pwln - sY)                                  # relations 
groebner_basis(I)                                                  # solve 
# Gröbner basis with elements
#   1: 2*a₂₁ - 1
#   2: 4*a₁₁ - 3
# with respect to the ordering
#   degrevlex([a₁₁, a₂₁])

#################
# B_22(2,1) = 2 #
#################
R, a = polynomial_ring_sig_transform(2,1)                          # d = 2, m = 1
TTS1 = trunc_tensor_seq(R,2,1)                                     # k = 2
TTS2 = trunc_tensor_seq(R,2,2)                                     # k = 2
sigX1 = matrix_tensorSeq_congruence([QQ(1,2) QQ(1,4);1 0],sig_axis(TTS2))
sigX2 = matrix_tensorSeq_congruence([1;QQ(1,2)],sig_axis(TTS1))
sY = bary([sigX1,sigX2])
s_pwln = matrix_tensorSeq_congruence(a,sig_axis(TTS1))             # Ansatz
I = ideal_of_entries(s_pwln - sY)                                  # relations 
dim(I)                                                  # solve 
#-infinity
R, a = polynomial_ring_sig_transform(2,2)                          # d = 2, m = 2
TTS1 = trunc_tensor_seq(R,2,1)                                     # k = 2
TTS2 = trunc_tensor_seq(R,2,2)                                     # k = 2
sigX1 = matrix_tensorSeq_congruence([QQ(1,2) QQ(1,4);1 0],sig_axis(TTS2))
sigX2 = matrix_tensorSeq_congruence([1;QQ(1,2)],sig_axis(TTS1))
sY = bary([sigX1,sigX2])
s_pwln = matrix_tensorSeq_congruence(a,sig_axis(TTS2))             # Ansatz
I = ideal_of_entries(s_pwln - sY) + ideal(R,[a[1,1]-QQ(1,2)])      # relations  
groebner_basis(I)                                                  # solve 
# Gröbner basis with elements
#   1: 7*a₂₂ - 2
#   2: 2*a₁₂ - 1
#   3: 4*a₂₁ + 4*a₂₂ - 3
#   4: 8*a₁₁ - 3
# with respect to the ordering
#   degrevlex([a₁₁, a₂₁, a₁₂, a₂₂])
normal_form.(a, Ref(I))
# 2×2 Matrix{QQMPolyRingElem}:
#  1//2  3//8
#  4//7  5//28

#################
# B_23(1,1) = 3 #
#################
R, a = polynomial_ring_sig_transform(2,2);                 # d=2, m=2
TTSm = trunc_tensor_seq(R,3,2);                            # k=3
TTSd = trunc_tensor_seq(R,3,2);
TTS1 = trunc_tensor_seq(R,3,1);
s1 = matrix_tensorSeq_congruence([1;QQ(1,2)],sig_axis(TTS1))
s2 = matrix_tensorSeq_congruence([1;-QQ(1,2)],sig_axis(TTS1))
sY = bary([s1,s2])
s_pwln = matrix_tensorSeq_congruence(a,sig_axis(TTSm))
I = ideal_of_entries(s_pwln - sY);
dim(I)
#-infinity
R, a = polynomial_ring_sig_transform(2,3);                 # d=2, m=3
TTSm = trunc_tensor_seq(R,3,3);                            # k=3
TTSd = trunc_tensor_seq(R,3,2);
TTS1 = trunc_tensor_seq(R,3,1);
s1 = matrix_tensorSeq_congruence([1;QQ(1,2)],sig_axis(TTS1))
s2 = matrix_tensorSeq_congruence([1;-QQ(1,2)],sig_axis(TTS1))
sY = bary([s1,s2])
s_pwln = matrix_tensorSeq_congruence(a,sig_axis(TTSm))
I = ideal_of_entries(s_pwln - sY) + ideal(R,[a[2,1]-QQ(1,4)]);
normal_form.(a, Ref(I))
# 2×3 Matrix{QQMPolyRingElem}:
#  1     -1    1
#  1//4  1//4  -1//2
