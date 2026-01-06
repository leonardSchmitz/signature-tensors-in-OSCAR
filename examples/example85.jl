# Example 8.5

R, a = polynomial_ring_sig_transform(2,3)                       # d=2, m=3
TTSm = trunc_tensor_seq(R,3,3)                                  # k=3
TTSd = trunc_tensor_seq(R,3,2) 
sY = bary([sig_segment_standard_direction(TTSd,j) for j in (1,2)]) 
s_pwln = matrix_tensorSeq_congruence(a,sig_axis(TTSm)) 
I = ideal_of_entries(s_pwln - sY) + ideal(R,[a[1,1]-QQ(3,4)])   # a[1,1]=3/4 

normal_form.(a, Ref(I))
# 2×3 Matrix{QQMPolyRingElem}:
#  3//4  -1//4  0
#  1//4  -3//4  1
