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

R, a = polynomial_ring_sig_transform(2,3);                 # d=2, m=3
TTSm = trunc_tensor_seq(R,3,3);                            # k=3
TTSd = trunc_tensor_seq(R,3,2);
TTS1 = trunc_tensor_seq(R,3,1);
A = Array(QQ[1 1; QQ(1,2) -QQ(1,2)])
s1 = matrix_tensorSeq_congruence(A,sig_segment_standard_direction(TTSd,1))
s2 = matrix_tensorSeq_congruence(A,sig_segment_standard_direction(TTSd,2))
sY = bary([s1,s2])
AsEl = matrix_tensorSeq_congruence(A,bary([sig_segment_standard_direction(TTSd,j) for j in (1,2)]))
@assert sY == AsEl


