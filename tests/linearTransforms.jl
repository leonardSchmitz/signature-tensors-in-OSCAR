############################################################
# test equivariance of linear transforms and barycenters   #
# used in the proof of Proposition 6.5                     #
############################################################
m = 10; d =3;
m1 = 5;
m2 = 4;
k = 3;
R, a = polynomial_ring_sig_transform(d,m);
TTSm = trunc_tensor_seq(R,k,m);
Sis = [sig_segment_standard_direction(TTSm,i) for i in (1:m)];
S1 = prod([Sis[i] for i in (1:m1)]);
S2 = prod([Sis[i] for i in (m1+1:m1+m2)]);
BS1S2 = bary([S1,S2]);
aS1 = matrix_tensorSeq_congruence(a,S1);
aS2 = matrix_tensorSeq_congruence(a,S2);
BaS1aS2 = bary([aS1,aS2]);
aBS1S2 = matrix_tensorSeq_congruence(a,BS1S2);
@assert BaS1aS2==aBS1S2

###########################################################
# test equivariance of evaluations and linear transforms  #
# see Lemma 3.9                                           #
###########################################################
k = 4; N = 2                                                 # truncation level,  number of samples
T,s = free_trunc_sig_alg_multiv(k,N);
for f in [s[1,2]*s[1,1] - s[2,2]*s[2,1] + s[3,1]*s[1,2] - s[2,1] + s[4,1] + s[4,2], 
          s[1,2] + s[2,2]*s[2,2] - s[1,1] + s[3,1] + QQ(1,2)*s[1,1]*s[2,2]]
  d = 5; n = 4;
  R, a = polynomial_ring_sig_transform(d,n);      # base ring (QQ would be desireble) 
  TTSn = trunc_tensor_seq(R,k,n);                 # struct of truncated tensor sequences (a group)
  TTSd = trunc_tensor_seq(R,k,d);
  S_axis_n = sig_axis(TTSn);
  #S_axis_n = bary([sig_axis(TTSn),one(TTSn)]);
  S_mono_n = sig_mono(TTSn);
  S_pwln = matrix_tensorSeq_congruence(a,S_axis_n);
  S_poly = matrix_tensorSeq_congruence(a,S_mono_n);
  sample = [S_axis_n,S_mono_n];
  Asample = [S_pwln,S_poly];
  fs = evaluate(f,sample);
  Afs = matrix_tensorSeq_congruence(a,fs);
  fAs = evaluate(f,Asample);
  @assert Afs == fAs
end 


