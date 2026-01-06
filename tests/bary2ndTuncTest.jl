d = 4;  # path dimension
k = 2;  # truncation level
N = 100; # sample size (for this test we assume N>3)
alpha = rand((4:10),N);
R, a = polynomial_ring(QQ, [:a]);
TTSd = trunc_tensor_seq(R,k,d);
A = [R.(rand(0:100, d, alpha[i])) for i in (1:N)];
axP = [sig_axis(trunc_tensor_seq(R,k,alpha[i])) for i in (1:N)];
sX = [matrix_tensorSeq_congruence(A[i],axP[i]) for i in (1:N)];
sX1 = sX[1]; sX2 = sX[2]; sX3 = sX[3]; sX4 = sX[4];

# test Theorem 4.3
@assert sX4*bary([sX1,sX2])*sX3 == bary([sX4*sX1*sX3,sX4*sX2*sX3])
@assert bary([sX1,sX2,sX3])*sX4 == bary([sX1*sX4,sX2*sX4,sX3*sX4])
@assert sX4*bary([sX1,sX2,sX3]) == bary([sX4*sX1,sX4*sX2,sX4*sX3])
@assert bary([inv(sX1),inv(sX2),inv(sX3),inv(sX4)]) == inv(bary([sX1,sX2,sX3,sX4]))

# basic arithmetic tests
@assert exp(log(sX1)) == sX1
@assert inv(inv(sX1)) == sX1
@assert sX1*inv(sX1) == one(TTSd) 
@assert inv(sX1)*sX1 == one(TTSd) 
@assert sX1^3 == sX1*sX1*sX1
@assert exp(log(sX1)+log(sX1)) == sX1^2

# test Proposition 4.4
@assert bary([sX1,sX2]) == sX1*exp(QQ(1,2)*log(inv(sX1)*sX2))
@assert bary([one(TTSd),sX2]) == exp(QQ(1,2)*log(sX2))
@assert bary_2samples(sX1,sX2) == bary([sX1,sX2])

# test Proposition 4.5 / Example 4.6
@assert bary([sX1^2,one(TTSd)])==sX1
@assert bary([sX1,inv(sX1)])==one(TTSd)
@assert bary([sX1^6,sX1^2,sX1])==sX1^3

# test Proposition 4.8
# i)
@assert bary([sX1,sX1])==sX1
@assert bary([sX1,sX1,sX1])==sX1
@assert bary([sX1,sX1,sX1,sX1])==sX1
# ii)
@assert bary([sX1,sX2,sX3,sX4])==bary([sX4,sX1,sX2,sX3])
@assert bary([sX1,sX2,sX3,sX4])==bary([sX4,sX3,sX2,sX1])
# iii)
@assert bary_Nminus1_samples_fixed_inverse([sX1,sX2,sX3],bary(sX[1:4])) == sX4

# test Remark 4.10
@assert bary(sX) == bary_2nd_trunc(sX)
# test Proposition 4.11
@assert bary(sX) == bary_2nd_trunc_closedform(sX)

##################### timinings for README.md #############################
#  @time bary_2nd_trunc(sX)
#    # 0.016592 seconds (129.30 k allocations: 7.413 MiB)
#
#  @time bary_2nd_trunc_closedform(sX);
#    # 0.365100 seconds (1.81 M allocations: 117.283 MiB, 30.82% gc time)
#
#  @time bary(sX)
#    # 1.938035 seconds (29.74 M allocations: 1.906 GiB, 35.24% gc time)
############################################################################

