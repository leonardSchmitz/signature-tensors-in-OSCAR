#####################################
# truncated tensor sequences 
#####################################

#construction
d = 2; k = 3; n = 3;                            # d dimension, k truncation level, n core dimension 
R, a = polynomial_ring_sig_transform(d,n);      # base ring (QQ would be desireble) 
TTSd = trunc_tensor_seq(R,k,d);                 # struct of truncated tensor sequences (a group)

#group elements:
one(TTSd)                                       # 1 wrt multiplication
sig_axis(TTSd)                                  # signature of n-dim axis path
sig_mono(TTSd)                                  # signature of n-dim monomial curve
sig_segment_standard_direction(TTSd,1)          # signature of segment in direction e_1

#------------------------------------
# what I would like to do (with 'type union' I had problems)
# TTSn = trunc_tensor_seq(QQ,k,d);  
#matrix parts:
QQMatrix(one(TTSd))                             # types are not correct here 
QQMatrix(sig_axis(TTSd))
QQMatrix(sig_mono(TTSd))
QQMatrix(sig_segment_standard_direction(TTSd,1))
#------------------------------------

# Chens identity 
S = sig_axis(TTSd)
@assert S == prod([sig_segment_standard_direction(TTSd,i) for i in (1:d)])
@assert inv(S)*S == one(TTSd)

# equivariance
TTSn = trunc_tensor_seq(R,k,n);                 # struct of truncated tensor sequences (a group)
a
S_axis_n = sig_axis(TTSn)
S_pwln = matrix_tensorSeq_congruence(a,S_axis_n)
@assert Matrix(S_pwln) == a*Matrix(S_axis_n)*transpose(a)
@assert parent(S_pwln) != parent(sig_axis(TTSn))
@assert parent(S_pwln) == parent(one(TTSd))


#####################################
# Learning path from signature tensors 
# Pfeffer, Seigal, Sturmfels '19
#####################################

A = a + rand(0:100, d, n) - a                  # random transform 
S_data = matrix_tensorSeq_congruence(A,S_axis_n)
temp = (tensor_sequence(S_pwln - S_data))[2:k+1];
gens4learn = vcat([vec(temp[i]) for i in (1:length(temp))]...);
I = ideal(R,gens4learn);
dim(I)


#####################################
# Lie group barycenter
# joint /w Amendola 
#####################################

sample = [one(TTSd),sig_axis(TTSd)];
S_bary = bary(sample)                           # group barycenter of 2 samples 

temp = (tensor_sequence(S_pwln - S_bary))[2:k+1];
gens4learn = vcat([vec(temp[i]) for i in (1:length(temp))]...);
I = ideal(R,gens4learn);
dim(I)       


#####################################
# Free turncated signature algebra (multivariate)
#####################################

k = 3; N = 2                                                 # truncation level,  number of samples
T,s = free_trunc_sig_alg_multiv(k,N,true);                   # last Bool includes variables for the barycenter (syn. sugar)

# algebra elements
one(T)
f = [free_sig_from_sample(i,T) for i in (1:N)]
free_bary_2sampl = f[1]*exp(QQ(1,2)*log(inv(f[1])*f[2]))     # closed form special for N=2 samples using geodesics

bary_defining_polynomial_system(k,N)                         # for N>2 we solve the system in the y's 

@assert evaluate(free_bary_2sampl,sample) == bary(sample)  


