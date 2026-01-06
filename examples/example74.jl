# Example 7.4

d = 8 

Pd = QQMatrixUdmUdTNFtrafo(d)
# [1    0   0    0   0    0   0   0]
# [0    1   0    0   0    0   0   0]
# [1   -1   1    0   0    0   0   0]
# [1   -1   0    1   0    0   0   0]
# [1   -1   1   -1   1    0   0   0]
# [1   -1   1   -1   0    1   0   0]
# [1   -1   1   -1   1   -1   1   0]
# [1   -1   1   -1   1   -1   0   1]

i = 1
P1 = identity_matrix(QQ,d) + sum([QQMatrixE(d,j,2*(i-1)+1)-QQMatrixE(d,j,2*i) for j in (2*i+1:d)])
# [1    0   0   0   0   0   0   0]
# [0    1   0   0   0   0   0   0]
# [1   -1   1   0   0   0   0   0]
# [1   -1   0   1   0   0   0   0]
# [1   -1   0   0   1   0   0   0]
# [1   -1   0   0   0   1   0   0]
# [1   -1   0   0   0   0   1   0]
# [1   -1   0   0   0   0   0   1]

i = 2 
P2 = identity_matrix(QQ,d) + sum([QQMatrixE(d,j,2*(i-1)+1)-QQMatrixE(d,j,2*i) for j in (2*i+1:d)])
# [1   0   0    0   0   0   0   0]
# [0   1   0    0   0   0   0   0]
# [0   0   1    0   0   0   0   0]
# [0   0   0    1   0   0   0   0]
# [0   0   1   -1   1   0   0   0]
# [0   0   1   -1   0   1   0   0]
# [0   0   1   -1   0   0   1   0]
# [0   0   1   -1   0   0   0   1]

i = 3
P3 = identity_matrix(QQ,d) + sum([QQMatrixE(d,j,2*(i-1)+1)-QQMatrixE(d,j,2*i) for j in (2*i+1:d)])
# [1   0   0   0   0    0   0   0]
# [0   1   0   0   0    0   0   0]
# [0   0   1   0   0    0   0   0]
# [0   0   0   1   0    0   0   0]
# [0   0   0   0   1    0   0   0]
# [0   0   0   0   0    1   0   0]
# [0   0   0   0   1   -1   1   0]
# [0   0   0   0   1   -1   0   1]

@assert Pd == P3*P2*P1
NF = QQMatrixUdmUdTNF(d)
C = QQMatrixU(d) - transpose(QQMatrixU(d))
@assert Pd*C*transpose(Pd) == NF





