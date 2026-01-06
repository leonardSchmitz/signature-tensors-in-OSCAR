from sage.algebras.lie_algebras.bch import bch_iterator
import numpy as np
import re
import time
import random

def random_color():
    return "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])

channels = 2  # d
depth = 4  # k
sample_size = 2 # N
segments = 2 # m


Lbl = LieAlgebra(QQ, channels, step=depth)
basislength = len(Lbl.gens())
vars_path = channels*segments
vars_all = sample_size*vars_path
R = PolynomialRing(QQ, vars_all + basislength, 'a')
a = R.gens()
L = LieAlgebra(R, channels, step=depth)
L.inject_variables()

logSigPath1 = L.from_vector(vector(R, a[:2]+(0,0,0,0,0,0)))  # 5 = basis_length - channels
logSigPath2 = L.from_vector(vector(R, a[2:4]+(0,0,0,0,0,0)))

LynCorPath1 = list(logSigPath1.monomial_coefficients().values())
LynCorPath2 = list(logSigPath2.monomial_coefficients().values())

y = -L.from_vector(vector(R, a[8:basislength+8]))
S = sum([Z for Z in bch_iterator(y, logSigPath1)]) + sum([Z for Z in bch_iterator(y, logSigPath2)])
#S = 3/4*sum([Z for Z in bch_iterator(y, logSigPath1)]) + 1/4*sum([Z for Z in bch_iterator(y, logSigPath2)])

coef_S = vector(list(S.monomial_coefficients().values()))
res0 = 1/2*(coef_S + 2*vector(a[8:8+basislength]))
#res0 = coef_S + vector(a[8:8+basislength])
rep1 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res0[0],res0[1],a[10],a[11],a[12],a[13],a[14],a[15]])
res1 = vector([rep1(f) for f in list(res0)])
rep2 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res1[0],res1[1],res1[2],a[11],a[12],a[13],a[14],a[15]])
res2 = vector([rep2(f) for f in list(res1)])
rep3 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res2[0],res2[1],res2[2],res2[3],res2[4],a[13],a[14],a[15]])
res = vector([rep3(f) for f in list(res2)])


## ========== explicit sample "hat_origin"

origin = vector([0,0])
path1Seg1 = vector([1/2,1/2])
path2Seg1 = vector([1/2,-1/2])
#path2Seg1 = vector([1/2,1/4])  no solution????

path1 = [origin,path1Seg1]
path2 = [origin,path2Seg1]


applyBaryOnLyndCoo = R.hom([path1Seg1[0],path1Seg1[1],
                            path2Seg1[0],path2Seg1[1],
                            0,0,
                            0,0,
                            0,0,0,0,0,0,0,0]) # last number are not evaluated anyhow

evalLynCorPath1 =  [applyBaryOnLyndCoo(f) for f in LynCorPath1]
evalLynCorPath2 =  [applyBaryOnLyndCoo(f) for f in LynCorPath2]

evalLynCorBary =  [applyBaryOnLyndCoo(f) for f in res]


#### 3 segments: 2 complex and no real solution:

logSigBarySeg1 = L.from_vector(vector(R, a[:2]+(0,0,0,0,0,0)))  # 3 = basis_length - channels
logSigBarySeg2 = L.from_vector(vector(R, a[2:4]+(0,0,0,0,0,0)))
logSigBarySeg3 = L.from_vector(vector(R, a[4:6]+(0,0,0,0,0,0)))  

logSigBarySeg12 = sum([Z for Z in bch_iterator(logSigBarySeg1, logSigBarySeg2)])
logSigBary = sum([Z for Z in bch_iterator(logSigBarySeg12, logSigBarySeg3)])

LynCorBary = list(logSigBary.monomial_coefficients().values())
gens = [LynCorBary[i] - applyBaryOnLyndCoo(res[i]) for i in range(basislength)]
J = R.ideal(gens) #dim(J)=0 with 2 complex solutions and no real solution 


#### 4 segments: infinityly many complex and no real solution ? 
logSigBarySeg4 = L.from_vector(vector(R, a[6:8]+(0,0,0,0,0,0)))  
logSigBarySeg123 = sum([Z for Z in bch_iterator(logSigBarySeg12, logSigBarySeg3)])
logSigBary = sum([Z for Z in bch_iterator(logSigBarySeg123, logSigBarySeg4)])

LynCorBary = list(logSigBary.monomial_coefficients().values())
gens = [LynCorBary[i] - applyBaryOnLyndCoo(res[i]) for i in range(basislength)]
J = R.ideal(gens)
# sage: J.elimination_ideal([a[0],a[1],a[2],a[3],a[4],a[5]])
# Ideal (a6^2 + 2*a7^2) of Multivariate Polynomial Ring
# sage: J1 = R.ideal(gens+[a[7]])
# sage: J1.groebner_basis()
# [a5^2 + 1/8, a0 - 1/2, a1 - a5, a2 + 1/2, a3 + 2*a5, a4 - 1/2, a6, a7]


