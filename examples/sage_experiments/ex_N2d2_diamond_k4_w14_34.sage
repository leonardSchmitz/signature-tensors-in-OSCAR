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

logSigPath1Seg1 = L.from_vector(vector(R, a[:2]+(0,0,0,0,0,0)))  # 5 = basis_length - channels
logSigPath1Seg2 = L.from_vector(vector(R, a[2:4]+(0,0,0,0,0,0)))
logSigPath2Seg1 = L.from_vector(vector(R, a[4:6]+(0,0,0,0,0,0)))
logSigPath2Seg2 = L.from_vector(vector(R, a[6:8]+(0,0,0,0,0,0)))  # 8 = vars_all, (d,k,N)=(2,3,2)

logSigPath1 = sum([Z for Z in bch_iterator(logSigPath1Seg1, logSigPath1Seg2)])
logSigPath2 = sum([Z for Z in bch_iterator(logSigPath2Seg1, logSigPath2Seg2)])

LynCorPath1 = list(logSigPath1.monomial_coefficients().values())
LynCorPath2 = list(logSigPath2.monomial_coefficients().values())

y = -L.from_vector(vector(R, a[8:basislength+8]))
#S = sum([Z for Z in bch_iterator(y, logSigPath1)]) + sum([Z for Z in bch_iterator(y, logSigPath2)])
#S = 3/4*sum([Z for Z in bch_iterator(y, logSigPath1)]) + 1/4*sum([Z for Z in bch_iterator(y, logSigPath2)])
S = 1/4*sum([Z for Z in bch_iterator(y, logSigPath1)]) + 3/4*sum([Z for Z in bch_iterator(y, logSigPath2)])

coef_S = vector(list(S.monomial_coefficients().values()))
#res0 = 1/2*(coef_S + 2*vector(a[8:8+basislength]))
res0 = coef_S + vector(a[8:8+basislength])
rep1 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res0[0],res0[1],a[10],a[11],a[12],a[13],a[14],a[15]])
res1 = vector([rep1(f) for f in list(res0)])
rep2 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res0[0],res0[1],res1[2],a[11],a[12],a[13],a[14],a[15]])
res2 = vector([rep2(f) for f in list(res1)])
rep3 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res0[0],res0[1],res1[2],res2[3],res2[4],a[13],a[14],a[15]])
res = vector([rep2(f) for f in list(res1)])


## ========== explicit sample "diamond" 

origin = vector([0,0])
path1Seg1 = vector([1/2,1/2])
path1Seg2 = vector([1/2,-1/2])
path2Seg1 = vector([1/2,-1/2])
path2Seg2 = vector([1/2,1/2])

path1 = [origin,path1Seg1,path1Seg1+path1Seg2]
path2 = [origin,path2Seg1,path2Seg1+path2Seg2]


# ========== evaluation homomorphism

applyBaryOnLyndCoo = R.hom([path1Seg1[0],path1Seg1[1],
                            path1Seg2[0],path1Seg2[1],
                            path2Seg1[0],path2Seg1[1],
                            path2Seg2[0],path2Seg2[1],
                            0,0,0,0,0,0,0,0]) # last number are not evaluated anyhow

evalLynCorPath1 =  [applyBaryOnLyndCoo(f) for f in LynCorPath1] 
evalLynCorPath2 =  [applyBaryOnLyndCoo(f) for f in LynCorPath2]

evalLynCorBary =  [applyBaryOnLyndCoo(f) for f in res]


## recovery of sample: 
gens_sample1 = [LynCorBary[i] - evalLynCorPath1[i] for i in range(basislength)]
J_sample1 = R.ideal(gens_sample1)   
# LynCorPath1


## four segments:

logSigBarySeg1 = L.from_vector(vector(R, a[:2]+(0,0,0,0,0,0)))  # 3 = basis_length - channels
logSigBarySeg2 = L.from_vector(vector(R, a[2:4]+(0,0,0,0,0,0)))
logSigBarySeg3 = L.from_vector(vector(R, a[4:6]+(0,0,0,0,0,0)))
logSigBarySeg4 = L.from_vector(vector(R, a[6:8]+(0,0,0,0,0,0)))

logSigBarySeg12 = sum([Z for Z in bch_iterator(logSigBarySeg1, logSigBarySeg2)])
logSigBarySeg123 = sum([Z for Z in bch_iterator(logSigBarySeg12, logSigBarySeg3)])
logSigBary = sum([Z for Z in bch_iterator(logSigBarySeg123, logSigBarySeg4)])

LynCorBary = list(logSigBary.monomial_coefficients().values())
gens = [LynCorBary[i] - applyBaryOnLyndCoo(res[i]) for i in range(basislength)]
J = R.ideal(gens)   

#sage: J.elimination_ideal([a[0],a[1],a[2],a[3],a[4],a[5],a[6]])
#Ideal (960*a7^4 - 1408*a7^3 + 560*a7^2 - 16*a7 - 21) of Multiv
#sage: factor(960*a[7]^4 - 1408*a[7]^3 + 560*a[7]^2 - 16*a[7] - 21)
#(24*a7^2 - 16*a7 - 3) * (40*a7^2 - 32*a7 + 7)

fakt_g = 24*a[7]^2 - 16*a[7] - 3
J1 = R.ideal(gens+[fakt_g])
#sage: J1.groebner_basis()
#[a7^2 + 2/3*a7 - 1/8, a0 + a7 + 5/4, a1 - a7, a2 - 5/2*a7 - 11/4, 
#a3 + a7 - 1/4, a4 + 5/2*a7 + 13/4, a5 + a7 + 1/4, a6 - a7 - 11/4]
nfs = [J1.reduce(a[i]) for i in range(basislength)]
#sage: nfs
#[-a7 - 5/4,
# a7,
# 5/2*a7 + 11/4,
# -a7 + 1/4,
# -5/2*a7 - 13/4,
# -a7 - 1/4,
# a7 + 11/4,
# a7]
root1 = 1/3 + sqrt(17/2)/6
applyBaryOnLyndCoo = R.hom([0,0,0,0,0,0,0,root1,
                            0,0,0,0,0,0,0,0]) # last number are not evaluated anyhow
a_nfs = [applyBaryOnLyndCoo(f) for f in nfs]
bary_Seg1 = vector([a_nfs[0],a_nfs[1]])
bary_Seg2 = vector([a_nfs[2],a_nfs[3]])
bary_Seg3 = vector([a_nfs[4],a_nfs[5]])
bary_Seg4 = vector([a_nfs[6],a_nfs[7]])

bary = [origin,bary_Seg1,
               bary_Seg1+bary_Seg2,
               bary_Seg1+bary_Seg2+bary_Seg3,
               bary_Seg1+bary_Seg2+bary_Seg3+bary_Seg4]

P = list_plot(path1, plotjoined=True, color='blue')
P += list_plot(path2, plotjoined=True, color='red')

P += list_plot(bary, plotjoined=True, color = 'purple')

root2 = 1/3 - sqrt(17/2)/6
applyBaryOnLyndCoo = R.hom([0,0,0,0,0,0,0,root2,
                            0,0,0,0,0,0,0,0]) # last number are not evaluated anyhow
a_nfs = [applyBaryOnLyndCoo(f) for f in nfs]
bary_Seg1 = vector([a_nfs[0],a_nfs[1]])
bary_Seg2 = vector([a_nfs[2],a_nfs[3]])
bary_Seg3 = vector([a_nfs[4],a_nfs[5]])
bary_Seg4 = vector([a_nfs[6],a_nfs[7]])

bary = [origin,bary_Seg1,
               bary_Seg1+bary_Seg2,
               bary_Seg1+bary_Seg2+bary_Seg3,
               bary_Seg1+bary_Seg2+bary_Seg3+bary_Seg4]

P += list_plot(bary, plotjoined=True, color = 'violet')


### regarding complex solutions 

fakt_g2 = 40*a[7]^2 - 32*a[7] + 7
J2 = R.ideal(gens+[fakt_g2])
nf2s = [J2.reduce(a[i]) for i in range(basislength)]
