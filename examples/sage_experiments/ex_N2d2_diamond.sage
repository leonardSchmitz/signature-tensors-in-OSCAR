from sage.algebras.lie_algebras.bch import bch_iterator
import numpy as np
import re
import time
import random

def random_color():
    return "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])

channels = 2  # d
depth = 3  # k
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

logSigPath1Seg1 = L.from_vector(vector(R, a[:2]+(0,0,0)))  # 3 = basis_length - channels
logSigPath1Seg2 = L.from_vector(vector(R, a[2:4]+(0,0,0)))
logSigPath2Seg1 = L.from_vector(vector(R, a[4:6]+(0,0,0)))
logSigPath2Seg2 = L.from_vector(vector(R, a[6:8]+(0,0,0)))  # 8 = vars_all, (d,k,N)=(2,2,2)

logSigPath1 = sum([Z for Z in bch_iterator(logSigPath1Seg1, logSigPath1Seg2)])
logSigPath2 = sum([Z for Z in bch_iterator(logSigPath2Seg1, logSigPath2Seg2)])

LynCorPath1 = list(logSigPath1.monomial_coefficients().values())
LynCorPath2 = list(logSigPath2.monomial_coefficients().values())

y = -L.from_vector(vector(R, a[8:basislength+8]))
S = sum([Z for Z in bch_iterator(y, logSigPath1)]) + sum([Z for Z in bch_iterator(y, logSigPath2)])

coef_S = vector(list(S.monomial_coefficients().values()))
res0 = 1/2*(coef_S + 2*vector(a[8:8+basislength]))
rep1 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res0[0],res0[1],a[10],a[11],a[12]])
res1 = vector([rep1(f) for f in list(res0)])
rep2 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res0[0],res0[1],res1[2],a[11],a[12]])
res = vector([rep2(f) for f in list(res1)])

## ========== explicit sample "diamond" 

origin = vector([0,0])
path1Seg1 = vector([1/2,1/2])
path1Seg2 = vector([1/2,-1/2])
path2Seg1 = vector([1/2,-1/2])
path2Seg2 = vector([1/2,1/2])

path1 = [origin,path1Seg1,path1Seg1+path1Seg2]
path2 = [origin,path2Seg1,path2Seg1+path2Seg2]

#P = list_plot(path1, plotjoined=True, color='blue')
#P += list_plot(path2, plotjoined=True, color='red')

## evaluation homomorphism
applyBaryOnLyndCoo = R.hom([path1Seg1[0],path1Seg1[1],
                            path1Seg2[0],path1Seg2[1],
                            path2Seg1[0],path2Seg1[1],
                            path2Seg2[0],path2Seg2[1],
                            0,0,0,0,0]) # last number are not evaluated anyhow

evalLynCorPath1 =  [applyBaryOnLyndCoo(f) for f in LynCorPath1] 
evalLynCorPath2 =  [applyBaryOnLyndCoo(f) for f in LynCorPath2]

evalLynCorBary =  [applyBaryOnLyndCoo(f) for f in res]


## ========== recovery of the barycenter
## two-segement path has no solution:

#gens = [LynCorPath1[i] - applyBaryOnLyndCoo(res2[i]) for i in range(basislength)]
#J = R.ideal(gens) 
# ---->  dim(J)=-1

## three segments:
logSigBarySeg1 = L.from_vector(vector(R, a[:2]+(0,0,0)))  # 3 = basis_length - channels
logSigBarySeg2 = L.from_vector(vector(R, a[2:4]+(0,0,0)))
logSigBarySeg3 = L.from_vector(vector(R, a[4:6]+(0,0,0)))

logSigBarySeg12 = sum([Z for Z in bch_iterator(logSigBarySeg1, logSigBarySeg2)])
logSigBary = sum([Z for Z in bch_iterator(logSigBarySeg12, logSigBarySeg3)])

LynCorBary = list(logSigBary.monomial_coefficients().values())
gens = [LynCorBary[i] - applyBaryOnLyndCoo(res[i]) for i in range(basislength)]
J = R.ideal(gens)   
#dim(J)-13+6 ---> dim(J) = 1
#J.groebner_basis()  ----> [a3*a5 + a5^2 + 1/4, a0 - 1, a1 + a3 + a5, a2 + 1, a4 - 1]


## a[1] = 1/2
J_a1_1 = R.ideal(gens+[a[1]-1/2])

bary_a1_1_Seg1 = vector([J_a1_1.reduce(a[0]),J_a1_1.reduce(a[1])])
bary_a1_1_Seg2 = vector([J_a1_1.reduce(a[2]),J_a1_1.reduce(a[3])])
bary_a1_1_Seg3 = vector([J_a1_1.reduce(a[4]),J_a1_1.reduce(a[5])])

bary_a1_1 = [origin,bary_a1_1_Seg1, 
                    bary_a1_1_Seg1+bary_a1_1_Seg2, 
                    bary_a1_1_Seg1+bary_a1_1_Seg2+bary_a1_1_Seg3]


P += list_plot(bary_a1_1, plotjoined=True, color='purple')



## a[1] = 1/3
J_a1_1 = R.ideal(gens+[a[1]-1/3])

bary_a1_1_Seg1 = vector([J_a1_1.reduce(a[0]),J_a1_1.reduce(a[1])])
bary_a1_1_Seg2 = vector([J_a1_1.reduce(a[2]),J_a1_1.reduce(a[3])])
bary_a1_1_Seg3 = vector([J_a1_1.reduce(a[4]),J_a1_1.reduce(a[5])])

bary_a1_1 = [origin,bary_a1_1_Seg1, 
                    bary_a1_1_Seg1+bary_a1_1_Seg2, 
                    bary_a1_1_Seg1+bary_a1_1_Seg2+bary_a1_1_Seg3]


P += list_plot(bary_a1_1, plotjoined=True, color='pink')


## four segments:

logSigBarySeg1 = L.from_vector(vector(R, a[:2]+(0,0,0)))  # 3 = basis_length - channels
logSigBarySeg2 = L.from_vector(vector(R, a[2:4]+(0,0,0)))
logSigBarySeg3 = L.from_vector(vector(R, a[4:6]+(0,0,0)))
logSigBarySeg4 = L.from_vector(vector(R, a[6:8]+(0,0,0)))

logSigBarySeg12 = sum([Z for Z in bch_iterator(logSigBarySeg1, logSigBarySeg2)])
logSigBarySeg123 = sum([Z for Z in bch_iterator(logSigBarySeg12, logSigBarySeg3)])
logSigBary = sum([Z for Z in bch_iterator(logSigBarySeg123, logSigBarySeg4)])

LynCorBary = list(logSigBary.monomial_coefficients().values())
gens = [LynCorBary[i] - applyBaryOnLyndCoo(res[i]) for i in range(basislength)]
J = R.ideal(gens)   


#picture 1:
#J1 = R.ideal(gens + [a[0] - 1/3, a[1] - 1/3, a[2] -1/6])
#J2 = R.ideal(gens + [a[0] - 1/3, a[1] - 1/3, a[2]])
#J3 = R.ideal(gens + [a[0] - 1/3, a[1] - 1/3, a[2] + 1/6])
#J7 = R.ideal(gens + [a[0] - 1/2, a[1] - 1/2, a[2]])
#J8 = R.ideal(gens + [a[0] - 2/5, a[1] - 2/5, a[2]])
#J9 = R.ideal(gens + [a[0] - 2/5, a[1] - 3/7, a[2]])
#J4 = R.ideal(gens + [a[0] - 1/4, a[1] - 1/4, a[2] -1/4])
#J5 = R.ideal(gens + [a[0] - 1/4, a[1] - 1/4, a[2] -2/3])
#J6 = R.ideal(gens + [a[0] - 1/4, a[1] - 1/4, a[2] -2/5])
#J10 = R.ideal(gens + [a[0] - 1/5, a[1] - 1/5, a[2] -5/9])


#picture 2:
#J1 = R.ideal(gens + [a[0] - 1/5, a[1] - 1/6, a[2] -1/2])
#J2 = R.ideal(gens + [a[0] - 1/5, a[1] - 1/6, a[2] -3/5])
#J3 = R.ideal(gens + [a[0] - 1/5, a[1] - 1/6, a[2] -4/7])
#J4 = R.ideal(gens + [a[0] - 1/5, a[1] - 1/6, a[2] -5/8])
#J10 = R.ideal(gens + [a[0] - 1/5, a[1] - 1/5, a[2] -5/9])
#J5 = R.ideal(gens + [a[0] - 1/4, a[1] - 1/4, a[2] -1/4])
#J6 = R.ideal(gens + [a[0] - 1/4, a[1] - 1/4, a[2] -1/3])
#J7 = R.ideal(gens + [a[0] - 1/4, a[1] - 1/4, a[2] -2/5])
#J8 = R.ideal(gens + [a[0] - 1/4, a[1] - 1/5, a[2] -1/2])
#J9 = R.ideal(gens + [a[0] - 1/4, a[1] - 1/5, a[2] -3/7])


#picture 3:
J1 = R.ideal(gens + [a[0] - 2/3, a[1] - 1/3, a[2] +1/2])
J4 = R.ideal(gens + [a[0] - 2/3, a[1] - 1/3, a[2] +1/3])
J2 = R.ideal(gens + [a[0] - 4/5, a[1] - 1/5, a[2] +1/2])
J3 = R.ideal(gens + [a[0] - 4/5, a[1] - 1/5, a[2] +2/3])
J5 = R.ideal(gens + [a[0] - 1/2, a[1] - 1/2, a[2] +1/3])
J6 = R.ideal(gens + [a[0] - 1/2, a[1] - 1/2, a[2] +1/4])
J7 = R.ideal(gens + [a[0] - 1/2, a[1] - 1/2, a[2] +2/5])
J8 = R.ideal(gens + [a[0] - 2/3, a[1] - 1/4, a[2] +1/3])
J9 = R.ideal(gens + [a[0] - 2/3, a[1] - 1/3, a[2] +1/3])
J10 = R.ideal(gens + [a[0] - 1/2, a[1] - 1/3, a[2] +1/3])

Js = [J1,J2,J3,J4,J5,J6,J7,J8,J9,J10]

P = list_plot(path1, plotjoined=True, color='blue')
P += list_plot(path2, plotjoined=True, color='red')

for J in Js:
   bary_Seg1 = vector([J.reduce(a[0]),J.reduce(a[1])])
   bary_Seg2 = vector([J.reduce(a[2]),J.reduce(a[3])])
   bary_Seg3 = vector([J.reduce(a[4]),J.reduce(a[5])])
   bary_Seg4 = vector([J.reduce(a[6]),J.reduce(a[7])])
   bary = [origin,bary_Seg1,
                  bary_Seg1+bary_Seg2,
                  bary_Seg1+bary_Seg2+bary_Seg3,
                  bary_Seg1+bary_Seg2+bary_Seg3+bary_Seg4]
   P += list_plot(bary, plotjoined=True, color = random_color())


## five segments:

logSigBarySeg1 = L.from_vector(vector(R, a[:2]+(0,0,0)))  # 3 = basis_length - channels
logSigBarySeg2 = L.from_vector(vector(R, a[2:4]+(0,0,0)))
logSigBarySeg3 = L.from_vector(vector(R, a[4:6]+(0,0,0)))
logSigBarySeg4 = L.from_vector(vector(R, a[6:8]+(0,0,0)))
logSigBarySeg5 = L.from_vector(vector(R, a[8:10]+(0,0,0)))

logSigBarySeg12 = sum([Z for Z in bch_iterator(logSigBarySeg1, logSigBarySeg2)])
logSigBarySeg123 = sum([Z for Z in bch_iterator(logSigBarySeg12, logSigBarySeg3)])
logSigBarySeg1234 = sum([Z for Z in bch_iterator(logSigBarySeg123, logSigBarySeg4)])
logSigBary = sum([Z for Z in bch_iterator(logSigBarySeg1234, logSigBarySeg5)])

LynCorBary = list(logSigBary.monomial_coefficients().values())
gens = [LynCorBary[i] - applyBaryOnLyndCoo(res[i]) for i in range(basislength)]
J = R.ideal(gens)

# picture 1:
J1 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] , a[3] +2/3, a[4] - 1/3])
J2 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] , a[3] +2/3, a[4] - 1/2])
J3 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] , a[3] +2/3, a[4] - 2/5])
J4 = R.ideal(gens + [a[0] - 1/4, a[1] -1/4, a[2] , a[3] +1/2, a[4] - 1/2])
J8 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] +1/10 , a[3] , a[4] - 4/6])
J7 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] +1/10 , a[3] , a[4] - 4/6])
J5 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2]  , a[3]+1/2,  a[4] - 1/5])
J6 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2]  , a[3]+2/3, a[4] - 2/6])
J10 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2], a[3] + 2/3, a[4] - 1/2])
J9 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2]  , a[3] +1, a[4] - 2/6])

# picture 2  (understanding J9 a bit better)
#J1 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2]  , a[3] +1, a[4] - 2/6])
#J2 = R.ideal(gens + [a[0] - 2/5, a[1] -2/5, a[2]  , a[3] +4/5, a[4] + 2/10])
#J3 = R.ideal(gens + [a[0] - 2/5, a[1] -2/5, a[2]  , a[3] +4/5, a[4] + 2/10])
#J4 = R.ideal(gens + [a[0] - 2/5, a[1] -2/5, a[2]  , a[3] +4/5, a[4] - 1/10])
#J5 = R.ideal(gens + [a[0] - 2/5, a[1] -2/5, a[2]  , a[3] +4/5, a[4] - 2/10])
#J6 = R.ideal(gens + [a[0] - 2/5, a[1] -2/5, a[2]  , a[3] +4/5, a[4] - 3/10])
#J7 = R.ideal(gens + [a[0] - 2/5, a[1] -2/5, a[2]  , a[3] +4/5, a[4] - 4/10])
#J8 = R.ideal(gens + [a[0] - 2/5, a[1] -2/5, a[2]  , a[3] +4/5, a[4] - 5/10])
#J9 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2]  , a[3] +1, a[4] - 2/6])
#J10 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2]  , a[3] +1, a[4] - 2/6])

# picture 3 (understadning J1)
J1 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] , a[3] +2/3, a[4] - 1/3])
J2 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] , a[3] +2/3, a[4] - 1/2])
J3 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] , a[3] +2/3, a[4] - 2/5])
J4 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] + 1/9 , a[3] +2/3, a[4] - 1/3])
J5 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] + 1/9, a[3] +2/3, a[4] - 1/2])
J6 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] + 1/9, a[3] +2/3, a[4] - 2/5])
J7 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] - 1/9, a[3] +2/3, a[4] - 1/3])
J8 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] - 1/9, a[3] +2/3, a[4] - 1/2])
J9 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] - 1/9, a[3] +2/3, a[4] - 2/5])
J10 = R.ideal(gens + [a[0] - 1/3, a[1] -1/3, a[2] , a[3] +2/3, a[4] - 1/3])

Js = [J1,J2,J3,J4,J8,J5,J6,J9,J10,J7]

P = list_plot(path1, plotjoined=True, color='blue')
P += list_plot(path2, plotjoined=True, color='red')

for J in Js:
   bary_Seg1 = vector([J.reduce(a[0]),J.reduce(a[1])])
   bary_Seg2 = vector([J.reduce(a[2]),J.reduce(a[3])])
   bary_Seg3 = vector([J.reduce(a[4]),J.reduce(a[5])])
   bary_Seg4 = vector([J.reduce(a[6]),J.reduce(a[7])])
   bary_Seg5 = vector([J.reduce(a[8]),J.reduce(a[9])])
   bary = [origin,bary_Seg1,
                  bary_Seg1+bary_Seg2,
                  bary_Seg1+bary_Seg2+bary_Seg3,
                  bary_Seg1+bary_Seg2+bary_Seg3+bary_Seg4,
                  bary_Seg1+bary_Seg2+bary_Seg3+bary_Seg4+bary_Seg5]
   P += list_plot(bary, plotjoined=True, color = random_color())

