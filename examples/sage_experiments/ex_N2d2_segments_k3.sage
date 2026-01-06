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

logSigPath1 = L.from_vector(vector(R, a[:2]+(0,0,0)))  # 5 = basis_length - channels
logSigPath2 = L.from_vector(vector(R, a[2:4]+(0,0,0)))

LynCorPath1 = list(logSigPath1.monomial_coefficients().values())
LynCorPath2 = list(logSigPath2.monomial_coefficients().values())

## ========== explicit sample "half_diamond"

origin = vector([0,0])
path1Seg1 = vector([1/2,1/2])
path2Seg1 = vector([1/2,-1/2])

path1 = [origin,path1Seg1]
path2 = [origin,path2Seg1]


applyBaryOnLyndCoo = R.hom([path1Seg1[0],path1Seg1[1],
                            path2Seg1[0],path2Seg1[1],
                            0,0,
                            0,0,
                            0,0,0,0,0]) # last number are not evaluated anyhow


evalLynCorPath1 =  [applyBaryOnLyndCoo(f) for f in LynCorPath1]
evalLynCorPath2 =  [applyBaryOnLyndCoo(f) for f in LynCorPath2]

w1 = 1/2
y = -L.from_vector(vector(R, a[8:basislength+8]))
S = w1*sum([Z for Z in bch_iterator(y, logSigPath1)]) + (1-w1)*sum([Z for Z in bch_iterator(y, logSigPath2)])

 coef_S = vector(list(S.monomial_coefficients().values()))
res0 = coef_S + vector(a[8:8+basislength])
rep1 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res0[0],res0[1],a[10],a[11],a[12]])
res1 = vector([rep1(f) for f in list(res0)])
rep2 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res0[0],res0[1],res1[2],a[11],a[12]])
res = vector([rep2(f) for f in list(res1)])

evalLynCorBary =  [applyBaryOnLyndCoo(f) for f in res]
logSigBarySeg1 = L.from_vector(vector(R, a[:2]+(0,0,0)))  # 3 = basis_length - channels
logSigBarySeg2 = L.from_vector(vector(R, a[2:4]+(0,0,0)))
logSigBarySeg3 = L.from_vector(vector(R, a[4:6]+(0,0,0)))
logSigBarySeg12 = sum([Z for Z in bch_iterator(logSigBarySeg1, logSigBarySeg2)])
logSigBary = sum([Z for Z in bch_iterator(logSigBarySeg12, logSigBarySeg3)])
LynCorBary = list(logSigBary.monomial_coefficients().values())
gens = [LynCorBary[i] - applyBaryOnLyndCoo(res[i]) for i in range(basislength)]
J = R.ideal(gens)

P = list_plot(path1, plotjoined=True, color='blue')
P += list_plot(path2, plotjoined=True, color='red')
for a1 in [1/i for i in range(1,10)]+[-1/i for i in range(1,10)]:
  J1 = R.ideal(gens + [a[1]-a1])
  nf1s = [J1.reduce(a[i]) for i in range(6)]
  
  bary_Seg1 = vector([nf1s[0],nf1s[1]])
  bary_Seg2 = vector([nf1s[2],nf1s[3]])
  bary_Seg3 = vector([nf1s[4],nf1s[5]])
  bary = [origin,bary_Seg1,
                 bary_Seg1+bary_Seg2,
                 bary_Seg1+bary_Seg2+bary_Seg3]
  P += list_plot(bary, plotjoined=True, color = random_color())




## ========== explicit sample with first coordinate 1/2 

origin = vector([0,0])
path1Seg1 = vector([1/2,1/2])
path2Seg1 = vector([1/2,0])

path1 = [origin,path1Seg1]
path2 = [origin,path2Seg1]


applyBaryOnLyndCoo = R.hom([path1Seg1[0],path1Seg1[1],
                            path2Seg1[0],path2Seg1[1],
                            0,0,
                            0,0,
                            0,0,0,0,0]) # last number are not evaluated anyhow


evalLynCorPath1 =  [applyBaryOnLyndCoo(f) for f in LynCorPath1]
evalLynCorPath2 =  [applyBaryOnLyndCoo(f) for f in LynCorPath2]

w1 = 1/2
y = -L.from_vector(vector(R, a[8:basislength+8]))
S = w1*sum([Z for Z in bch_iterator(y, logSigPath1)]) + (1-w1)*sum([Z for Z in bch_iterator(y, logSigPath2)])

coef_S = vector(list(S.monomial_coefficients().values()))
res0 = coef_S + vector(a[8:8+basislength])
rep1 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res0[0],res0[1],a[10],a[11],a[12]])
res1 = vector([rep1(f) for f in list(res0)])
rep2 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res0[0],res0[1],res1[2],a[11],a[12]])
res = vector([rep2(f) for f in list(res1)])

evalLynCorBary =  [applyBaryOnLyndCoo(f) for f in res]
logSigBarySeg1 = L.from_vector(vector(R, a[:2]+(0,0,0)))  # 3 = basis_length - channels
logSigBarySeg2 = L.from_vector(vector(R, a[2:4]+(0,0,0)))
logSigBarySeg3 = L.from_vector(vector(R, a[4:6]+(0,0,0)))
logSigBarySeg12 = sum([Z for Z in bch_iterator(logSigBarySeg1, logSigBarySeg2)])
logSigBary = sum([Z for Z in bch_iterator(logSigBarySeg12, logSigBarySeg3)])
LynCorBary = list(logSigBary.monomial_coefficients().values())
gens = [LynCorBary[i] - applyBaryOnLyndCoo(res[i]) for i in range(basislength)]
J = R.ideal(gens)

P = list_plot(path1, plotjoined=True, color='blue')
P += list_plot(path2, plotjoined=True, color='red')
for a1 in [1/i+1/4 for i in range(1,10)]+[-1/i+1/4 for i in range(1,10)]:
  J1 = R.ideal(gens + [a[1]-a1])
  nf1s = [J1.reduce(a[i]) for i in range(6)]

  bary_Seg1 = vector([nf1s[0],nf1s[1]])
  bary_Seg2 = vector([nf1s[2],nf1s[3]])
  bary_Seg3 = vector([nf1s[4],nf1s[5]])
  bary = [origin,bary_Seg1,
                 bary_Seg1+bary_Seg2,
                 bary_Seg1+bary_Seg2+bary_Seg3]
  P += list_plot(bary, plotjoined=True, color = random_color())





## ========== generic sample

origin = vector([0,0])
path1Seg1 = vector([1/2,1/3])
path2Seg1 = vector([1/4,1/5])

path1 = [origin,path1Seg1]
path2 = [origin,path2Seg1]


applyBaryOnLyndCoo = R.hom([path1Seg1[0],path1Seg1[1],
                            path2Seg1[0],path2Seg1[1],
                            0,0,
                            0,0,
                            0,0,0,0,0]) # last number are not evaluated anyhow


evalLynCorPath1 =  [applyBaryOnLyndCoo(f) for f in LynCorPath1]
evalLynCorPath2 =  [applyBaryOnLyndCoo(f) for f in LynCorPath2]

w1 = 1/2
y = -L.from_vector(vector(R, a[8:basislength+8]))
S = w1*sum([Z for Z in bch_iterator(y, logSigPath1)]) + (1-w1)*sum([Z for Z in bch_iterator(y, logSigPath2)])

coef_S = vector(list(S.monomial_coefficients().values()))
res0 = coef_S + vector(a[8:8+basislength])
rep1 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res0[0],res0[1],a[10],a[11],a[12]])
res1 = vector([rep1(f) for f in list(res0)])
rep2 = R.hom([a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],res0[0],res0[1],res1[2],a[11],a[12]])
res = vector([rep2(f) for f in list(res1)])

evalLynCorBary =  [applyBaryOnLyndCoo(f) for f in res]
logSigBarySeg1 = L.from_vector(vector(R, a[:2]+(0,0,0)))  # 3 = basis_length - channels
logSigBarySeg2 = L.from_vector(vector(R, a[2:4]+(0,0,0)))
logSigBarySeg3 = L.from_vector(vector(R, a[4:6]+(0,0,0)))
logSigBarySeg12 = sum([Z for Z in bch_iterator(logSigBarySeg1, logSigBarySeg2)])
logSigBary = sum([Z for Z in bch_iterator(logSigBarySeg12, logSigBarySeg3)])
LynCorBary = list(logSigBary.monomial_coefficients().values())
gens = [LynCorBary[i] - applyBaryOnLyndCoo(res[i]) for i in range(basislength)]
J = R.ideal(gens)

P = list_plot(path1, plotjoined=True, color='blue')
P += list_plot(path2, plotjoined=True, color='red')
for a1 in [1/i for i in range(5,100)]:
  J1 = R.ideal(gens + [a[1]-a1])
  nf1s = [J1.reduce(a[i]) for i in range(6)]

  bary_Seg1 = vector([nf1s[0],nf1s[1]])
  bary_Seg2 = vector([nf1s[2],nf1s[3]])
  bary_Seg3 = vector([nf1s[4],nf1s[5]])
  bary = [origin,bary_Seg1,
                 bary_Seg1+bary_Seg2,
                 bary_Seg1+bary_Seg2+bary_Seg3]
  P += list_plot(bary, plotjoined=True, color = random_color())
