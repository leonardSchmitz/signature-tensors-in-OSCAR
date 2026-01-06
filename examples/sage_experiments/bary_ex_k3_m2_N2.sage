from sage.algebras.lie_algebras.bch import bch_iterator
import numpy as np
import re
import time

channels = 2  # d
depth = 3  # k
sample_size = 2 # N
segments = 2 # m


#### compute signatures of 2 paths with 2 segments each
Lbl = LieAlgebra(QQ, channels, step=depth)
basislength = len(Lbl.gens())
vars_path = channels*segments
vars_all = sample_size*vars_path
R_a = PolynomialRing(QQ, vars_all, 'a')
a = R.gens()
L_a = LieAlgebra(R_a, channels, step=depth)
L_a.inject_variables()

logSigPath1Seg1 = L_a.from_vector(vector(R_a, a[:2]+(0,0,0)))  # 3 = basis_length - channels
logSigPath1Seg2 = L_a.from_vector(vector(R_a, a[2:4]+(0,0,0)))
logSigPath2Seg1 = L_a.from_vector(vector(R_a, a[4:6]+(0,0,0)))
logSigPath2Seg2 = L_a.from_vector(vector(R_a, a[6:8]+(0,0,0)))  # 8 = vars_all, (d,k,N)=(2,2,2)

logSigPath1 = sum([Z for Z in bch_iterator(logSigPath1Seg1, logSigPath1Seg2)])
logSigPath2 = sum([Z for Z in bch_iterator(logSigPath2Seg1, logSigPath2Seg2)])

LynCorPath1 = list(logSigPath1.monomial_coefficients().values()) 
LynCorPath2 = list(logSigPath2.monomial_coefficients().values()) 



########  try a compleatly new approach (this works best for the moment): 
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
res2 = vector([rep2(f) for f in list(res1)])

# insert 'generic' sample 

## uniwue solution segement path [1/4,1/2] -> [0,1/2]   leaves the polytope
path1Seg1 = [0,0]
path1Seg2 = [0,0]
path2Seg1 = [1/2,1/2]
path2Seg2 = [1/2,-1/2]
#J.groebner_basis()
#[a0 - 1/4, a1 - 1/2, a2 - 1/4, a3 + 1/2]

# generically no solution
#path1Seg1 = [1/2,1/44]
#path1Seg2 = [1/33,-1/2]
#path2Seg1 = [1/2,-1/23]
#path2Seg2 = [1/2,3/23]
#dim(J)
#-1

applyBaryOnLyndCoo = R.hom([path1Seg1[0],path1Seg1[1],
                            path1Seg2[0],path1Seg2[1],
                            path2Seg1[0],path2Seg1[1],
                            path2Seg2[0],path2Seg2[1],
                            0,0,0,0,0]) # last number are not evaluated anyhow

gens = [LynCorPath1[i] - applyBaryOnLyndCoo(res2[i]) for i in range(basislength)]
J = R.ideal(gens)


#### compute polynomials for the bary map (mostly old code) 

R = PolynomialRing(QQ, 2*basislength, 'm')
m = R.gens()
L = LieAlgebra(PolynomialRing(QQ, 2*basislength, 'm'), channels, step=depth)
L.inject_variables()
logSXn = L.from_vector(vector(R, m[basislength:]))
y = -L.from_vector(vector(R, m[:basislength]))
S = sum([Z for Z in bch_iterator(y, logSXn)])
Sdic = S.monomial_coefficients()
Sval = list(Sdic.values())                # derive polynomials 
for i in range(len(Sval)):
    Sval[i] += m[i]

def string_preproc(st_poly):
    result = re.sub(r'\^', '**',   st_poly)
    result = re.sub(r'/', '/',   result)
    result = re.sub(r'()m([0-9]*)()', r'\1m[\2]\3', result)
    #result = re.sub(r'/', '//',   st_poly)
    return "lambda m : "+result

outp_0 = [string_preproc(str(a)) for a in Sval]
outp = [eval(a) for a in outp_0]


#### aBCH
Ltemp = LieAlgebra(QQ,['B','C'])
print(Ltemp)
Lyn = Ltemp.Lyndon()
B,C = Lyn.gens()
bch = bch_iterator(B,C)
Stemp = Ltemp();
for i in range(floor((depth + 1)/ 2)) :
  next(bch)
  T3list = next(bch).list()
  Stemp  = Stemp + sum([Lyn({item[0][1]: 2*item[1]}) for item in T3list])

def apply_bracket(bracket, Bsub, Csub, Liealg) :
  if (str(bracket) == 'B'): return Bsub
  if (str(bracket) == 'C'): return Csub
  return Liealg.bracket(apply_bracket(bracket[0], Bsub,Csub, Liealg), apply_bracket(bracket[1], Bsub,Csub, Liealg))

Snew = sum([item[1] * (apply_bracket(item[0],y,logSXn,L)) for item in Stemp.list()])
SvalNew = list(Snew.monomial_coefficients().values())
outp_New_0 = [string_preproc(str(a)) for a in SvalNew]
outp_New = [eval(a) for a in outp_New_0]



#### bring everything together

sample_logSig = matrix([LynCorPath1,LynCorPath2])  # N x lambda , here 2 x 5

res = [0*f for f in LynCorPath1]
#sage: type([0*f for f in LynCorPath1][1])
#<class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
for k in range(basislength):
    res[k] = sum([outp[k](list(res)+[sample_logSig[n,_k] for _k in range(basislength) ]) for n in range(sample_size)])/sample_size

# try again:
res = [0*f for f in LynCorPath1]
for k in range(basislength):
    temp = res[0]*0
    for n in range(sample_size):
        temp = temp + outp[k](list(res)+[sample_logSig[n,_k] for _k in range(basislength) ])
        print(k,n)
        print(list(res)+[sample_logSig[n,_k] for _k in range(basislength) ])
        print(outp[k](list(res)+[sample_logSig[n,_k] for _k in range(basislength) ]))
        print()
    res[k] = temp/sample_size

res_New = [0*f for f in LynCorPath1]
#sage: type([0*f for f in LynCorPath1][1])
#<class 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
for k in range(basislength):
    res_New[k] = sum([outp_New[k](list(res_New)+[sample_logSig[n,_k] for _k in range(basislength) ]) for n in range(sample_size)])/sample_size

# try again:
res_New = [0*f for f in LynCorPath1]
for k in range(basislength):
    temp = res_New[0]*0
    for n in range(sample_size):
        temp = temp + outp_New[k](list(res_New)+[sample_logSig[n,_k] for _k in range(basislength) ])
        print(k,n)
        print(list(res_New)+[sample_logSig[n,_k] for _k in range(basislength) ])
        print(outp_New[k](list(res_New)+[sample_logSig[n,_k] for _k in range(basislength) ]))
        print()
    res_New[k] = temp/sample_size

def string_preproc_var_a(st_poly):
    result = re.sub(r'\^', '**',   st_poly)
    result = re.sub(r'/', '//',   result)
    result = re.sub(r'()a([0-9]*)()', r'\1a[\2]\3', result)
    #result = re.sub(r'/', '//',   st_poly)
    return "lambda a : "+result


