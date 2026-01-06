from sage.algebras.lie_algebras.bch import bch_iterator
import numpy as np
import re
import time
import random

def random_color():
    return "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])

d = 2  # d
trunc_level = 3  # k
sample_size = 2 # N
m_axis = 3
m_learn = 5

# compute log_sig axis path in lyndon coordinates

L = LieAlgebra(QQ, m_axis, step=trunc_level)
basislength = len(L.gens())
L.inject_variables()
LynBas = L.gens()
logSigAxis = LynBas[0]
for i in range(1,m_axis):  
  logSigAxis = sum([Z for Z in bch_iterator(logSigAxis, LynBas[i])])
logBarySigAxisOrigin = 1/2*logSigAxis
LynCorAxis = list(logBarySigAxisOrigin.monomial_coefficients().values())


# path learning

R = PolynomialRing(QQ, m_learn*m_axis, 'a')
a = R.gens()
L_learn = LieAlgebra(R, m_axis, step=trunc_level)
L_learn.inject_variables()
LynBas_learn = L_learn.gens()
logSig_pwln_learn = sum(a[i]*LynBas_learn[i] for i in range(0,m_axis))
for j in range(1,m_learn):   
  print(j)
  temp = sum(a[i+j*m_axis]*LynBas_learn[i] for i in range(0,m_axis))
  logSig_pwln_learn = sum([Z for Z in bch_iterator(logSig_pwln_learn, temp)])
LynCor_logSig_pwln_learn = list(logSig_pwln_learn.monomial_coefficients().values())

gens = [LynCor_logSig_pwln_learn[i] - LynCorAxis[i] for i in range(basislength)]
J = R.ideal(gens)

for i in range(0,m_learn*m_axis):
  a_no_ai = list(a)
  a_no_ai.remove(a[i])
  print(Jp.elimination_ideal(tuple(a_no_ai)))

#Ideal (a0 - 1) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (a1) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (2*a2 + 1) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (1024*a3^3 + 5184*a3^2 + 4416*a3 + 1039) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (1024*a4^3 + 5056*a4^2 - 3712*a4 + 633) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (128*a5^3 + 176*a5^2 - 160*a5 + 29) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (2048*a6^3 - 5824*a6^2 - 9488*a6 + 927) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (2048*a7^3 - 6976*a7^2 - 25040*a7 + 9441) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (256*a8^3 - 528*a8^2 - 104*a8 + 269) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (2048*a9^3 + 2752*a9^2 - 6512*a9 - 2199) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (2048*a10^3 + 1856*a10^2 - 26160*a10 + 15399) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (256*a11^3 + 272*a11^2 - 664*a11 + 219) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (1024*a12^3 - 2112*a12^2 - 672*a12 - 47) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (1024*a13^3 - 4032*a13^2 - 608*a13 - 9) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field
#Ideal (128*a14^3 - 432*a14^2 + 368*a14 - 93) of Multivariate Polynomial Ring in a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14 over Rational Field

