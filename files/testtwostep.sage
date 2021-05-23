#!/usr/bin/env sage

from twostepRMFE import *
from twostepinstance import *

################################ TEST FOR 2step RMFE ###########################################

################################ PARAMETERS #################################### 

#Instance parameters (k1,e1,k2,e2) where k1=2 or 3, 2*k1<=e1+1, k2<=2^e1+1, 2*k2-1<=e2 
instance=twostepinstance(3,5,33,65)

k1=instance.k1
k2=instance.k2
e1=instance.e1
e2=instance.e2
k=instance.k
e=instance.e
m1=instance.m1
m2=instance.m2
F=instance.F
a=instance.a
H=instance.H
c=instance.c
P=instance.P
X=instance.X
R=instance.R
Y=instance.Y
f=instance.f
h=instance.h
g=instance.g

input_length=5*k-4
################################ START TEST #################################### 
print ('TEST FOR 2step RME')

#Creating random vectors
v=[]
w=[]

for i in range(input_length): 
	v.append(GF(2).random_element())
	w.append(GF(2).random_element())
        
#Applying phi(first map from the RMFE) to each vector. 
f=phi_RMFE(v,instance)
g=phi_RMFE(w,instance)
print ('Initial vector', v, 'is embedded in', f, 'whose elements belong to', (f[0]).parent())
print ('Initial vector', w, 'is embedded in', g, 'whose elements belong to', (g[0]).parent())
h=[f[i]*g[i] for i in range(len(f))]
print ('Product', h)

#The vector profile will tell psi how to parse the output: this is because, since the implementation of phi allows vectors of arbitrary length and parses the input (splits the input in sub-blocks) accordingly (as a vector in (F_2)^k x...x (F_2)^k x F_2^(k'), k'<=k ), psi needs to produce the same parsing, which is input as a vector profile=(s,s,...,s,s')
profile=[]
number_blocks=((len(v)-1)//k)+1
for j in range(number_blocks-1):
	profile.append(k)
profile.append(len(v)-(number_blocks-1)*k)

#Applying psi
m=psi_RMFE(h,profile,instance)

#Computing the coordinatewise product of original inputs, for testing.
m_sup=[v[i]*w[i] for i in range(len(v))]

#Checking both values are the same.
print ('Is m as supposed (coordinatewise product of initial vectors)?', m==m_sup)
print ('Computed', m)
print ('Supposed', m_sup)

#Timing:
#Timing phi:
#t=timeit('f=phi_RMFE(v,instance)', seconds=True)
#Timing psi:
#t1=timeit('m=psi_RMFE(h,profile,instance)', seconds=True)
#print (t, t1)


#########################################



