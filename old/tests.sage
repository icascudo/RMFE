#!/usr/bin/env sage

from RMFE_23 import *
from RMFE_24 import *
from RMFE_35 import *

#Tests for RMFE_23, RMFE_24, RMFE_35. The tests generate two random inputs v,w in F_2^k, apply the phi function for the corresponding RMFE, obtaining vectors f,g in a larger field H (each coordinate of f is phi applied to a subvector of v, same with g and w); the test then computes the product h in H (coordinatewise with the field product), applies psi to the result and checks that the result really corresponds to the coordinatewise products of v,w.



################################ TEST FOR RMFE_23 ###########################################

print ('TEST FOR RMFE_23')

#Input length (can be changed)
input_length=31

v=[]
w=[]


#Creating random vectors.
for i in range(input_length): 
	v.append(GF(2).random_element())
	w.append(GF(2).random_element())
        
#Applying phi(first map from the RMFE) to each vector. 
f=phi_RMFE23(v)
g=phi_RMFE23(w)
#print ('Initial vector', v, 'embedded in', f, 'whose elements belong to', (f[0]).parent())
#print ('Initial vector', w, 'embedded in', g, 'whose elements belong to', (g[0]).parent())
h=[f[i]*g[i] for i in range(len(f))]
#print ('Product', h)

#The vector profile will tell psi how to parse the output: this is because, since the implementation of phi allows vectors of arbitrary length and parses the input (splits the input in sub-blocks) accordingly (as a vector in (F_2)^s x...x (F_2)^s x F_2^(s'), s=18, s'<=18 ), psi needs to produce the same parsing, which is input as a vector profile=(s,s,...,s,s')
profile=[]
number_blocks=len(v)//18+1
for j in range(number_blocks-1):
	profile.append(18)
profile.append(len(v)%18)


#Applying psi
m=psi_RMFE23(h,profile)

#Computing the coordinatewise product of original inputs, for testing.
m_sup=[v[i]*w[i] for i in range(len(v))]

#Checking both values are the same.
print ('Is m as supposed (coordinatewise product of initial vectors)?', m==m_sup)


#Timing:
#Timing phi:
#t=timeit('f=phi_RMFE23(v)', seconds=True)
#Timing psi:
#t1=timeit('m=psi_RMFE23(h,profile)', seconds=True)
#print (t, t1)

##################################TEST FOR RMFE_24#############################################################

print ('TEST FOR RMFE_24')

#Input length (can be changed)
input_length=123

v=[]
w=[]

#Creating random vectors.
for i in range(input_length): 
	v.append(GF(2).random_element())
	w.append(GF(2).random_element())

#Applying phi(first map from the RMFE) to each vector. 
f=phi_RMFE24(v)
g=phi_RMFE24(w)
#print ('Initial vector', v, 'embedded in', f, 'whose elements belong to', (f[0]).parent())
#print ('Initial vector', w, 'embedded in', g, 'whose elements belong to', (g[0]).parent())
h=[f[i]*g[i] for i in range(len(f))]
#print ('Product', h)


#The vector profile will tell psi how to parse the output: this is because, since the implementation of phi allows vectors of arbitrary length and parses the input (splits the input in sub-blocks) accordingly (as a vector in (F_2)^s x...x (F_2)^s x F_2^(s'), s=32, s'<=32 ), psi needs to produce the same parsing, which is input as a vector profile=(s,s,...,s,s')
profile=[]
number_blocks=len(v)//32+1
for j in range(number_blocks-1):
	profile.append(32)
profile.append(len(v)%32)

#Applying psi
m=psi_RMFE24(h,profile)

#Computing the coordinatewise product of original inputs, for testing.
m_sup=[v[i]*w[i] for i in range(len(v))]

#Checking both values are the same.
print ('Is m as supposed (coordinatewise product of initial vectors)?', m==m_sup)


#Timing:
#Timing phi:
#t=timeit('f=phi_RMFE24(v)', seconds=True)
#Timing psi:
#t1=timeit('m=psi_RMFE24(h,profile)', seconds=True)
#print (t, t1)

################################ TEST FOR RMFE_35 ###########################################

print ('TEST FOR RMFE_35')
#Input length (can be changed)
input_length=512

v=[]
w=[]

#Creating random vectors.
for i in range(input_length):
	v.append(GF(2).random_element())
	w.append(GF(2).random_element())

#Applying phi(first map from the RMFE) to each vector. 
f=phi_RMFE35(v)
g=phi_RMFE35(w)
#print ('Initial vector', v, 'embedded in', f, 'whose elements belong to', (f[0]).parent())
#print ('Initial vector', w, 'embedded in', g, 'whose elements belong to', (g[0]).parent())
h=[f[i]*g[i] for i in range(len(f))]
#print ('Product', h)

#The vector profile will tell psi how to parse the output: this is because, since the implementation of phi allows vectors of arbitrary length and parses the input (splits the input in sub-blocks) accordingly (as a vector in (F_2)^s x...x (F_2)^s x F_2^(s'), s=99, s'<=99 ), psi needs to produce the same parsing, which is input as a vector profile=(s,s,...,s,s')
profile=[]
number_blocks=len(v)//99+1
for j in range(number_blocks-1):
	profile.append(99)
profile.append(len(v)%99)

#Applying psi
m=psi_RMFE35(h,profile)

#Computing the coordinatewise product of original inputs, for testing.
m_sup=[v[i]*w[i] for i in range(len(v))]

#Checking both values are the same.
print ('Is m as supposed (coordinatewise product of initial vectors)?', m==m_sup)


#Timing:
#Timing phi:
#t=timeit('f=phi_RMFE35(v)', seconds=True)
#Timing psi:
#t1=timeit('m=psi_RMFE35(h,profile)', seconds=True)
#print (t, t1)


