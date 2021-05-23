
# This file was *autogenerated* from the file RMFE_24.sage
from sage.all_cmdline import *   # import sage library

_sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_4 = Integer(4); _sage_const_128 = Integer(128); _sage_const_16 = Integer(16); _sage_const_32 = Integer(32)#!/usr/bin/env sage

from FFTpreproc import *
from FFTa import * 
from field_iso import *

##Implementation of a (32,128)_2 RMFE, as a concatenation of a (2,4)_2 RMFE and a (8,32)_16 RMFE. Functionality is extended to deal with (k,128)_2 RMFEs meaning that: if k<32, phi can be applied to an input (F_2)^k by first padding with zeros to get a vector in (F_2)^32; and if k>32, then the input vector to phi is splitted in subvectors of length 32 (and a remainder block of length k') and the map phi is applied to each vector. In the case of the map psi, it receives a vector of elements in the larger field H=F_(2^128) and applies psi:H->F_2^k to each coordinate, where the k's are specified by the user (they need to be k<=32) and can be different.

##Functions:
#map23: map phi of the (2,4)_2-RMFE: takes as input a vector of (F_2)^2, outputs a single element in the field F=F_16. 
#invmap23: map psi of (2,4)_2-RMFE: takes as input a single element in the field F=F_16, outputs a vector in (F_2)^2.
#phi_RMFE23: map phi of (k,128)_2-RMFE with "extended functionality": takes as input a vector in F_2^k, outputs a vector in H^m, where H=F_{2^128} (field of 2^128 elements) and m=k/32+1.
#psi_RMFE23: map psi of (k,128)_2-RMFE with "extended functionality": takes as input a vector in H^m, where H=F_{2^128} (field of 2^128 elements), and a vector of integers (k_1,...,k_m) with k_i<=32, and outputs a vector in F_2^(k_1+k_2+...+k_m), where each block is the result of applying psi to the i-th component of the input H-vector.  


H = GF(_sage_const_2 **_sage_const_128 , modulus="primitive", names=('c',)); (c,) = H._first_ngens(1)
h=H.modulus()

F = GF(_sage_const_2 **_sage_const_4 , names=('a',)); (a,) = F._first_ngens(1)
f=F.modulus()

P = PolynomialRing(GF(_sage_const_2 ), names=('X',)); (X,) = P._first_ngens(1)
R = PolynomialRing(F, names=('Y',)); (Y,) = R._first_ngens(1)

g=R(h).factor()[_sage_const_0 ][_sage_const_0 ]






#map24: phi-map of (2,4)_2-RMFE

def map24(v): 
    return v[_sage_const_0 ]+(v[_sage_const_0 ]+v[_sage_const_1 ])*a




#invmap24: psi-map of (2,4)_2-RMFE

def invmap24(d): 
    if d!=_sage_const_0 :
	p=d.polynomial(X)
	return [p(_sage_const_0 ),p(_sage_const_1 )]
    else:
	return [_sage_const_0 ,_sage_const_0 ]
    return D




#phi_RMFE24: If given a binary vector of length k<=32, computes phi-map of (k,128)_2-RMFE. Else, split in blocks of length 32 (and a remainder block of length k') and compute phi-map of (32,128)_2-RMFE on each block (and (k',128)_2-RMFE on the last one). Outputs a list of elements in GF(2^128)

def phi_RMFE24(v): 
   
    if len(v)>_sage_const_32 :
	w=[]
        number_blocks=len(v)//_sage_const_32 +_sage_const_1 
	for i in range(number_blocks-_sage_const_1 ):
		w.append(v[_sage_const_32 *i:_sage_const_32 *i+_sage_const_32 ])
	w.append(v[_sage_const_32 *(number_blocks-_sage_const_1 ):])
	res=[]
        for j in range(len(w)):
		res=res+phi_RMFE24(w[j])
	return res

    #First step, split the binary vector in blocks of two cordinates (fill in left-over block with zeros), apply (2,4)_2-RMFE to each block. 
    else: 
	k=len(v)
    	odd=k%_sage_const_2 
    
    	if odd:
        	v.append(_sage_const_0 )
		l=_sage_const_1 
    	else:
		l=_sage_const_0 

    	v1=[]
    	for i in range((k+l)//_sage_const_2 ):
        	t=map24(v[_sage_const_2 *i:_sage_const_2 *i+_sage_const_2 ])
        	v1.append(t)            
    
    	#Second step, apply (k',32)_16-RMFE to result v1. Here k'=length(v1)<=16. Apply inverse FFT to v1, obtaining an interpolating polynomial of degree <=15. Map the <=15-degree polynomial into an element of F_(16^32) and represent it as an element of F_(2^128) via field_iso_desc 
    	while len(v1)<_sage_const_16 :
    		v1.append(_sage_const_0 )
    	

    	#We generate the preprocessing data for the FFT (TODO: having this as a preprocessing would only be useful if the precomputation would be used for more than one evaluation of phi, not the case currently).
    	B=[a**i for i in range(_sage_const_4 )]
    	data=FFTpreproc(_sage_const_4 ,B)

    	#Apply inverse FFT
    	v2=invbinaryFFT(v1,_sage_const_4 ,B,data)
    	
    
    	#Represent result as polynomial
    	m=_sage_const_0 
    	for i in range(len(v2)):
        	m+=v2[i]*Y**i

     
    	#Map the <=15-degree polynomial from F_16[X] into an element of F_(2^128) via field_iso_desc  (by implicitely first mapping into an element of F_(16^32) and then changing to a representation in F_2^(128)).    
    	r=field_iso_desc(m,_sage_const_4 ,g,h,F,H,P,R)
    
    	
    return [r]




#psi_RMFE24: Given a list of elements w of F_(2^128), and a vector of values k<=32, computes psi-map of (k,128)_2-RMFE on each element and outputs the concatenation of the resulting vectors

def psi_RMFE24(w,k):
   
    if len(w)!=len(k):
	raise Exception("inputs to psi_RMFE24 must be of same length")
    for i in range(len(k)): 
	if k[i]>_sage_const_32 :
		raise Exception("every coordinate on second input of psi_RMFE24 needs to be at most 32")
    B=[a**i for i in range(_sage_const_4 )]
    data=FFTpreproc(_sage_const_4 ,B)
    
    res=[]
    for j in range(len(w)):
    
    	#First change field representation to represent input as element of F_(32^65) and hence as a polynomial in F_32[X] of degree at most 64.
    	m=field_iso_asc(w[j],_sage_const_4 ,g,R)
    	m=list(m)
    
    	#Before applying the FFT we need to a polynomial of degree <=15. For this we take modulo X^16+X, as this does not modify evaluation in points of F_16:
    	hred=listsum(m[_sage_const_0 :_sage_const_16 ],[_sage_const_0 ]+m[_sage_const_16 :])
        
    	
        
    	#Apply FFT
    	w1=binaryFFT(hred,_sage_const_4 ,B,data)
    
    	
    
    	#Based on value of k, we adjust size of the output. 
    	upper=(k[j]+_sage_const_1 )//_sage_const_2 	
    	del w1[upper:]



    	#Apply psi from (2,4)_2-RMFE to each element of resulting vector.
    	r=[]
    	for i in range(len(w1)):
    		r=r+invmap24(w1[i])
    
    	#Adjust size of output.
    	del r[k[j]:]

	#Concatenate this to global vector.
        res=res+r
    return res
    

