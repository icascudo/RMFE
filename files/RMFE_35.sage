#!/usr/bin/env sage

from FFTpreproc import *
from FFTa import * 
from field_iso import *

##Implementation of a (99,325)_2 RMFE, as a concatenation of a (3,5)_2 RMFE and a (33,65)_32 RMFE. Functionality is extended to deal with (k,325)_2 RMFEs meaning that: if k<99, phi can be applied to an input (F_2)^k by first padding with zeros to get a vector in (F_2)^99; and if k>32, then the input vector to phi is splitted in subvectors of length 99 (and a remainder block of length k') and the map phi is applied to each vector. In the case of the map psi, it receives a vector of elements in the larger field H=F_(2^325) and applies psi:H->F_2^k to each coordinate, where the k's are specified by the user (they need to be k<=99) and can be different.

##Functions:
#map23: map phi of the (3,5)_2-RMFE: takes as input a vector of (F_2)^3, outputs a single element in the field F=F_32. 
#invmap23: map psi of (3,5)_2-RMFE: takes as input a single element in the field F=F_32, outputs a vector in (F_2)^3.
#phi_RMFE23: map phi of (k,325)_2-RMFE with "extended functionality": takes as input a vector in F_2^k, outputs a vector in H^m, where H=F_{2^325} (field of 2^325 elements) and m=k/99+1.
#psi_RMFE23: map psi of (k,325)_2-RMFE with "extended functionality": takes as input a vector in H^m, where H=F_{2^325} (field of 2^325 elements), and a vector of integers (k_1,...,k_m) with k_i<=99, and outputs a vector in F_2^(k_1+k_2+...+k_m), where each block is the result of applying psi to the i-th component of the input H-vector.  



H.<c>=GF(2^325, modulus="primitive")
h=H.modulus()

F.<a>=GF(2^5)
f=F.modulus()

P.<X>=PolynomialRing(GF(2))
R.<Y>=PolynomialRing(F)

g=R(h).factor()[0][0]






#map35: phi-map of (3,5)_2-RMFE

def map35(v): 
    return v[0]+(v[0]+v[1]+v[2])*a+v[2]*a**2




#invmap35: psi-map of (3,5)_2-RMFE

def invmap35(d): 
    if d!=0:
	p=d.polynomial(X)
	return [p(0),p(1),p[4]]
    else:
	return [0,0,0]
    return D




#phi_RMFE35: If given a binary vector of length k<=99, computes phi-map of (k,325)_2-RMFE. Else, split in blocks of length 99 (and a remainder block of length k') and compute phi-map of (99,325)_2-RMFE on each block (and (k',325)_2-RMFE on the last one). Outputs a list of elements in GF(2^325)

def phi_RMFE35(v): 
    ###print 'Applying phi to', v
    if len(v)>99:
	w=[]
        number_blocks=len(v)//99+1
	for i in range(number_blocks-1):
		w.append(v[99*i:99*i+99])
	w.append(v[99*(number_blocks-1):])
	res=[]
        for j in range(len(w)):
		res=res+phi_RMFE35(w[j])
	return res

    #First step, split the binary vector in blocks of three cordinates (fill in left-over block with zeros), apply (3,5)_2-RMFE to each block. 
    else:
    	k=len(v)
    	modthree=k%3
    
    	if modthree==1:
        	v=v+[0,0]
		l=2
    	elif modthree==2:
		v.append(0)
		l=1
    	else:
		l=0

    	v1=[]
    	for i in range((k+l)//3):
        	t=map35(v[3*i:3*i+3])
        	v1.append(t)            
    
    	#Second step, apply (k',65)_32-RMFE to result v1. Here k'=length(v1)<=33. Apply inverse FFT to first 32 elements, obtaining an interpolating polynomial of degree <=31. If input has a 33-th coordinate, adjust result by     summing an appropriate multiple of X^32-X. Map the <=32-degree polynomial into an element of F_(32^65) and represent it as an element of F_(2^325) via field_iso_desc 
    	while len(v1)<32:
        	v1.append(0)
    	

    	#We generate the preprocessing data for the FFT (TODO: having this as a preprocessing would only be useful if the precomputation would be used for more than one evaluation of phi, not the case currently).
    	B=[a^i for i in range(5)]
    	data=FFTpreproc(5,B)

    	#Apply inverse FFT
    	v2=invbinaryFFT(v1,5,B,data)
    	
    
    	#Represent result as polynomial
    	m=0
    	for i in range(len(v2)):
        	m+=v2[i]*Y^i

    	#Adjust evaluation at point at infinity if necessary. 
    	if len(v1)==33:
		m+=v1[32]*(Y^32-Y)
    
    	#Map the <=32-degree polynomial from F_32[X] into an element of F_(2^325) via field_iso_desc  (by implicitely first mapping into an element of F_(32^65) and then changing to a representation in F_2^(325)).    
    	r=field_iso_desc(m,5,g,h,F,H,P,R)
    
    	return [r]



#psi_RMFE35: Given a list of elements w of F_(2^325), and a vector of values k<=99, computes psi-map of (k,325)_2-RMFE on each element and outputs the concatenation of the resulting vectors

def psi_RMFE35(w,k):
    
    if len(w)!=len(k):
	raise Exception("inputs to psi_RMFE24 must be of same length")
    for i in range(len(k)): 
	if k[i]>99:
		raise Exception("every coordinate on second input of psi_RMFE24 needs to be at most 99")
    B=[a^i for i in range(5)]
    data=FFTpreproc(5,B)
    
    res=[]
    for j in range(len(w)): 
    	#First change field representation to represent input as element of F_(32^65) and hence as a polynomial in F_32[X] of degree at most 64.
    	m=field_iso_asc(w[j],5,g,R)
    	m=list(m)
    
    	#Before applying the FFT we need to a polynomial of degree <=31. For this we take modulo X^32+X, as this does not modify evaluation in points of F_32:
    	hred=listsum(m[0:32],[0]+m[32:63])
    	hred=listsum(hred,[0]+m[63:])
        
    	#Apply FFT
    	w1=binaryFFT(hred,5,B,data)
    
    	 
    	#Based on value of k, we adjust size of the output. If k is between 97 and 99 we need to add the evaluation in the point at infinity.

    	if k[j]>=97:
		if len(m)>64:
		    	w1.append(m[64])
		else:
			w1.append(0)    
    	else:
		upper=(k[j]+2)//3	
		del w1[upper:]
   	
    	#Apply psi from (3,5)_2-RMFE to each element of resulting vector.
    	r=[]
    	for i in range(len(w1)):
    		r=r+invmap35(w1[i])
    
    	#Adjust size of output.
    	del r[k[j]:]

	#Concatenate this to global vector.
        res=res+r
    return res
    	
