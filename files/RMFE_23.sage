#!/usr/bin/env sage

from FFTpreproc import *
from FFTa import * 
from field_iso import *


##Implementation of a (18,51)_2 RMFE, as a concatenation of a (2,3)_2 RMFE and a (9,17)_8 RMFE. Functionality is extended to deal with (k,51)_2 RMFEs meaning that: if k<18, phi can be applied to an input (F_2)^k by first padding with zeros to get a vector in (F_2)^18; and if k>18, then the input vector to phi is splitted in subvectors of length 18 (and a remainder block of length k') and the map phi is applied to each vector. In the case of the map psi, it receives a vector of elements in the larger field H=F_(2^51) and applies psi:H->F_2^k to each coordinate, where the k's are specified by the user (they need to be k<=18) and can be different.

##Functions:
#map23: map phi of the (2,3)_2-RMFE: takes as input a vector of (F_2)^2, outputs a single element in the field F=F_8. 
#invmap23: map psi of (2,3)_2-RMFE: takes as input a single element in the field F=F_8, outputs a vector in (F_2)^2.
#phi_RMFE23: map phi of (k,51)_2-RMFE with "extended functionality": takes as input a vector in F_2^k, outputs a vector in H^m, where H=F_{2^51} (field of 2^51 elements) and m=k/18+1.
#psi_RMFE23: map psi of (k,51)_2-RMFE with "extended functionality": takes as input a vector in H^m, where H=F_{2^51} (field of 2^51 elements), and a vector of integers (k_1,...,k_m) with k_i<=18, and outputs a vector in F_2^(k_1+k_2+...+k_m), where each block is the result of applying psi to the i-th component of the input H-vector.  


H.<c>=GF(2^51,modulus="primitive")
h=H.modulus()

F.<a>=GF(2^3)
f=F.modulus()

P.<X>=PolynomialRing(GF(2))
R.<Y>=PolynomialRing(F)

g=R(h).factor()[0][0]


#map23: phi-map of (2,3)_2-RMFE

def map23(v): 
    return v[0]+(v[0]+v[1])*a

#invmap23: psi-map of (2,3)_2-RMFE

def invmap23(d): 
    if d!=0:
	p=d.polynomial(X)
	return [p(0),p(1)]
    else:
	return [0,0]
    return D


#phi_RMFE23: If given a binary vector of length k<=18, computes phi-map of (k,51)_2-RMFE. Else, split in blocks of length 18 (and a remainder block of length k') and compute phi-map of (18,51)_2-RMFE on each block (and (k',51)_2-RMFE on the last one). Outputs a list of elements in GF(2^51)

def phi_RMFE23(v): 
    if len(v)>18:
        #If the vector is of length >18, then split the vector in blocks of 18, apply map to each block.
	w=[]
        number_blocks=len(v)//18+1
	for i in range(number_blocks-1):
		w.append(v[18*i:18*i+18])
	w.append(v[18*(number_blocks-1):])
	res=[]
        for j in range(len(w)):
		res=res+phi_RMFE23(w[j])
	return res
    
    else:
    	#First step, split the binary vector in blocks of two cordinates (fill in left-over block with zeros), apply (2,3)_2-RMFE to each block. 
    	k=len(v)
    	odd=k%2
    
    	if odd:
    		v.append(0)
		l=1
    	else:
		l=0

    	v1=[]
    	for i in range((k+l)//2):
        	t=map23(v[2*i:2*i+2])
        	v1.append(t)            
    
    	#Second step, apply (k',17)_8-RMFE to result v1. Here k'=length(v1)<=9. Apply inverse FFT to the first 8 coordinates of v1 (filling in with zeros if necessary), obtaining an interpolating polynomial of degree <=7. 	    If length(v1)=9, then sum a multiple of X^8+X so that the X^8 coeficient is the 9-th coordinate. Note this 	operation does not change the evaluation of the polynomial in F_8. Map the <=8-degree polynomial into an      element of F_(8^17) and represent it as an element of F_(2^51) via field_iso_desc 
        while len(v1)<8:
        	v1.append(0)
    	
    	#We generate the preprocessing data for the FFT (TODO: having this as a preprocessing would only be useful if the precomputation would be used for more than one evaluation of phi_RMFE23, not the case currently).
    	B=[a^i for i in range(3)]
    	data=FFTpreproc(3,B)

    	#Apply inverse FFT
    	v2=invbinaryFFT(v1,3,B,data)
    	    
    	#Represent result as polynomial
    	m=0
    	for i in range(len(v2)):
        	m+=v2[i]*Y^i

    	#Adjust evaluation at point at infinity if necessary. 
    	if len(v1)==9:
		m+=v1[8]*(Y^8-Y)
     		
 
    	#Map the <=15-degree polynomial from F_16[X] into an element of F_(2^128) via field_iso_desc  (by implicitely first mapping into an element of F_(16^32) and then changing to a representation in F_2^(128)).    
    	r=field_iso_desc(m,3,g,h,F,H,P,R)
            
    	
    	return [r]





#psi_RMFE23: Given a list of elements w of F_(2^51), and a vector of values k<=18, computes psi-map of (k,51)_2-RMFE on each element and outputs the concatenation of the resulting vectors

def psi_RMFE23(w,k):
    ###print ('lengths', len(w), len(k))
    if len(w)!=len(k):
	raise Exception("inputs to psi_RMFE23 must be of same length")
    for i in range(len(k)): 
	if k[i]>18:
		raise Exception("every coordinate on second input of psi_RMFE need to be at most 18")
    B=[a^i for i in range(3)]
    data=FFTpreproc(3,B)
    ###print ('Applying psi to', w)
    res=[]
    for j in range(len(w)):
 	#First change field representation to represent input as element of F_(8^17) and hence as a polynomial in F_8[X] of degree at most 16.
    	m=field_iso_asc(w[j],3,g,R)
    	m=list(m)
    	###print ('After translating to F_8', m)
  
    	#Before applying the FFT we need to a polynomial of degree <=8. For this we take modulo X^8+X, as this does not modify evaluation in points of F_8:
    	hredi=listsum(m[0:8],[0]+m[8:15])
    	hred=listsum(hredi[:],[0]+m[15:])
    	###print ('After reduction', hred) 
    
    	#Apply FFT
    	w1=binaryFFT(hred,3,B,data)
    
        ###print ('After applying FFT we get', w1)
    
    	#Based on value of k, we adjust size of the output. If k is 17 or 18 we need to add the evaluation in the point at infinity.

    	if k[j]>=17:
		if len(m)>=17:
		    	w1.append(m[16])
		else:
			w1.append(0)
    	else:
		upper=(k[j]+1)//2	
		del w1[upper:]
    	###print ('After adjusting we get', w1)

    	#Apply psi from (2,3)_2-RMFE to each element of resulting vector.
    	r=[]
    	for i in range(len(w1)):
    		r=r+invmap23(w1[i])
    
    	#Adjust size of output.
    	del r[k[j]:]
        #Concatenate this to global vector.
        res=res+r
    return res


