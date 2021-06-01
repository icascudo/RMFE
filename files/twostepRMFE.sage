#!/usr/bin/env sage

from FFTpreproc import *
from FFTa import * 
from field_iso import *
from twostepinstance import *
from polyintereval import *

#Computation of the composition of two RMFEs over F_2. Takes an instance (k_1,e_1,k_2,e_2) and computes diverse maps of the composition of (k_1,e_1)_2, (k_2,e_2)_{2^e_1} RMFEs, where k_1 = 2 or 3, 2k_1-1 <= e_1, k_2 <= 2^(e_2)+1, 2k_2-1<=e_2. The composition is a (k,e)_2-RMFE, where k=k_1k_2, e=e_1e_2
#phi_first applies phi-map of (k_1,e_1)_2-RMFE to a vector v of F_2^(k_1)
#psi_first applies psi-map of (k_1,e_1)_2-RMFE to an element in F=F_{2^e_1}=F_{m1}
#phi_RMFE_exact applies phi-map of the composition (the (k,e)_2-RMFE) to a vector with coordinates in F_2^(k).The output is an element in H=F_{2^e}.

#phi_RMFE applies phi-map of the composition (the (k,e)_2-RMFE) to a vector with coordinates in F_2^(k'). Here k' does not need to be k. If k'<k, the function fills in input with zeros until it is in F_2^k and applies phi. Ik k'>k, the function splits the vectir in blocks of k coordinates (plus possibly one vector of k'<k coordinates that it fills with zeros as above) and applies phi to each block. The output therefore is a vector of elements in H=F_{2^e}, represented as a list
#psi_RMFE applies psi-map of the composition (the (k,e)_2-RMFE) to a vector of elements with coordinates in F_{2^e}. It will apply the map to each coordinate and concatenate the results in a vector over F_2, given as a list. However, we introduce an input vk which is a list of numbers <=k of the same length of the input vector. The i-th coordinate of vk tells us the desired dimension of the image of psi on the i-th coordinate of the input. The reason is simply to make it consistent with phi_RMFE, meaning that we intend to apply it as vk=(k,k,....,k,k') where k'<=k. If we didn't define this vk, we would always get an output in F_2^k x...x F_2^k, while our phi_RMFE is more general.

#phi_first: phi-map of (k1,e1)_2-RMFE

def phi_first(v,instance):
	a=instance.a
	if instance.k1==2: 	 
  		return v[0]+(v[0]+v[1])*a
	elif instance.k1==3:
		return v[0]+(-v[0]+v[1]-v[2])*a+v[2]*a^2

#psi_first: psi-map of (k1,e1)_2-RMFE

def psi_first(d,instance):
	P=instance.P
	X=instance.X
	if d!=0:
		p=d.polynomial(X)
		r=[p(0),p(1)]
		if instance.k1==3:
			r.append(p[4])
		return r		
	else:
		return ([0]*instance.k1)



#phi_RMFE_exact: Phi map for (k,e)_2-RMFE when input is of size exactly k

def phi_RMFE_exact(v,instance): 
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

    	#Split the binary vector in blocks of k1 cordinates, apply (k1,e1)_2-RMFE to each block. 
 	v1=[]
   	
	for i in range(k2):
		t=phi_first(v[k1*i:k1*i+k1],instance)
		v1.append(t)            
   	#Second step, apply (k',e2)_2^e1-RMFE to result v1.
    	
    	#Interpolate
	p=interpol(v1,instance)
    	    
    	#Now we have a polynomial in F[X]=F_{2^e1}[X] of degree <=e2. We map this into an element of F_(2^e) via field_iso_desc  (by implicitely first mapping into an element of F_(m1^e2) and then changing to a representation in F_(2^e)).  Recall m1=2^e1
	r=field_iso_desc(p,instance)
        
    	
	return r


#phi_RMFE: Applies phi-map RMFE to a vector of arbitrary length. Output is *always* a list.

def phi_RMFE(v,instance):
	k=instance.k
	w=[]
	number_blocks=(len(v)-1)//k+1
	for i in range(number_blocks-1):
		w.append(v[k*i:k*i+k])
	w.append(v[k*(number_blocks-1):]+[0]*(k*number_blocks-len(v)))
	res=map(lambda x:phi_RMFE_exact(x,instance),w)
	return res
    
	




#psi_RMFE: Psi map for (k,e)_2-RMFE

def psi_RMFE(w,vk,instance):
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
    ###print ('lengths', len(w), len(k))
	if len(w)!=len(vk):
		raise Exception("inputs to psi_RMFE must be of same length")
	for i in range(len(vk)): 
		if vk[i]>k:
			raise Exception("every coordinate on second input of psi_RMFE need to be at most k")
	B=[a^i for i in range(e1)]
	data=FFTpreproc(e1,B)
    	
	res=[]
	for j in range(len(w)):
 	#First change field representation to represent input as element of F_(m1^e2) and hence as a polynomial in F[X] of degree at most e2. Recall m1=2^e1
		p=field_iso_asc(w[j],instance)
		p=list(p)
    		#print ('After translating to F', p)
		
    		#Before applying the FFT we need to a polynomial of degree <=m1. For this we take modulo X^m1+X, as this does not modify evaluation in points of F=F_m1:
		hredi=listsum(p[0:m1],[0]+p[m1:2*m1-1])
		hred=listsum(hredi[:],[0]+p[2*m1-1:])
    		#print ('After reduction', hred)
		
    
    		#Apply FFT
		w1=eval(hred,instance)
    
        	#print ('After applying FFT we get', w1)
    
    		#At this point we have a vector w1 of m1 elements in F and want to apply psi_first to each of these coordinates (obtaining k1 coordinates in F2 at each of these m1). If we apply it directly we get an output of m1*k1 coordinates, buit this may not coincide with vk[j], because: first, the theory also allows to have one more point, the point at infinity so in fact vk[j] may be as large as (m1+1)*k1; second, it also maybe that vk[j]<m1*k1. So based on value of vk, we adjust size of the output. If vk is larger than m1*k1 we need to add the evaluation in the point at infinity. If it is smaller we adjust so that the output will have vk[j] coordinates

		if vk[j]>m1*k1:
			if len(p)>2*m1:
				w1.append(p[2*m1])
			else:
				w1=w1+([0]*(vk[j]-len(p)))
		else:
			
			upper=((vk[j]-1)//k1) + 1	
			del w1[upper:]
		
    	###print ('After adjusting we get', w1)

    	#Apply psi_first to each element of resulting vector.
		r=[]
		for i in range(len(w1)):
			r=r+psi_first(w1[i],instance)
    
    	#Adjust size of output.
		del r[vk[j]:]
        #Concatenate this to global vector.
		res=res+r
	return res


