#!/usr/bin/env sage

from sys import *
from lowlevel import *
from taylor import *
from pypolyfunctions import *


#Implementation in Sage of Gao-Mateer10 first additive binary FFT algorithm, and its inverse.

def binaryFFT(f,m,B,data): #f is the polynomial of degree <n=2^m=q, B is a basis of Fq over F2 
	f=list(f)
	m1=len(data[-1])
	if m==1:
		return (f[0],f[0]+B[0]*f[1])
	else:
		G,D=data[-m]
		#binv=data[-1][m1-m+1]
		K=1<<m-1
		M=K<<1
        
		g=comp_cX(f,B[-1])
	                                                              
		b=Taylor(g,M)
		gzero=map(lambda x:x[0], b) 
		gone=map(lambda x:x[1], b) 
        
		u=binaryFFT(gzero,m-1,D,data)
		v=binaryFFT(gone,m-1,D,data)
		
		Li=lin(G)                                                                                        
		w=map(lambda x,y,z:x+y*z, u,Li,v)
		w1=map(lambda x,y:x+y, w,v) 

		return w+w1 

def invbinaryFFT(w,m,B,data):
#Inverse of previous algorithm.
	if m==1:
		if len(w)==1:
			return [w[0]]
		else:
			return [w[0],(w[0]+w[1])*B[0]^(-1)]
	else:
		G,D=data[-m]
		Li=lin(G)
		K=1<<m-1
		M=K<<1
		
		v=map(lambda x,y:x-y, w[K:2*K],w[0:K])
		u=map(lambda x,y,z: x-y*z, w[0:K],Li,v)		
		
		gzero=invbinaryFFT(u,m-1,D,data)
		gone=invbinaryFFT(v,m-1,D,data)
		
		f=combine(gzero,gone,B[-1]^(-1))
		
		return f+[0]*(M-len(f))



