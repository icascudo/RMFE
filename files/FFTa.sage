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
        binv=data[-1][m1-m+1]
        K=1<<m-1
        M=K<<1
        
        g=comp_cX(f,B[-1])                                                                            
        b=Taylor(g,M)
        gzero=[b[i][0] for i in range(len(b))]
        gone=[b[i][1] for i in range(len(b))]
        
	u=binaryFFT(gzero,m-1,D,data)
        v=binaryFFT(gone,m-1,D,data)
        Li=lin(G)
                                                                                             
        w=[u[i]+Li[i]*v[i] for i in range(K)]
	w1=[w[i]+v[i] for i in range(K)]

        return w+w1 

def invbinaryFFT(w,m,B,data):
#Inverse of previous algorithm.
    if m==1:
        if len(w)==1:
            return [w[0],0]
        else:
            return [w[0],(w[0]+w[1])*B[0]^(-1)]
    else:
        r=1
        G,D=data[-m]
        Li=lin(G)
	K=1<<m-1
        M=K<<1
 	u=[]
        v=[]

	for i in range(K):
		v.append(w[K+i]-w[i])
		u.append(w[i]-Li[i]*v[i])        

        gzero=invbinaryFFT(u,m-1,D,data)
        gone=invbinaryFFT(v,m-1,D,data)

        f=combine(gzero,gone,B[-1]^(-1))
	f=list(f)
	for i in range(M-len(f)):
		f.append(0)

        return f


