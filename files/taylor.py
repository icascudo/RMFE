
import sys
from lowlevel import *

#Implementation in Python of Taylor's-type algorithm in Gao-Mateer 10

#Taylor-This is case t=2, of GM10's algorithm. Takes list of coefficients f representing a polynomial and n power of two larger than its length, and decomposes degree-n polynomial f as f=f0+(X^2-X)*f1+(X^2-X)^2*f2+... ;Outputs  list [f0,f1,f2,...]
def Taylor(f,n):  
	V=[]
	if n<=2:
		f=f+[0]*(len(f)==1)
		V.append(f)
	else:
		[K,M]=bitsizet(n,2)
		a=M-K                                                                              
		f0=f[0:M]                                                                  
		f1=f[M:a+M]                                                                   
		f2=f[a+M:]                                                                    
		h=listsum(f1,f2)                                                                             
		g0=listsum(f0,[0]*K+h)                                                                       
		g1=listsum(h,[0]*a+f2)                                                                      
		V0=Taylor(g0,M)                                                                      
		V1=Taylor(g1,n-M)
		V=V+V0+V1
	return V
#Compl. S(n)=n/2+2*S(n/2), S(2)=0,   S(n)=n/2*log(n)





