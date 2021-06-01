
from lowlevel import *


#Some operations on polynomials, implemented fully in python for efficiency.


#comp_cX- Given a polynomial f represented as a list, and a field element c, outputs the list of coefficients of polynomial f(c*X), the composition of f and c*X.
def comp_cX(f,c):
	n=nextpowertwo(len(f))
	if n==0:
		return []
	elif n==1:
		return [f[0]]
	else:
		g1=f[0:n//2] 
		g2=f[n//2:]
		h1=comp_cX(g1,c) 
		h2=comp_cX(g2,c) 
		d=c**(n//2)
		r=[0]*(n//2)+map(lambda x:d*x, h2)
		s=listsum(h1, r)
		while s!=[] and s[-1]==0:
			s.pop(-1)
	return s

# comp_XX2- Given a polynomial f represented as a list, computes composition of f with (c*X)^2+c*X=c^2*X^2+c*X. Outputs the list of coefficients of the polynomial
#This is done recursively using that if n/2 is a power of two, then (c^2X^2+cX)^(n/2)=c^n*X^n+c^(n/2)*X^(n/2) 

def comp_XX2(f,c): 
	n=nextpowertwo(len(f))
	if n==0:
		return []
	elif n==1:
		return [f[0]]
	else:
		g1=f[0:n//2] 
		g2=f[n//2:]
		h1=comp_XX2(g1,c) #h1 and h2 have length ?
		h2=comp_XX2(g2,c) 
		d=c**(n//2)
		e=d**2
		r1=[0]*(n//2)+map(lambda x:d*x, h2)
#[d*h2[i] for i in range (len(h2))]
		r2=[0]*n+map(lambda x:e*x, h2)
#[e*h2[i] for i in range (len(h2))]
		s=listsum(h1, listsum(r1,r2))
		while s!=[] and s[-1]==0:
			s.pop(-1)
		return s


#Combine- Given two polynomials g0, g1 represented as lists, and a field element c, computes a list representing polynomial g0((c*X)^2+(c*X))+(c*X)*g1((c*X)^2+(c*X))
def combine(g0,g1,c):
	n=len(g0)
	m=len(g1)
	m0=comp_XX2(g0,c)
	m1=comp_XX2(g1,c)
	m1c=[0]+map(lambda x:c*x, m1)      
	result=listsum(m0,m1c)
	return result

