import sys

#listsum- Takes two lists v,w, computes the componentwise sum. If v and w are of unequal length, it completes the smallest one with zeroes.
def listsum(v,w): 
	m=len(v)
	n=len(w)
	if n<m:
		r=map(lambda x,y:x+y, v, w+[0]*(m-n))
	else:
		r=map(lambda x,y:x+y, w, v+[0]*(n-m))
	return r

def conc(v,w):
	return v[:]+w[:]
	
def shift(v,m):
	return [0]*m+v

#nextpowertwo- Takes an integer and outputs the smallest power of two, which is larger than n.
def nextpowertwo(n): 
	return (n>1)*2**(len(bin(n-1))-2)+(n==1)

#bitsizet- Given integers n, t, it computes K power of two such that tK<n<=2*tK. Outputs K and M=tK.
def bitsizet(n,t):
	c=(n-1)//t
	K=(c>0)*2**(len(bin(c))-3)
	M=t*K                                                                      
	return K, M

#lin- Given a list H of length l, returns [sum H[j]*i[j]: i in {0,...,2^l-1}] where i[j] is the j-th bit of the binary decomposition of i.
def lin(H): 
	l=len(H)
	if l==1:
		return [0,H[0]]
	else:
		V=lin(H[:l-1])                                                                
		h=H[-1]
		W=map(lambda x: h+x, V) 		                                     
        	return V+W
#complexity 2^l sums    




