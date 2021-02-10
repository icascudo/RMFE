import sys

#listsum- Takes two lists v,w, computes the componentwise sum. If v and w are of unequal length, it completes the samllest one with zeroes.
def listsum(v,w): 
    m=len(v)
    n=len(w)
    if n<m:
        r=[v[i]+w[i] for i in range(n)]+[v[j] for j in range(n,m)] 
    else:
        r=[v[i]+w[i] for i in range(m)]+[w[j] for j in range(m,n)]
    return r

def conc(v,w): #Concatenation of lists
    n=len(w)
    for i in range(n):
        v.append(w[i])
    return v

def shift(v,m): #Creates the list 0^m||v
    w=[]
    for i in range(m):
        w.append(0)
    for j in range(len(v)):
        w.append(v[j])
    return w

#nextpowertwo- Takes an integer and outputs the smallest power of two, which is largest than n.
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
        V=lin(H[:l-1])                                                                   #S(l-1)       =S(L/2)  
        W=[H[-1]+V[j] for j in range(len (V))]                                           #2^(l-1) sums =L/2
        return V+W
    #complexity 2^l sums    




