#!/usr/bin/env sage

#Computes a list data with the following structure: i-th component is [G,D] where G and D are the bases in Algorithm 2 of Gao-Mateer10 at recursion level m-i, i=0,...,m-2. The last component contains all [beta_j^{-1}] for levels j=m,...,2.
 
def FFTpreproc(m,B):
    data=[]
    Binv=[]
    assert (B[0].parent()).cardinality()==2^m
    C=B[:]
    for level in range(0,m-1):
        cc=C[-1]^(-1)
	Binv.append(cc)
	G=[cc*C[j] for j in range(m-level-1)]
	D=[G[j]**2-G[j] for j in range(m-level-1)]
                                                                       
        G1=G[:]
        D1=D[:] 
        datacurrent=[G1,D1]
        data.append(datacurrent)
        C=D
    data.append(Binv)         
    return data
