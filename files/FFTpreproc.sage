#!/usr/bin/env sage

def FFTpreproc(m,B):
    data=[]
    Binv=[]
    assert (B[0].parent()).cardinality()==2^m
    C=B[:]
    for level in range(0,m-1):
        G=[]
        D=[]
        c=C[-1]^(-1)
        Binv.append(c)                                                                        
        for j in range(m-level-1):
            gam=c*C[j]                                                                      
            G.append(gam)     
            D.append(gam**2-gam)
        G1=G[:]
        D1=D[:] 
        datacurrent=[G1,D1]
        data.append(datacurrent)
        C=D
    data.append(Binv)         
    return data
