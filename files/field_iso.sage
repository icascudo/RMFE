
#H.<c>=GF(2^51, modulus="primitive")
#h=H.modulus()

#F.<a>=GF(2^3)
#f=F.modulus()

#P.<X>=PolynomialRing(GF(2))
#R.<Y>=PolynomialRing(F)

#g=R(h).factor()[0][0]


#field_iso_desc- Given m in Fq[X]/(g), where q=2^d, outputs element in F2[X]/(h) via the "canonical" isomorphism. It is required that g is irreducible in Fq[X], h is irreducible in F2[X] and g is one of the d factors of the factorization of h over Fq[X]. For technical reasons, indeterminate variables in F2[X] and Fq[X] are called differently: X and Y respectively. R needs to be defined as Fq[Y], P as F2[X], c is a root of h defining the field F2[X]/(h).
#Isomorphism is computed via the CRT isomorphism between Fq[X]/(h) and Fq[X]/(g)xFq[X]/(g_1)...xFq[X]/(g_d-1) where the g_i's are the other factors of h, which are the conjugates of g by the Frobenius action on Fq
#We use the fact that the image of F2[X]/(h) via the isomorphism consists of an element and its Frobenius conjugates. So in order to compute the preimage of m, we compute the Frobenius actions m_i on it, and do CRT on (m,m_1,...,m_d-1).
#However, we use the fact that this is the sum of the CRT-preimage of (m,0,...,0) and its Frobenius-conjugates.


#Trying to decrease complexity for the case of RMFE_23, by precomputing part of the .
#H51.<c51>=GF(2^51, modulus="primitive")
#h51=H51.modulus()




def field_iso_desc(m,d,g,h,F,H,P,R):
    hR=R(h)
    Y=R.gen()
    a=R.base_ring().gen()
    #if h==h51:
    #    s=R(a^2*Y^49 + a*Y^47 + a*Y^46 + a^2*Y^45 + a^2*Y^44 + a*Y^43 + a^2*Y^42 + a*Y^41 + a*Y^40 + a*Y^39 + a^2*Y^38 + a^2*Y^37 + a*Y^35 + a*Y^34 + (a^2 + a)*Y^33 + a*Y^32 + a^2*Y^31 + a^2*Y^30 + a*Y^29 + (a^2 + a)*Y^26 + a*Y^22 + a^2*Y^21 + a^2*Y^19 + a*Y^18 + a^2*Y^15 + a^2*Y^14 + (a^2 + a)*Y^13 + (a^2 + a)*Y^11 + (a^2 + a)*Y^10 + a*Y^9 + a^2*Y^6 + a*Y^5 + a^2*Y^4 + (a^2 + a)*Y^3 + a^2*Y^2 + a^2 + 1)
	   
    p=R(hR/g)
    	#p=prod (g.map_coefficients(lambda z:z.frobenius(i)) for i in range(1,d)) #Computing product of Frobenius-conjugates of g, not including g.
    s=xgcd(p,g)[1]*p
    
    q=s*m                                                       #Computing CRT-preimage of (m,0,0,...,0)
    r=q.map_coefficients(lambda z:z.trace())  				     #Computing result as sum of Frobenius-conjugates mod h (i.e. coefficients of q are mapped to their trace)
    s=r.mod(hR)
    X=P.variable_name()
    t=P(s.change_variable_name(X))                                           #Horrible hack to get Sage to see this polynomial as a F2[X]-poly
    c=H.gen()
    return t(c)

def field_iso_asc(u,d,g,R):
    Y=R.variable_name()
    return R(u.polynomial(Y)).mod(g)


