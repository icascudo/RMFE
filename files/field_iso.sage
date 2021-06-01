
import sys

#field_iso_desc- Given m in Fq[X]/(g), where q=2^d, outputs element in F2[X]/(h) via the "canonical" isomorphism. It is required that g is irreducible in F[X], h is irreducible in F2[X] and g is one of the d factors of the factorization of h over Fq[X]. For technical reasons, indeterminate variables in F2[X] and Fq[X] are called differently: X and Y respectively. R needs to be defined as Fq[Y], P as F2[X], c is a root of h defining the field F2[X]/(h).
#Isomorphism is computed via the CRT isomorphism between Fq[X]/(h) and Fq[X]/(g)xFq[X]/(g_1)...xFq[X]/(g_d-1) where the g_i's are the other factors of h, which are the conjugates of g by the Frobenius action on Fq
#We use the fact that the image of F2[X]/(h) via the isomorphism consists of an element and its Frobenius conjugates. So in order to compute the preimage of m, we compute the Frobenius actions m_i on it, and do CRT on (m,m_1,...,m_d-1).
#However, we use the fact that this is the sum of the CRT-preimage of (m,0,...,0) and its Frobenius-conjugates.



def field_iso_desc(m,instance):
	g=instance.g
	c=instance.c
	P=instance.P
	hR=instance.hR
	p=instance.isopol
  	R=instance.R	
	
#Computing CRT-preimage q of (m,0,0,...,0)
	q=xgcd(p,g)[1]*p*m
	
#Computing result as sum of Frobenius-conjugates mod h (i.e. coefficients of q are mapped to their trace)		                                                      
	r=q.map_coefficients(lambda z:z.trace())  				 
	s=r.mod(hR)
	t=P(s)                                           
	return t(c)


#field_iso_asc is the inverse of the function above, which is simply defined as u mod g, when u (originally an element of H, and hence representable as a polynomial in F_2[X]) is seen as a polynomial in F[Y].
def field_iso_asc(u,instance):
	return instance.R(u.polynomial(instance.Y)).mod(instance.g)


############UNUSED
#Alternative to field_iso_desc, seems slower 
def field_iso_desc_other(m,instance):
	X=instance.X
	g=instance.g
	c=instance.c
	P=instance.P
	hR=instance.hR
	p=instance.isopol
  	R=instance.R
	e1=instance.e1	


	Moduli=[g.map_coefficients(lambda z:z.frobenius(i)) for i in range(e1)]
	q=CRT([m]+[R(0)]*(e1-1), Moduli)	
	
#Computing result as sum of Frobenius-conjugates mod h (i.e. coefficients of q are mapped to their trace)		                                                      
	r=q.map_coefficients(lambda z:z.trace())  				 
	s=r.mod(hR)
	t=P(s.change_variable_name(X))                                           
	return t(c)

