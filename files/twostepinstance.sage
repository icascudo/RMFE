#!/usr/bin/env sage

#We only deal for now with (k,e)_2-RMFEs constructed as the composition of a (k_1,e_1)_2-RMFE and a (k_2,e_2)_(2^e1)-RMFE.
#Note k=k1*k2, e=e1*e2
#Notations:
#m1=2^e1, m2=2^e 
#F is the intermediate field of size 2^e1, H is the final field of size 2^e.
#f is the irreducible polynomial in F_2[X] used for representing F, h is the irreducible polynomial in F_2[X] used for representing H.
#a is the class of X in F, i.e., a in F with f(a)=0. c is the class of X in H, i.e., c in H with h(c)=0
#P is the polynomial ring F_2[X], R is the polynopmial ring F[Y], note the different variable name, otherwise Sage is confused.
#g is the factorization of h in the polynomial ring R=F[Y] of the intermediate field (where h is not irreducible anymore).

class twostepinstance:
	def __init__(self,k1,e1,k2,e2):
		if k1!=2 and k1!=3:
			raise Exception ("First argument (k_1) must be 2 or 3")
		if e1<2*k1-1:
			raise Exception ("First argument (k_1) and second argument (e_1) must satisfy e_1>=2k_1-1")
		if k2>2**e1+1:
			raise Exception ("Second argument (e_1) and third argument (k_2) must satisfy k_2<=2^e_1+1")
		if e2<2*k2-1:
			raise Exception ("Third argument (k_2) and fourth argument (e_2) must satisfy e_2>=2k_2-1")
		self.k1 = k1
		self.k2 = k2
		self.e1 = e1
		self.e2 = e2
		self.k = k1*k2
		self.e = e1*e2
		self.m1 = 2**e1
		self.m2 = 2**(self.e) 
		self.F = GF(self.m1, names=('a',))
		(self.a,) = self.F._first_ngens(1)
		self.B = [self.a^i for i in range(self.e1)]       #Only needed if computing evaluation and interpolation via additive FFT.		
		#self.H, self.isoFH = self.F.extension(self.e2, map=True)   
                self.H = GF(self.m2, modulus="primitive", names=('c',))
		(self.c,) = self.H._first_ngens(1)
		self.P = PolynomialRing(GF(2), names=('X',))
		(self.X,) = self.P._first_ngens(1)
		self.R = PolynomialRing(self.F, names=('Y',))
		(self.Y,) = self.R._first_ngens(1)
		self.f = self.F.modulus()
		self.h = self.H.modulus()

#Next three atributes are needed for computation of field representation change (at field_iso).
		self.hR = self.R(self.h)                          #See h, initially a polynomial in P=F_2[X], as a polynomial in R=F[Y], where F is the intermediate field. 
		self.g = self.hR.factor()[0][0]			  #One factor of h as a polynomial in F[Y] (where it is not irreducible)
		self.isopol= self.R(self.hR/self.g)	          #Computing product of Frobenius-conjugates of g, not including g. 

		
