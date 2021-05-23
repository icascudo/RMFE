#!/usr/bin/env sage

from a import *
from twostepinstance import *


#Generates the generator "matrix" (in this case a vector) over H where the i-th element is the image of the i-th unit vector. Elements are given in the basis 1,c,c^2,... where h(c)=0, for the representing polynomial h=instance.h
def generatormatrixphi_element(instance):
	Result=[]
	for i in range(instance.k):
		v=[GF(2)(0)]*(i)+[GF(2)(1)]+[GF(2)(0)]*(instance.k-i-1)
		Result.append(phi_RMFE(v,instance))
	return Result


# Generates the generator matrix of phi, for a given instance. The matrix is a k x e binary matrix, with respect to the representation of GF(2^e) by polynomial instance.h

def generatormatrixphi_bits(instance):
	Result=[]
	genmat=generatormatrixphi_element(instance)
	return matrix(GF(2), instance.k, instance.e, lambda i, j: genmat[i][0]._vector_()[j])

def gener	
