
# This file was *autogenerated* from the file RMFEmapsfrommatrices.sage
from sage.all_cmdline import *   # import sage library

_sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0)
from twostepRMFE import *
from twostepinstance import *
from generatormatrix import *


#On input data and instance returns file "output_data/RMFE_k_e/data_k_e.txt"
def file_from_data(data, instance):
	return 'output_data/'+'RMFE'+'_'+str(instance.k)+'_'+str(instance.e)+'/'+data+'_'+str(instance.k)+'_'+str(instance.e)+'.txt'


#Apply phi to vector v, from the generator matrix N of phi. The input length is n, and k is as in instance
def phi_from_gen_matrix(v,N,n,k):
	if n%k!=_sage_const_0 :
		v=v+[_sage_const_0 ]*(k-(n%k))
	L=[sum ((v[a+j]==_sage_const_1 )*N[j] for j in range(k)) for a in range(_sage_const_0 ,n,k)]
	return L


#Applies psi to the vectors of a list V(vectors have length n), from the generator matrix N of phi. Used in [CG21] to compute Phi row-wise on matrix.
def phi_of_list_from_gen_matrix(V,N,n,k):
	return [phi_from_gen_matrix(V[i],N,n,k) for i in range(len(V))]


#Applies psi to all unit vectors of length n. Used in [CG21] to compute Phi row-wise on identity matrix
def phi_of_unit_vectors(N,n,k):
	result=[]
	for i in range(n):
		u=[GF(_sage_const_2 )(_sage_const_0 )]*i+[GF(_sage_const_2 )(_sage_const_1 )]+[GF(_sage_const_2 )(_sage_const_0 )]*(n-i-_sage_const_1 )
		result.append(phi_from_gen_matrix(u,N,n,k))
	return result


#Apply psi to vector v, from the generator matrix M of psi. The *output* length is n.
def psi_from_gen_matrix(v,M,n): 
	res=[]
	for i in range(len(v)):
		vt=v[i]._vector_()
		res=res+list(sum ((vt[j]==_sage_const_1 )*M[j] for j in range(len(vt))))
	return res[_sage_const_0 :n]


#Samples linear GF(2)-combination of coordinates vector"
def samplematrix(v):
	return sum (GF(_sage_const_2 ).random_element()*v[i] for i in range(len(v)))


#Samples at random from Ker psi [CG21]
def sampleKerpsi(instance):		
	return samplematrix(Kerpsi(instance))	

	
#Samples at random from Ker S o psi [CG21] 
def sampleKerS_o_psi(instance):		
	return samplematrix(KerS_o_psi(instance))


		




