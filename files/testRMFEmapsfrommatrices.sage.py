
# This file was *autogenerated* from the file testRMFEmapsfrommatrices.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_10 = Integer(10); _sage_const_6 = Integer(6); _sage_const_5 = Integer(5); _sage_const_4 = Integer(4); _sage_const_9 = Integer(9); _sage_const_8 = Integer(8); _sage_const_10000 = Integer(10000); _sage_const_65 = Integer(65); _sage_const_17 = Integer(17); _sage_const_16 = Integer(16); _sage_const_33 = Integer(33); _sage_const_32 = Integer(32); _sage_const_50 = Integer(50)
from twostepRMFE import *
from twostepinstance import *
from generatormatrix import *
from RMFEmapsfrommatrices import *


#Tests for RMFEmapsfrommatrices.sage

#########TESTS###############


#####TEST that phi_from_gen_matrix, psi_from_matrix give correct results when used with stored generator matrices 'phi_matrix_elements_H' and 'psi_matrix' at files/outputdata/[subfolder corresponding to the given instance]

#Admits an instance as input. Syntax: "sage RMFEmapsfrommatrices.sage k1 e1 k2 e2" . Requires that the corresponding RMFE matrices are stored at output_data/RMFE_k_e
if len(sys.argv)==_sage_const_5 :
	inp=map(int, sys.argv[_sage_const_1 :_sage_const_5 ])
	instance=twostepinstance(inp[_sage_const_0 ],inp[_sage_const_1 ],inp[_sage_const_2 ],inp[_sage_const_3 ])
else:
	instance=twostepinstance(_sage_const_3 ,_sage_const_6 ,_sage_const_16 ,_sage_const_32 )

k=instance.k

#Input length, adjustable
n=_sage_const_50 

#Loading matrices
filename=file_from_data('phi_matrix_elements_H',instance)
filephi=open(filename, "r")
L=filephi.read()
N=vector(map(instance.H,L[_sage_const_2 :-_sage_const_2 ].split("], [")))
filephi.close()

filename=file_from_data('psi_matrix',instance)
filepsi=open(filename, "r")
L=filepsi.readlines()
Mrows=[]
v=L[_sage_const_0 ][_sage_const_1 :-_sage_const_2 ].split(" ")
for i in range(len(L)-_sage_const_1 ):
	Mrows.append(vector(map(GF(_sage_const_2 ),L[i][_sage_const_1 :-_sage_const_2 ].split(" "))))
Mrows.append(vector(map(GF(_sage_const_2 ),L[-_sage_const_1 ][_sage_const_1 :-_sage_const_1 ].split(" "))))
M=matrix(Mrows)
filepsi.close()

#Create random vectors of n bits. 
v=[GF(_sage_const_2 ).random_element() for i in range(n)]
w=[GF(_sage_const_2 ).random_element() for i in range(n)]

print ('Input vector v:')
print (v)
print ('Input vector w:')
print (w)

#Apply Phi to a and b
a= phi_from_gen_matrix(v, N, n,k)
b= phi_from_gen_matrix(w, N, n,k)

print ('v embedded in', a)
print ('w embedded in', b)

#Multiply images in H.
d=map(lambda x,y:x*y, a,b)
print ('Field product of embeddings:')
print (d)

#Apply Psi to product of images.
z=psi_from_gen_matrix(d,M,n)
print ('Psi of product:')
print (z)


#Check that this is indeed equal to the componentwise product zsup of v and w. 
zsup=map(lambda x,y:x*y, v,w)

print ('Supposed result(coordinatewise product of inputs):')
print (zsup)
print ('Is the obtained result as supposed?', z==zsup)

print '\n'

#Timing of application of phi:
n=_sage_const_10000 
print 'Timings for applications of phi to vectors of length', n
List_instances=[twostepinstance(_sage_const_2 ,_sage_const_4 ,_sage_const_4 ,_sage_const_8 ),twostepinstance(_sage_const_2 ,_sage_const_3 ,_sage_const_8 ,_sage_const_16 ),twostepinstance(_sage_const_2 ,_sage_const_3 ,_sage_const_9 ,_sage_const_17 ),twostepinstance(_sage_const_2 ,_sage_const_4 ,_sage_const_8 ,_sage_const_16 ),twostepinstance(_sage_const_2 ,_sage_const_4 ,_sage_const_16 ,_sage_const_32 ),twostepinstance(_sage_const_3 ,_sage_const_5 ,_sage_const_16 ,_sage_const_32 ),twostepinstance(_sage_const_3 ,_sage_const_6 ,_sage_const_16 ,_sage_const_32 ),
twostepinstance(_sage_const_3 ,_sage_const_8 ,_sage_const_16 ,_sage_const_32 ),twostepinstance(_sage_const_3 ,_sage_const_5 ,_sage_const_33 ,_sage_const_65 )]

v=[GF(_sage_const_2 ).random_element() for i in range(n)]

for instance in List_instances:
	filename=file_from_data('phi_matrix_elements_H',instance)
	filet=open(filename, "r")
	L=filet.read()
	N=vector(map(instance.H,L[_sage_const_2 :-_sage_const_2 ].split("], [")))
	t1=timeit('f=phi_from_gen_matrix(v, N, n, instance.k)', seconds=True)
	print 'time for (', instance.k, ',', instance.e, ')-RMFE ', t1
	filet.close()
#TODO: Timing for psi.


print '\n'
print '\n'

#####TEST that sampling with Kerpsi and KerS_o_psi (with sampling functions at files/RMFEmapsfrommatrices) give indeed values in the kernels of these functions.  

print 'Tests for sampling from kernels'

instance=twostepinstance(_sage_const_2 ,_sage_const_4 ,_sage_const_4 ,_sage_const_8 )

M=generatormatrixpsi(instance)


for i in range(_sage_const_10 ):
	
#Test 1: Sampling from Ker psi, and applying psi to the result. Should give 0 vector.
	v=sampleKerpsi(instance)
	print 'Sampling from Ker psi', v
	y=(psi_from_gen_matrix([v],M,instance.k))
	print 'Applying psi to result', y	
	print 'Is result 0 vector?', y==([_sage_const_0 ]*instance.k) 

#Test 2: Sampling from Ker S o psi, and applying psi to the result. Should give 0 value.
	w=sampleKerS_o_psi(instance)
	print 'Sampling from KerS o psi', w
	u=psi_from_gen_matrix([w],M,instance.k)
	print 'Applying psi to result', u
	#Applying S (sum of components) to result.
	z=(sum (u[i] for i in range(len(u))))
	print 'Applying S to that', z	
	print 'Is result 0?', z==_sage_const_0 
	print '\n'

#TODO: Timings, and reading Ker psi and Ker S o psi from text files at files/output_data.


