from twostepRMFE import *
from twostepinstance import *
from generatormatrix import *
from RMFEmapsfrommatrices import *


#Tests for RMFEmapsfrommatrices.sage

#########TESTS###############


#####TEST that phi_from_gen_matrix, psi_from_matrix give correct results when used with stored generator matrices 'phi_matrix_elements_H' and 'psi_matrix' at files/outputdata/[subfolder corresponding to the given instance]

#Admits an instance as input. Syntax: "sage RMFEmapsfrommatrices.sage k1 e1 k2 e2" . Requires that the corresponding RMFE matrices are stored at output_data/RMFE_k_e
if len(sys.argv)==5:
	inp=map(int, sys.argv[1:5])
	instance=twostepinstance(inp[0],inp[1],inp[2],inp[3])
else:
	instance=twostepinstance(3,6,16,32)

k=instance.k

#Input length, adjustable
n=50

#Loading matrices
filename=file_from_data('phi_matrix_elements_H',instance)
filephi=open(filename, "r")
L=filephi.read()
N=vector(map(instance.H,L[2:-2].split("], [")))
filephi.close()

filename=file_from_data('psi_matrix',instance)
filepsi=open(filename, "r")
L=filepsi.readlines()
Mrows=[]
v=L[0][1:-2].split(" ")
for i in range(len(L)-1):
	Mrows.append(vector(map(GF(2),L[i][1:-2].split(" "))))
Mrows.append(vector(map(GF(2),L[-1][1:-1].split(" "))))
M=matrix(Mrows)
filepsi.close()

#Create random vectors of n bits. 
v=[GF(2).random_element() for i in range(n)]
w=[GF(2).random_element() for i in range(n)]

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
n=10000
print 'Timings for applications of phi to vectors of length', n
List_instances=[twostepinstance(2,4,4,8),twostepinstance(2,3,8,16),twostepinstance(2,3,9,17),twostepinstance(2,4,8,16),twostepinstance(2,4,16,32),twostepinstance(3,5,16,32),twostepinstance(3,6,16,32),
twostepinstance(3,8,16,32),twostepinstance(3,5,33,65)]

v=[GF(2).random_element() for i in range(n)]

for instance in List_instances:
	filename=file_from_data('phi_matrix_elements_H',instance)
	filet=open(filename, "r")
	L=filet.read()
	N=vector(map(instance.H,L[2:-2].split("], [")))
	t1=timeit('f=phi_from_gen_matrix(v, N, n, instance.k)', seconds=True)
	print 'time for (', instance.k, ',', instance.e, ')-RMFE ', t1
	filet.close()
#TODO: Timing for psi.


print '\n'
print '\n'

#####TEST that sampling with Kerpsi and KerS_o_psi (with sampling functions at files/RMFEmapsfrommatrices) give indeed values in the kernels of these functions.  

print 'Tests for sampling from kernels'

instance=twostepinstance(2,4,4,8)

M=generatormatrixpsi(instance)


for i in range(10):
	
#Test 1: Sampling from Ker psi, and applying psi to the result. Should give 0 vector.
	v=sampleKerpsi(instance)
	print 'Sampling from Ker psi', v
	y=(psi_from_gen_matrix([v],M,instance.k))
	print 'Applying psi to result', y	
	print 'Is result 0 vector?', y==([0]*instance.k) 

#Test 2: Sampling from Ker S o psi, and applying psi to the result. Should give 0 value.
	w=sampleKerS_o_psi(instance)
	print 'Sampling from KerS o psi', w
	u=psi_from_gen_matrix([w],M,instance.k)
	print 'Applying psi to result', u
	#Applying S (sum of components) to result.
	z=(sum (u[i] for i in range(len(u))))
	print 'Applying S to that', z	
	print 'Is result 0?', z==0
	print '\n'

#TODO: Timings, and reading Ker psi and Ker S o psi from text files at files/output_data.


