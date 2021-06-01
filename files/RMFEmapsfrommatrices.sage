from twostepRMFE import *
from twostepinstance import *
from generatormatrix import *

def file_from_data(data, instance):
	return 'output_data/'+'RMFE'+'_'+str(instance.k)+'_'+str(instance.e)+'/'+data+'_'+str(instance.k)+'_'+str(instance.e)+'.txt'


#Apply phi to vector, reading the generator matrix of phi from a previously saved file.
def phi_from_gen_matrix(v,N,n,k):
	if n%k!=0:
		v=v+[0]*(k-(n%k))
	L=[sum ((v[a+j]==1)*N[j] for j in range(k)) for a in range(0,n,k)]
	return L

def psi_from_gen_matrix(v,M,n):
	res=[]
	for i in range(len(v)):
		vt=v[i]._vector_()
		res=res+list(sum ((vt[j]==1)*M[j] for j in range(len(vt))))
	return res[0:n]






#########TESTS###############

#Admits an instance as input. Syntax: "sage RMFEmapsfrommatrices.sage k1 e1 k2 e2" . Requires that the corresponding RMFE matrices are stored at output_data/RMFE_k_e
if len(sys.argv)==5:
	inp=map(int, sys.argv[1:5])
	instance=twostepinstance(inp[0],inp[1],inp[2],inp[3])
else:
	instance=twostepinstance(3,6,16,32)

k=instance.k

#Input length, adjustable
n=3*k+2

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


'''
#Timing:
List_instances=[twostepinstance(2,4,4,8),twostepinstance(2,3,8,16),twostepinstance(2,3,9,17),twostepinstance(2,4,8,16),twostepinstance(2,4,16,32),twostepinstance(3,5,16,32),twostepinstance(3,6,16,32),
twostepinstance(3,8,16,32),twostepinstance(3,5,33,65)]

for instance in List_instances:
	filename=file_from_data('phi_matrix_elements_H',instance)
	filet=open(filename, "r")
	L=filet.read()
	N=vector(map(instance.H,L[2:-2].split("], [")))
	k=instance.k
	t1=timeit('f=phi_from_gen_matrix(v, N, n,k)', seconds=True)
	print t1
	filet.close()
'''
