#!/usr/bin/env sage

from generatormatrix import *

import os
import errno


def print_to_file(instance,data):
	switcher={'phi_matrix_elements_H': generatormatrixphi(instance), 'phi_matrix_bits':generatormatrixphi_element(instance), 'psi_matrix':generatormatrixpsi(instance), 'S_o_psi_matrix':generatormatrixS_o_psi(instance), 'h_poly':instance.h, 'phi_matrix_bits':generatormatrixphi_bits(instance), 'Ker_S_o_psi_basis': KerS_o_psi(instance), 'Ker_psi_basis':Kerpsi(instance), 'phi_of_one':phi_of_one(instance)}

	option=switcher.get(data, lambda: 'Second argument should be one of the following'+str(list(switcher)))

	filename = 'output_data/'+'RMFE'+'_'+str(instance.k)+'_'+str(instance.e)+'/'+data+'_'+str(instance.k)+'_'+str(instance.e)+'.txt'
	if not os.path.exists(os.path.dirname(filename)):
		try:
			os.makedirs(os.path.dirname(filename))
		except OSError as exc: 
			if exc.errno != errno.EEXIST:
				raise
	
	filet=open(filename, "w")
	
	filet.write(str(option))
		
	filet.close()	
		



#BUG TO BE CORRECTED! The following does not seem to work well when List_instances contains more than one instance, in some cases where k1=3. When List_instances only contains one instance, or several instances with k1=2, it seems to work well.

#List_instances=[twostepinstance(2,4,4,8),twostepinstance(2,3,8,16),twostepinstance(2,3,9,17),twostepinstance(2,4,8,16),twostepinstance(2,4,16,32),twostepinstance(3,5,16,32),twostepinstance(3,6,16,32),
#twostepinstance(3,8,16,32),twostepinstance(3,5,33,65)]

List_instances=[twostepinstance(3,5,33,65)]

List_data=['phi_matrix_bits', 'psi_matrix', 'S_o_psi_matrix', 'h_poly', 'phi_matrix_bits', 'Ker_S_o_psi_basis', 'Ker_psi_basis', 'phi_of_one', 'phi_matrix_elements_H']

for instance in List_instances:
	for data in List_data:
		print_to_file(instance,data)

