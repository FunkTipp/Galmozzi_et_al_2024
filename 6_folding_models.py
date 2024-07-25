import numpy as np
import pickle
import urllib.request
import json
import os
from scipy.spatial import distance
import matplotlib.pyplot as plt

def coord_array_maker(infile):
	coord_array = []
	pLDDT_array = []
	for D in infile[1:]:
		if D == '#':
			break
		E = D.split()
		aa_num = int(E[8])-1
		#atom coordinates
		X = float(E[10])
		Y = float(E[11])
		Z = float(E[12])
		coord_array.append([aa_num, X, Y, Z])
		pLDDT_array.append(float(E[14]))
	return np.array(coord_array), np.array(pLDDT_array)

def get_contacts(coord_dict_in, dist, gap, pLDDT, pLDDT_clip):
    #calculate distances between all atom
	all_dists = distance.squareform(distance.pdist(coord_dict_in[:,1:], metric='euclidean'))
    #calculate distances between residues	
	all_gaps = coord_dict_in[:,0] - coord_dict_in[:,0][np.newaxis].T
	all_pLDDT = np.amin(np.array(np.meshgrid(pLDDT,pLDDT)), axis=0)
	#sort out contacts of interest
	cond_ids = np.where((all_dists <= dist) & (all_gaps >= gap) & (all_pLDDT >= pLDDT_clip))
	res_array = np.array([coord_dict_in[cond_ids[0],0], coord_dict_in[cond_ids[1],0]]).T
	#sum up contacting atoms to residues
	set_res_array, count = np.unique(res_array, return_counts=True, axis=0)
	#array row = [resA, resB, strengh]
	return np.concatenate((set_res_array, count[np.newaxis].T), axis=1).astype(np.int_)

def cont_PAE(contact_map, PAE):
	out = np.zeros(contact_map.shape[0])
	for i in np.arange(contact_map.shape[0]):
		out[i] = PAE[contact_map[i, 0], contact_map[i, 1]]
	return out

def res_saturation(contact_map, last_emerged):
    #find contacts of category II and III
	stable_contact_ids = np.where((contact_map[:,0] <= last_emerged) & \
						   (contact_map[:,1] <= last_emerged))[0]
	stable_contact_ids_Kmodel2 = np.where((contact_map[:,0] <= last_emerged) & \
						   (contact_map[:,1] <= last_emerged) & \
						   (contact_map[:,3] != 1))[0]
	unsat_interD_contact_ids = np.where((contact_map[:,0] <= last_emerged) & \
						   (contact_map[:,1] > last_emerged) & \
						   (contact_map[:,3] == 1))[0]
	unsat_interD_contact_ids_Kmodel2addition = np.where((contact_map[:,0] <= last_emerged) & \
						   (contact_map[:,1] <= last_emerged) & \
						   (contact_map[:,3] == 1))[0]
	unsat_intraD_contact_ids = np.where((contact_map[:,0] <= last_emerged) & \
						   (contact_map[:,1] > last_emerged) & \
						   (contact_map[:,3] == 2))[0]
	#sum up interaction strength for every residues in each category
	#1: all stable contacts
	#2: all stable contacts for DnaK model 2
	#3: all interdomain contacts
	#4: additional interdomain conatcs for DnaK model 2
	#5: all intradomain contacts
	res_mat = np.zeros((last_emerged, 5))
	for res in np.arange(last_emerged):
		stable_res_ids = np.where(contact_map[stable_contact_ids,:2] == res)[0]
		res_mat[res, 0] = np.sum(contact_map[stable_res_ids, 2])
		stable_res_KM2_ids = np.where(contact_map[stable_contact_ids_Kmodel2,:2] == res)[0]
		res_mat[res, 1] = np.sum(contact_map[stable_res_KM2_ids, 2])
		unsat_inter_res_ids = np.where(contact_map[unsat_interD_contact_ids,:2] == res)[0]
		res_mat[res, 2] = np.sum(contact_map[unsat_inter_res_ids, 2])
		unsat_inter_res_KM2add_ids = np.where(contact_map[unsat_interD_contact_ids_Kmodel2addition,:2] == res)[0]
		res_mat[res, 3] = np.sum(contact_map[unsat_inter_res_KM2add_ids, 2])
		unsat_intra_res_ids = np.where(contact_map[unsat_intraD_contact_ids,:2] == res)[0]
		res_mat[res, 4] = np.sum(contact_map[unsat_intra_res_ids, 2])
	return res_mat

path = '/file/path/'
####get PAE data for all mmcif files####
for p, file in enumerate(os.listdir(path + '')):
	if file.endswith(".cif"):
		try:
			urllib.request.urlretrieve("https://alphafold.ebi.ac.uk/files/"+file[2:-13]+"-predicted_aligned_error_v2.json", \
					path + '/6_PAE_temp.json')
			A = open(path + '/6_PAE_temp.json')
			B = json.load(A)[0]
			C = np.zeros((max(B['residue1']), max(B['residue1'])), dtype=np.int16)
			C[np.array(B['residue1'])-1, np.array(B['residue2'])-1] = np.array(B['distance'])
			with open(path + "6_PAE_"+file[:-13]+".pick", "wb") as handle:
				pickle.dump(C, handle, protocol=pickle.HIGHEST_PROTOCOL)
		except:
			print(file)
			continue

####extract contact data from mmcif files####
#load CATH domain file
with open(path + "5_CATH_assigned_out.pick", "rb") as handle:
	CATH = pickle.load(handle)
contact_maps_out = {}
#iterate over all mmcif files
for p, file in enumerate(os.listdir(path + '')):
	if file.endswith(".cif"):
		#read the mmcif file
		inflie = open(path + '' + file).read()
		#get gene name
		gene_A = inflie.index('_ma_target_ref_db_details.gene_name')
		gene_B = inflie[gene_A:].split('\n')[0]
		gene = gene_B.split(' ')[-1]
		#get gene length
		gene_length_A = inflie.index('_ma_target_ref_db_details.seq_db_align_end')
		gene_length_B = inflie[gene_length_A:].split('\n')[0]
		gene_length = int(gene_length_B.split(' ')[-1])
		#get contacts
		coord_A = inflie.index('_atom_site.pdbx_sifts_xref_db_res')
		coord_B = inflie[coord_A:].split('\n')
		coord_array, pLDDT = coord_array_maker(coord_B)
		#contact requirements
		# <= 5 angstroms distance
		# >= 6 residues distance on the peptide chain
		# am pLDDT >= 70
		contacts_p = get_contacts(coord_array, 5, 6, pLDDT, 70)
		#implement domain structure
		# add another column of zeros to the contact map 
		contacts = np.c_[contacts_p, np.zeros(len(contacts_p), dtype=np.int8)]
		if gene not in CATH:
			continue	
		#sort domains from N to C terminus
		ddx = {}
		for k in CATH[gene]:
			ddx[k[1]] = [k[0], k[1], k[2]]
		#find contact residues within domains
		for p, (l, d) in enumerate(sorted(ddx.items())):
			for pp, (ll, dd) in enumerate(sorted(ddx.items())):
				A = np.isin(contacts[:, 0], np.arange(d[1], d[2]))
				B = np.isin(contacts[:, 1], np.arange(dd[1], dd[2]))
				C = A*B
				# inter domain contacts get the TAG 1 
				if p != pp:
					contacts[C, 3] = 1
				# intra domain contacts get the TAG 2
				if p == pp:
					contacts[C, 3] = 2
		#load predicted alginment error files (PAE)
		try:
			with open(path + "6_PAE_"+file[:-13]+".pick", "rb") as handle:
				PAE = pickle.load(handle)
			PAE[PAE<=5] = 1
			PAE[PAE>5] = 0
		except:
			PAE = np.ones((gene_length, gene_length))
		#exclude contacts with PAE > 5 angstrom
		PAE_mask = cont_PAE(contacts, PAE)
		conts_PAE = contacts[PAE_mask == 1, :]
		# add contact map and gene length to outputdict 
		contact_maps_out[gene] = [conts_PAE, gene_length]
#save contact maps
with open(path + "6_contact_maps.pick", "wb") as handle:
	pickle.dump(contact_maps_out, handle, protocol=pickle.HIGHEST_PROTOCOL)

####chaperone binding calculation####
#iterate over contact maps
master = {}
for g_count, (gene, contacts_length) in enumerate(contact_maps_out.items()):
    #get contacts and gene length
	contact_map = contacts_length[0]
	gene_length = contacts_length[1]
	#create chaperone binding outfile
	#row 0 = TF binding model
 	#row 1 = DnaK binding model 1
	#row 2 = DnaK binding model 2
	chaperone_binding_outfile = np.zeros((gene_length, 3))
	#iterate over translation states
	for ts in np.arange(gene_length):
		#correct for the tunnel
		last_emerged = ts-30
		if last_emerged < 0:
			continue
		#get emerged residue saturation
		residue_mat = res_saturation(contact_map, last_emerged)
		if len(residue_mat[:, 0]) == 0:
			continue

		####TF instant folding model####
		chaperone_binding_outfile[ts, 0] = np.sum(residue_mat[:, 4])

		####TF delayed folding model####
		#calculate unfolded C-terminus for the delayed folding model
		stable_res_A = residue_mat[:, 0]
		stable_res_A[stable_res_A < 20] = 0
		A = np.cumsum(stable_res_A[::-1])
		B = np.where((A[::-1] > 750))[0]
		if len(B) == 0:
			continue
		compact_ids_C = np.amax(B)
		if compact_ids_C == 0:
			continue
		chaperone_binding_outfile[ts, 1] = np.sum(residue_mat[:compact_ids_C, 4])
  
		####DnaK model 1####
		chaperone_binding_outfile[ts, 2] = np.sum(residue_mat[:compact_ids_C, 2])
  
		####DnaK model 2####
		K2_stable_res_A = residue_mat[:, 1]
		K2_stable_res_A[K2_stable_res_A < 20] = 0
		K2_A = np.cumsum(K2_stable_res_A[::-1])
		K2_B = np.where((K2_A[::-1] > 750))[0]
		if len(K2_B) == 0:
			continue
		K2_compact_ids_C = np.amax(K2_B)
		if K2_compact_ids_C == 0:
			continue
		chaperone_binding_outfile[ts, 3] = np.sum(residue_mat[:K2_compact_ids_C, 2])+np.sum(residue_mat[:K2_compact_ids_C, 3])
	master[gene] = chaperone_binding_outfile
with open(path + "6_chaperone_binding_out.pick", "wb") as handle:
	pickle.dump(master, handle, protocol=pickle.HIGHEST_PROTOCOL)

