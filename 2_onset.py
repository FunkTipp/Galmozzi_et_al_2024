import numpy as np
import pickle

def num_value(chaperone_CI):
    binding_HC = np.where(chaperone_CI[:, 0] > 1.5)[0]
    binding_LC = np.where((chaperone_CI[:, 1] < 1.5) & (chaperone_CI[:, 0] > 1.5))[0]
    no_binding = np.where(chaperone_CI[:, 0] < 1.5)[0]
    binding = np.zeros(chaperone_CI.shape[0])
    binding[binding_HC] = 1
    binding[binding_LC] = 0.2
    binding[no_binding] = 0
    return binding
 
####get chaperone CI files#####
#for simplicity the same demo dataset is used for both replicates
path = '/file/path/'
with open(path + "1_CI_chap_CI.pick", "rb") as handle:
	chaperone_replicate_1 = pickle.load(handle)
with open(path + "1_CI_chap_CI.pick", "rb") as handle:
	chaperone_replicate_2 = pickle.load(handle)

####iterate over genes#####
binding_onet = {}
for gene in chaperone_replicate_1:
	#calculate binding sites and onset
	num_value_repliacte_1 = num_value(chaperone_replicate_1[gene])
	num_value_repliacte_2 = num_value(chaperone_replicate_1[gene])
	num_value_average = np.mean([num_value_repliacte_1, num_value_repliacte_1], axis=0)
	binding_sites = np.where(num_value_average >= 0.6)[0]

	binding_patches = np.array_split(binding_sites, np.flatnonzero(np.diff(binding_sites)>1)+1)
	valid_patches = []
	for i in binding_patches:
		if len(i) >= 6:
			valid_patches.extend(list(i))
	if len(valid_patches) == 0:
		onset = None
	else:
		onset = np.amin(valid_patches)
	binding_onet[gene] = onset

####Save file#####
with open(path + '2_chaperone_onsets.pick', 'wb') as handle:
	pickle.dump(binding_onet, handle, protocol=pickle.HIGHEST_PROTOCOL)

