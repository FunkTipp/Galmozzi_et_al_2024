#!/usr/local/bin/python3

import numpy as np
import pickle
from tslearn.metrics import dtw
from scipy.cluster.hierarchy import linkage, dendrogram, set_link_color_palette
from scipy.spatial.distance import squareform
import sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

sys.setrecursionlimit(4000)

def similarity_matrix_maker(mat, end):
	size = mat.shape[0]
	out = np.zeros((size, size))
	for p in np.arange(size):
		aa = mat[p]
		a = aa[30:end[p]-10]
		for pp in np.arange(size):
			if p < pp:
				bb = mat[pp]
				b = bb[30:end[pp]-10]
				out[p, pp] = dtw(a, b)
	return out

####load data####
#for simplicity the same demo dataset is used for both replicates
path = '/file/path/'
with open(path + "1_CI_chap_CI.pick", "rb") as handle:
	chaperone_CI_r1 = pickle.load(handle)
with open(path + "1_CI_chap_CI.pick", "rb") as handle:
	chaperone_CI_r2 = pickle.load(handle)

####construct similarity matrix####
#exclude genes which CO-IP themself and with repeated codon sites
excl_genes = ['tufA', 'tufB', 'dnaK', 'groEL', '']
#create a gene list to know the gene position in the similariy matrix
gene_list = np.array([x for x in chaperone_CI_r1.keys() \
    			if x not in excl_genes and len(chaperone_CI_r1[x][:, 1]) >= 50])
with open(path + '4_gene_list.pick', 'wb') as handle:
	pickle.dump(gene_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
#iterate over replicates
for exp, name in zip([chaperone_CI_r1, chaperone_CI_r2], \
					['4_C1', '4_C2']):
	longest_gene = 2360
    #create one array with low CI borders of every gene (in_mat) and one array with gene lengths
	in_mat = np.zeros((len(gene_list), longest_gene))
	gene_size = np.zeros(len(gene_list), dtype=np.int_)
	for p, g in enumerate(gene_list):
		r = np.arange(len(exp[g][:, 1]))
		gene_size[p] = len(exp[g][:, 1])
		in_mat[p, r] = exp[g][:, 1]
	#perform dtw on every gene pair (takes ≈ 1h for the E.coli genome)
	start_time = time.time()
	sim_mat = similarity_matrix_maker(in_mat, gene_size)
	print(time.time()-start_time)
	#save similarity matrix
	with open(path + str(name) + '_mat.pick', 'wb') as handle:
		pickle.dump(sim_mat, handle, protocol=pickle.HIGHEST_PROTOCOL)

####cluster genes based on similarity matrices####
####load cluster maps and gene_list ####
#for simplicity the same demo dataset is used for both replicates
with open(path + "4_gene_list.pick", "rb") as handle:
	gene_list = pickle.load(handle)
with open(path + "4_C1_mat.pick", "rb") as handle:
	C1 = pickle.load(handle)
with open(path + "4_C1_mat.pick", "rb") as handle:
	C2 = pickle.load(handle)

#average similarity matrices
X = np.mean([C1, C2], axis=0)
X+=X.T
#clip similarity matrix values
#this avoids in the GroEL IP for sdhA and frdA to from a third independent cluster
X[X > 500] = 500
X_flat = squareform(X)
#hierachical clustering
W = linkage(X_flat, method='ward')

####cluster plot####
fig, ax = plt.subplots(1, 3, figsize=(4, 6), sharey=True, gridspec_kw={'width_ratios': [3, 1, 0.25], 'wspace':0.1})
#Dendogram
D = dendrogram(W, orientation='left', color_threshold=500, ax=ax[0], above_threshold_color='black', labels=gene_list)
set_link_color_palette(['grey', 'orange'])
ax[0].spines["right"].set_visible(False)
ax[0].spines["top"].set_visible(False)
ax[0].spines["left"].set_visible(False)
ax[0].set_xticks([0, 500, 1000])
ax[0].set_yticks([]) #avoids plotting 4000 gene names
ax[0].set_xlabel('distance', fontsize=12)
ax[0].tick_params(axis="both", which="major", labelsize=12)
#Score plot
image = np.array([[np.amax(chaperone_CI_r1[g][:, 1]), np.amax(chaperone_CI_r2[g][:, 1])]*10 \
    				for g in D['ivl']]).reshape(len(gene_list)*10, 2)
image[image <= 0.001] = 0.001
im_in = np.log2(image)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ['whitesmoke','orange'])
imb = ax[1].imshow(im_in, interpolation='nearest', aspect='auto', cmap=cmap, vmin=0, vmax=4)
ax[1].spines["right"].set_visible(False)
ax[1].spines["top"].set_visible(False)
ax[1].spines["left"].set_visible(False)
ax[1].set_xticks([0, 1])
ax[1].tick_params(axis="both", which="major", labelsize=12)
#colorbar
ax[2].axis('off')
cbar = plt.colorbar(imb, pad=0.1, ticks=[0, 1, 2, 3, 4], ax=ax[1])
cbar.ax.set_yticklabels(['< 1', '2', '4', '8', '> 16'])
cbar.ax.set_ylabel('score', rotation=90, labelpad=-3, fontsize=12)
cbar.ax.tick_params(axis="both", which="major", labelsize=12)
plt.savefig(path + "4_cluster_plot_loCI.pdf", dpi=300)
plt.show()
