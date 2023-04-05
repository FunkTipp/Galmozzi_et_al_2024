import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.ticker import ScalarFormatter
import matplotlib
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

def mCI(a, confidence=0.95):
	n = len(a)
	m, se = np.mean(a), stats.sem(a)
	h = se * stats.t.ppf((1 + confidence) / 2., n-1)
	return m, m-h, m+h

####set parameter####
bin_window = 6 #in codons
plot_size = 500 #in codons

####load data####
path = "/Users/.../"
with open(path + "1_data_files/" + 'T1_raw_codon_P15.pick', 'rb') as f:
	total_1 = pickle.load(f)
with open(path + "1_data_files/" + 'T2_raw_codon_P15.pick', 'rb') as f:
	total_2 = pickle.load(f)
with open(path + "1_data_files/" + 'IP1_raw_codon_P15.pick', 'rb') as f:
	ip_1 = pickle.load(f)
with open(path + "1_data_files/" + 'IP2_raw_codon_P15.pick', 'rb') as f:
	ip_2 = pickle.load(f)

####get total reads####
ip = np.sum([np.mean([np.sum(ip_1[x]), np.sum(ip_2[x])]) for x in ip_1.keys()])
total = np.sum([np.mean([np.sum(total_1[x]), np.sum(total_2[x])]) for x in total_1.keys()])
seq_norm = total/ipn

###from start coding####
meta_dict = {x:[] for x in np.arange(0, plot_size+5, bin_window)}
###iterate over genes####
for g in groups:
    #average replicates
	ip_av = np.mean([ip_1[g], ip_2[g]], axis=0)
	tn_av = np.mean([total_1[g], total_2[g]], axis=0)
	#select position values
	for p in np.arange(0, plot_size+5, bin_window):
		#exclude larger positions
		if p < len(total_1):
			#exclude empty bins to avoit division by zero
			if np.sum(total_1[p:p+bin_window]) != 0:
				# add normalized values to the position wise meta dict
				meta_dict[p].append((np.sum(ip_av[p:p+bin_window])/\
				np.sum(tn_av[p:p+bin_window ]))*seq_norm1)

###average values for every position and caluculate 95 % CI####
meta_plot_dict = {}
for i, j in meta_dict.items():
	meta_plot_dict[i] = mCI(j)

###plot metagene profile####
fig, ax = plt.subplots(1, 1, figsize=(3.5, 2.2))
ax.tick_params(axis="both", which="major", labelsize=12)
plt.grid(which="major", axis="y", ls="--", zorder=0)
plt.xlim(0, plot_size)
plt.ylim(0.4, 5)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_yscale('log', base=2)
ax.yaxis.set_major_formatter(ScalarFormatter())
ids = sorted(meta_plot_dict.keys())
plt.plot(ids, [meta_plot_dict[x][0] for x in ids], color='orange', lw=2)
plt.fill_between(ids, [meta_plot_dict[x][1] for x in ids], \
					[meta_plot_dict[x][2] for x in ids], \
					lw=0, alpha=0.5, color='orange')
plt.ylabel('enrichment', fontsize=12)
plt.xlabel('distance from start (codons)', fontsize=12)
plt.yticks([0.5, 1, 2, 4])
plt.tight_layout()
plt.savefig(path + "out.pdf", dpi=300)
plt.show()
plt.close()
