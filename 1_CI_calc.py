import numpy as np
import pickle
from statsmodels.stats.proportion import proportion_confint

####parameter####
window = 45 #in codons

####iterate over files#####
path = '/Users/.../'
input_files = [['total_translatome_replicate_1', 'ip_replicate_1'], \
				['total_translatome_replicate_2', 'ip_replicate_2']]
output_file_names = ['CI_chap_repliacte_1', 'CI_chap_repliacte_2']
for experiment, out_name in zip(input_files, output_file_names):
	#load files
	with open(path + "1_data_files/" + str(experiment[0]) + '_raw_codon_P15.pick', 'rb') as f:
		total = pickle.load(f)
	with open(path + "1_data_files/" + str(experiment[1]) + '_raw_codon_P15.pick', 'rb') as f:
		ip = pickle.load(f)
	#count all reads
	outfile_CI, outfile_score = {}, {}
	sum_total = np.sum([np.sum(x) for x in total.values()])
	sum_ip = np.sum([np.sum(x) for x in ip.values()])
	normalization = sum_ip/sum_total
	#iterate over genes
	for p, gene in enumerate(total.keys()):
		total_gene = total[gene]
		ip_gene = ip[gene]
		#apply window
		total_window, ip_window = [], []
		for pos in np.arange(len(total_gene)):
			if pos <= np.floor(window/2):
				start = 0
			else:
				start = int(pos-np.floor(window/2))
			if pos > len(total_gene)-np.ceil(window/2):
				stop = len(total_gene)
			else:
				stop = int(pos+np.ceil(window/2))
			total_window.append(np.sum(total_gene[start:stop]))
			ip_window.append(np.sum(ip_gene[start:stop]))
		#calculate low CI
		CI = proportion_confint(ip_window, np.add(ip_window, total_window), alpha=0.05, method='agresti_coull')
		odds_CI_low = np.divide(CI[0], np.subtract(1, CI[0])) / norm
		#calculate high CI
		#clip high CIs to avoid division by zero
		CI_add = CI[1].copy()
		CI_add[CI_add == 1] = 0.9
		odds_CI_high = np.divide(CI[1], np.subtract(1, CI_add)) / norm
		outfile_CI[gene] = np.array([odds_CI_high, odds_CI_low]).T

		#calculate score
		outfile_score[gene] = np.amax(odds_CI_low)
	####Save files#####
	with open(path + "CI_calc/" + str(out_name) + '_CI.pick', 'wb') as handle:
		pickle.dump(outfile_CI, handle, protocol=pickle.HIGHEST_PROTOCOL)
	with open(path + "CI_calc/" + str(out_name) + '_score.pick', 'wb') as handle:
		pickle.dump(outfile_score, handle, protocol=pickle.HIGHEST_PROTOCOL)
