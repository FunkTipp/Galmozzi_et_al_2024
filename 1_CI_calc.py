import numpy as np
import pickle
from statsmodels.stats.proportion import proportion_confint
import time
####parameter####
window_transcripts = 15 #in codons
window_score = 45 #in codons

####iterate over files#####
path = '/file/path/'
input_files = [['0_IP_out', '0_total_out']]
output_file_names = ['1_CI_chap']
start_time = time.time()
for experiment, out_name in zip(input_files, output_file_names):
	#load files
	with open(path + str(experiment[0]) + '_raw_codon_P15.pick', 'rb') as f:
		total = pickle.load(f)
	with open(path + str(experiment[1]) + '_raw_codon_P15.pick', 'rb') as f:
		ip = pickle.load(f)
	#count all reads
	outfile_CI, outfile_score = {}, {}
	sum_total = np.sum([np.sum(x) for x in total.values()])
	sum_ip = np.sum([np.sum(x) for x in ip.values()])
	normalization = sum_ip/sum_total
	#iterate over genes
	for p, gene in enumerate(total.keys()):
		gene_length = len(total[gene])
		ip_window_transcripts, total_window_transcripts = [], []
		ip_window_score, total_window_score = [], []
		for pos in np.arange(gene_length):
			#apply transcrip window
			if pos <= np.floor(window_transcripts/2):
				start = 0
			else:
				start = int(pos-np.floor(window_transcripts/2))
			if pos > gene_length-np.ceil(window_transcripts/2):
				stop = gene_length
			else:
				stop = int(pos+np.ceil(window_transcripts/2))
			total_window_transcripts.append(np.sum(total[gene][start:stop]))
			ip_window_transcripts.append(np.sum(ip[gene][start:stop]))
			#apply score window
			if pos <= np.floor(window_score/2):
				start = 0
			else:
				start = int(pos-np.floor(window_score/2))
			if pos > gene_length-np.ceil(window_score/2):
				stop = gene_length
			else:
				stop = int(pos+np.ceil(window_score/2))
			total_window_score.append(np.sum(total[gene][start:stop]))
			ip_window_score.append(np.sum(ip[gene][start:stop]))

		#calculate low CI
		CI_transcripts = proportion_confint(ip_window_transcripts, np.add(ip_window_transcripts, total_window_transcripts), alpha=0.05, method='agresti_coull')
		odds_CI_low_transcript = np.divide(CI_transcripts[0], np.subtract(1, CI_transcripts[0])) / normalization

		CI_score = proportion_confint(ip_window_score, np.add(ip_window_score, total_window_score), alpha=0.05, method='agresti_coull')
		odds_CI_low_score = np.divide(CI_score[0], np.subtract(1, CI_score[0])) / normalization

		if len(odds_CI_low_score) <= 40:
			continue

		#calculate high CI
		#clip high CIs to avoid division by zero
		CI_add = CI_transcripts[1].copy()
		CI_add[CI_add == 1] = 0.9
		odds_CI_high = np.divide(CI_transcripts[1], np.subtract(1, CI_add)) / normalization
		outfile_CI[gene] = np.array([odds_CI_high, odds_CI_low_transcript]).T

		#calculate score
		outfile_score[gene] = np.amax(odds_CI_low_score[30:-10])
  
	####Save files#####
	print((time.time()-start_time), len(outfile_score))
	with open(path + str(out_name) + '_CI.pick', 'wb') as handle:
		pickle.dump(outfile_CI, handle, protocol=pickle.HIGHEST_PROTOCOL)
	with open(path + str(out_name) + '_score.pick', 'wb') as handle:
		pickle.dump(outfile_score, handle, protocol=pickle.HIGHEST_PROTOCOL)
