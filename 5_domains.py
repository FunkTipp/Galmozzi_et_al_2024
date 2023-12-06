import pickle
import numpy as np
import os
import time
import json
#you need to install curl first


###CATH server domain search####
#dictionary of amino acid seqeunces
#str:str
#gene_name:'AMINOACIDSEQUENCEONELETTERCODE'
path = '/file/path/'
#with open(path + "aa_seq/aa_seq_MC4100.pick", "rb") as handle:
#	aa_seq = pickle.load(handle)

test_dict = {'aspC':'MFENITAAPADPILGLADLFRADERPGKINLGIGVYKDETGKTPVLTSVKKAEQYLLENETTKNYLGIDGIPEFGRCTQELLFGKGSALINDKRARTAQTPGGTGALRVAADFLAKNTSVKRVWVSNPSWPNHKSVFNSAGLEVREYAYYDAENHTLDFDALINSLNEAQAGDVVLFHGCCHNPTGIDPTLEQWQTLAQLSVEKGWLPLFDFAYQGFARGLEEDAEGLRAFAAMHKELIVASSYSKNFGLYNERVGACTLVAADSETVDRAFSQMKAAIRANYSNPPAHGASVVATILSNDALRAIWEQELTDMRQRIQRMRQLFVNTLQEKGANRDFSFIIKQNGMFSFSGLTKEQVLRLREEFGVYAVASGRVNVAGMTPDNMAPLCEAIVAVL'}

outfile_1 = {}
starttime = time.time()
#iterate over genes
for p, (gene, seq) in enumerate(test_dict.items()):
	#1. write the sequence into a fasta file
	A = open(path + '5_temp_1.fasta', 'w')
	A.write('fasta=>'+gene+'\n'+seq)
	A.close()
	#2. submit the job to the CATH server
	os.system('curl -o "'+path+'5_temp_2.json" -w "\n" -s -X POST -H "Accept: application/json" --data-binary "@'+path+'5_temp_1.fasta" http://www.cathdb.info/search/by_funfhmmer')
	#3. extract the jobnumber form the recieved job file
	with open(path + '5_temp_2.json') as json_file:
		B = json.load(json_file)
	id_ = B['task_id']
	#write the code for checking your submission status
	os_in_B = 'curl -o "'+path+'5_temp_3.json" -w "\n" -s -X GET -H "Accept: application/json" http://www.cathdb.info/search/by_funfhmmer/check/'+id_
	print(p, gene)
	#check every second if the sumbmission is complete
	#usually a gene takes 10-20 seconds
	c=0
	while True:
		time.sleep(10 - ((time.time() - starttime) % 10))
		c+=10
		if c%10 == 0:
			print('- waiting for ' + str(c) + ' seconds -')
		os.system(os_in_B)
		#some genes don't give a job status file (for what ever reason),
		#so this try/except command is a safety that the program runs trough all genes
		try:
			with open(path + '5_temp_3.json') as json_file:
				C = json.load(json_file)
		except:
			break
		if C['success'] == 1:
			break
	#get the results
	os_in_C = 'curl -o "'+path+'5_temp_4.json" -w "\n" -s -X GET -H "Accept: application/json" http://www.cathdb.info/search/by_funfhmmer/results/'+id_
	os.system(os_in_C)
	#save the results (only the Top 50 are given)
	#note that genes with no results don't give an outputfile
	try:
		with open(path + '5_temp_4.json') as json_file:
			D = json.load(json_file)
	except:
		continue
	outfile_1[gene] = D['funfam_scan']['results'][0]['hits']
#save outputfile 1 with all annotated domains
with open(path + "5_CATH_all_out.pick", "wb") as handle:
	pickle.dump(outfile_1, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open(path + "5_CATH_all_out.pick", "rb") as handle:
	outfile_1 = pickle.load(handle)
####domain sorting after alignment rank#####
outfile_2 = {}
for gene, data in outfile_1.items():
	temp = []
	#the gene data in the outfile_1 is already sorted after the highest ranked domain, so iteration starts\
	#with the best fit
	for p, i in enumerate(data):
		for pp, j in enumerate(i['hsps']):
			#exclude domains shorter or equal than 10 aa
			if j['length'] > 10:
				take = True
				query_range = np.arange(j['query_start'], j['query_end'])
				#only assign new domains if higher ranked domains are not jet assigend at these positions and \
				#allow the assignment of a new domain if the overlap to already assigned domains is below 5%
				for ppp in temp:
					match_range = np.arange(ppp[1], ppp[2])
					compare = np.count_nonzero(np.isin(query_range, match_range))
					if compare/len(match_range) > 0.05:
						take = False
				if take == True:
					#CATH ID, domain start, domain end
					temp.append([i['match_cath_id']['id'], j['query_start'], j['query_end']])
	if len(temp) > 0:
		sorted_temp = sorted(temp, key = lambda x: int(x[1]))
		outfile_2[gene] = sorted_temp
with open(path + "5_CATH_assigned_out.pick", "wb") as handle:
	pickle.dump(outfile_2, handle, protocol=pickle.HIGHEST_PROTOCOL)
