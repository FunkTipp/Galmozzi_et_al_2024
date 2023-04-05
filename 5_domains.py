import pickle
import os
import time
import json
#you need to install curl first


###CATH server domain search####
#dictionary of amino acid seqeunces
#str:str
#gene_name:'AMINOACIDSEQUENCEONELETTERCODE'
path = "/Users/.../"
with open(path + "aa_seq/aa_seq_MC4100.pick", "rb") as handle:
	aa_seq = pickle.load(handle)

outfile_1 = {}
starttime = time.time()
#iterate over genes
for p, (gene, seq) in enumerate(aa_seq.items()):
	#1. write the sequence into a fasta file
	A = open(path + 'temp_1.fasta', 'w')
	A.write('fasta=>'+gene+'\n'+seq)
	A.close()
	#2. submit the job to the CATH server
	os.system('curl -o "/Users/.../temp_2.json" -w "\n" -s -X POST -H "Accept: application/json" --data-binary "@/Users/.../temp_1.fasta" http://www.cathdb.info/search/by_funfhmmer')
	#3. extract the jobnumber form the recieved job file
	with open(path + 'temp_2.json') as json_file:
		B = json.load(json_file)
	id = B['task_id']
	#write the code for checking your submission status
	os_in_B = 'curl -o "/Users/.../temp_3.json" -w "\n" -s -X GET -H "Accept: application/json" http://www.cathdb.info/search/by_funfhmmer/check/'+id
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
			with open(path + 'temp_3.json') as json_file:
				C = json.load(json_file)
		except:
			break
		if C['success'] == 1:
			break
	#get the results
	os_in_C = 'curl -o "/Users/.../temp_4.json" -w "\n" -s -X GET -H "Accept: application/json" http://www.cathdb.info/search/by_funfhmmer/results/'+id
	os.system(os_in_C)
	#save the results (only the Top 50 are given)
	#note that genes with no results don't give an outputfile
	try:
		with open(path + '/temp_4.json') as json_file:
			D = json.load(json_file)
	except:
		continue
	#here I only save the CATH ID and the significance, but you can save all other data too
	outfile_1[gene] = []
	D_out = D['funfam_scan']['results'][0]['hits']
	for match in D['funfam_scan']['results'][0]['hits']:
		outfile_1[gene].append([match['match_cath_id']['id'], match['significance']])
#save outputfile 1 with all annotated domains
with open(path + "/CATH_all_out.pick", "wb") as handle:
	pickle.dump(outfile_1, handle, protocol=pickle.HIGHEST_PROTOCOL)

####domain sorting after alignment rank#####
outfile_2 = {}
for gene, data in outfile_1.items():
	temp = []
	#the gene data in the outfile_1 is already sorted after the highest ranked domain, so iteration starts\
	#with the best fit
	for p, i in enumerate(data[::-1]):
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
		outfile_2[gene] = temp
with open("/Users/.../CATH_assigned_out.pick", "wb") as handle:
	pickle.dump(outfile_2, handle, protocol=pickle.HIGHEST_PROTOCOL)
