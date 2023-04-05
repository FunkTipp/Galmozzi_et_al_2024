import numpy as np
import pickle
import time

####parameter####
P_site_off = 15
FP_size_range = [20, 40]

####files####
main_path = "/Users/.../"
in_path = main_path + "2_sam_files/"
out_path = main_path + "3_RP_files/"
GFF3_file = main_path + "/0_genome_data/GCF_000005845.2_ASM584v2_genomic.gff"

####create reference dict####
CDS_dict = {}
read_collector_root = {}
with open(GFF3_file) as infile:
	for line in infile:
		if line[0] == '#':
			continue
		A = line.split('\t')
		if len(A) != 9:
			continue
		if A[2] != 'CDS':
			continue
		id = A[8].index('gene=')
		gene_p = A[8][id:].split(';')[0]
		gene = gene_p[5:]
		CDS_dict.setdefault(gene, [])
		CDS_dict[gene].append([A[0], int(A[3])-1, int(A[4]), A[6]])

####iterate over files#####
input_file_names = ['sample1_in']
output_file_names = ['sample1_out']
for file_name, out_name in zip(input, output):
	start_time = time.time()
	#create read collector array
	read_collector = np.zeros((2, 5000000))
	c=0
	#read sam file
	with open(in_path + file_name + '.sam') as infile:
		for line in infile:
			if line[0] == '@':
				continue
			c+=1
			if c%1e6==0:
				print(file_name, c/1e6, ' million reads collected,', int(time.time()-start_time), 'sec')
			A = line.split("\t")
			#get footprint length
			FP_length = abs(int(A[5][:-1]))
			#plus strand
			if A[1] == '0':
				if FP_size_range[0] <= FP_length <= FP_size_range[1]:
					r_pos = int(A[3])+len(A[9])-P_site_off+1
					read_collector[0, r_pos] += 1
			#minus strand
			if A[1] == '16':
				FP_length = abs(int(A[5][:-1]))
				if FP_range_paired[0] <= FP_length <= FP_range_paired[1]:
					r_pos = int(A[3])+P_site_off-1
					read_collector[1, r_pos] += 1

	####CDS assignment#####
	print(file_name, 'CDS assignment', int(time.time()-start_time), 'sec')
	all_reads = np.sum(read_collector)
	CDS_reads = 0
	nt_raw_reads = {}
	for gene, CC in CDS_dict.items():
		raw_reads = []
		for C in CC:
			if C[3] == '+':
				raw_reads = np.concatenate((raw_reads_pre, read_collector[0, C[1]:C[2]]))
			if C[3] == '-':
				raw_reads = np.concatenate((raw_reads_pre, np.flip(read_collector[1, C[1]:C[2]])))
		nt_raw_reads[gene] = raw_reads
		CDS_reads += np.sum(raw_reads)

	####Codon assignment#####
	print(file_name, 'Codon assignment', int(time.time()-start_time), 'sec')
	out_raw, out_rpm, out_rpkm = {}, {}, {}
	for gene, raw_reads in nt_raw_reads.items():
		rpm_reads = (raw_reads*1e6)/CDS_reads
		if len(raw_reads)%3 == 0:
			raw_reads_codon = np.sum(raw_reads.reshape(-1, 3), axis=1)
			rpm_reads_codon = np.sum(rpm_reads.reshape(-1, 3), axis=1)
		else:
			raw_reads_codon = np.sum(np.append(raw_reads, np.zeros(3-len(raw_reads)%3)).reshape(-1, 3), axis=1)
			rpm_reads_codon = np.sum(np.append(rpm_reads, np.zeros(3-len(raw_reads)%3)).reshape(-1, 3), axis=1)
		out_raw[gene] = raw_reads_codon
		out_rpm[gene] = rpm_reads_codon

		sum_raw_read = np.sum(raw_reads)
		gene_length = len(raw_reads)
		rpkm = (sum_raw_read * 1e9) / (CDS_reads * gene_length)
		rpm = (sum_raw_read * 1e6) / CDS_reads
		out_rpkm[gene] = [rpkm, rpm, sum_raw_read, gene_length]

	####Save files#####
	for X, n in zip([out_raw, out_rpm, out_rpkm], ['_raw_codon_P15', '_rpm_codon_P15', '_rpkm']):
		with open(out_path + out_name + n + ".pick", "wb") as handle:
			pickle.dump(X, handle, protocol=pickle.HIGHEST_PROTOCOL)
	print(file_name, 'done', int(time.time()-start_time), 'sec')
