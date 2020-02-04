#!/usr/bin/python
# python 2.6 or 2.7 

# privacy-rarefaction.kmers.v1.py
# Mathias Scharmann
# 2018-11-13


import sys
import numpy
import time
import os
import argparse
import subprocess
import multiprocessing as mp
import collections
import random
import gzip

"""
example usage:
python /privacy-rarefaction.kmers.v3.py --dump_dir ../ --CPUs 20 --sex_list fufulist.txt --n_resampling 20


PSEUDOCODE:

quantitative assessment of sex-specificity:
for n in 1,2,3, ... [minimum sample size of the two sexes]: # subsamplings
	repeat n_resampling times: # bootstraps	
		1. take random subsamples of size n from each males and females
			A. get observed number of male-specific and female-specific loci
				- female specific count: count the loci that do map female reads but not male reads
				- male specific count: count the loci that do map male reads but not female reads
				- store the two counts
			B. get permuted number of male-specific and female-specific loci
				- permute the sexes: male and female are mixed
				- count the loci that do map group A reads but not group B reads
				- count the loci that do map group B reads but not group A reads
				- store the two counts
	collect the observed and permuted sex-specific counts and compare the distributions: the p-value indicates the proportion of permuted sex-specific 	counts that are equal to or larger than the mean of the observed sex-specific count distribution 

qualitative assessment of sex-specificity:
for n in 1,2,3, ... [minimum sample size of the two sexes]:
	repeat n_resampling times:	
		1. take random subsamples of size n from each males and females
		2. female specific set: find the set of loci that do map female reads but not male reads
		3. male specific set: find the set of loci that do map male reads but not female reads
		4. store these two sets of loci
	then qualitatively evaluate the generated sets of loci:
		1. count how many times each locus occured in the male specific sets
		2. count how many times each locus occured in the female specific sets
		This generates a bootstrap support value for each locus.
		
"""


# checks for file existence:
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format(x)
		exit()
	x = str(x)
	return x


# checks for non-UNIX linebreaks:
def linebreak_check(x):

	if "\r" in open(x, "rb").readline():
		print "Error: classic mac (CR) or DOS (CRLF) linebreaks in {0}".format(x)
		exit()
	
# parses command line arguments
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--dump_dir", required=True, help="path of directory with jellyfish dump files", metavar="DIRECTORY")
	parser.add_argument("--dump_suffix", required=True, help="suffix of the jellyfish dump files, e.g. -RG.dump; if suffix includes dash (-), please wrap it in quotes and add a terminal space character (known bug in argparse")
	parser.add_argument("--sex_list", required=True,
        dest="sexlistfile", type=extant_file,
		help="name and path of the sex_list file: 1st column barcode/sample name separated by tab from second column indicating the sex; males = 1 and females = 2", metavar="FILE")
	parser.add_argument("--CPUs", required=True, help="number of CPUs to use in multiprocessing/parallel parts of script", metavar="INT")
	parser.add_argument("--o", required=True, help="name for the output files", metavar="STRING")
	parser.add_argument("--n_resampling", nargs='?', help="number of resampled datasets to be drawn for jacknifing over sample size & number of permutations for sex vs. sample ID, default = 200", metavar="INT", default = "200")
	parser.add_argument("--min_support_to_report_kmers", nargs='?', help="minimum boostrap support that a kmer must reach to be reported in the qualitative results, default = 0.5", metavar="FLOAT", default = "0.5")
	parser.add_argument("--min_count", nargs='?', help="minimum count per kmer to accepts as 'present', below kmer kmer will be classified as 'absent', default = 2", metavar="INT", default = "2")
	parser.add_argument("--min_shared", nargs='?', help="minimum number of samples that must share a kmer to consider it, default = 6", metavar="INT", default = "6")
	parser.add_argument("--min_stringency", nargs='?', help="smallest number of individuals to compete (stringency), default = 6", metavar="INT", default = "6")
	
	args = parser.parse_args()
	
	linebreak_check(args.sexlistfile)
	
	return args

#######

def check_congruence (sexlistfile, dump_folder, dump_suffix):
	
	sexlistsamples = []
	with open(sexlistfile, "r") as INFILE:
		for line in INFILE:
			if len(line) > 1:
				fields = line.strip("\n").split("\t")
				sexlistsamples.append(fields[0])
	sexlistsamples = set(sexlistsamples)
	
	for sample in sexlistsamples:
		extant_file(dump_folder + sample + dump_suffix)
	print "all samples in sex_list also have jellyfish dump files in the dump_dir, good to go!"
	
		


#########
		 		 
def read_sexlist (sexlist_file):
	
	sexdict = {}
	with open(sexlist_file, "r") as INFILE:	
		for line in INFILE:
			if len(line) > 1:
				sample = line.strip("\n").split("\t")[0] 
				sex = line.strip("\n").split("\t")[1] 
				if sex == "1":
					gender = "male"
				elif sex == "2":
					gender = "female"			
				try:
					sexdict[gender].append(sample)
				except KeyError:
					sexdict[gender] = [sample]
	
	for k, v in sexdict.items():
		print k + ":		"  + v[0] + "\n\t\t" + "\n\t\t".join(v[1:]) + "\n"
		 
	
	return sexdict

###############

def MT_read_inputs(all_samples, CPUs, dump_dir, dump_suffix, min_count):

	results = {}
	
	pool_size = min( [len(all_samples), CPUs ] )
	
	pool = mp.Pool( pool_size ) #use all available cores, if fewer samples than cores only as many as there are samples
	for sample in all_samples: # there is one worker for each 
		thefile = dump_dir + sample + dump_suffix
		print thefile
		results[sample] = pool.apply_async(read_jellyfish_dump, args=(thefile,min_count,) )
	pool.close()
	pool.join()


	# Get process results from the output queue
	#print output
	results1 = {}
	for i, result in results.items():
		results1[i] = result.get()

	return results1
	

def read_jellyfish_dump( infile, min_count ):
	
	count_data = []
	if infile.endswith(".gz"):
		F = gzip.open(infile, "r")
	else:
		F = open(infile, "r")


	for line in F:
		if line.startswith(">"):
			cnt = int( line.rstrip("\n").lstrip(">"))
		else:
			tag = line.rstrip("\n")
			if cnt >= min_count:
				count_data.append( [tag, 1] )
	
	F.close()

	return count_data
	

def privacy_rarefaction_core(dump_dir, sexdict, n_resampling, CPUs, min_support_to_report_kmers, min_count, output_name, dump_suffix, min_shared, min_stringency):
	
	all_samples = sorted(sexdict["male"] + sexdict["female"])
	
	sexdict_idx = {}
	for sex in sexdict.keys():
		for sample in sexdict[sex]:
			idx = all_samples.index(sample)
			try:
				sexdict_idx[sex].append(idx)
			except KeyError:
				sexdict_idx[sex] = [idx]
#	print sexdict_idx	
	
	all_samples = sorted(sexdict["male"] + sexdict["female"])
	
	# read mapping info to memory; only once!
	print "reading kmer counts (multiprocessing)"
	mapping_data_raw = MT_read_inputs(all_samples, CPUs, dump_dir, dump_suffix, min_count)
	print "Done reading kmer counts"
	mapping_data = {}
	for idx,sample in enumerate(all_samples):
		sample_mapped = mapping_data_raw[sample]
		for kmer_entry in sample_mapped:
			try:
				mapping_data[kmer_entry[0]][idx] = kmer_entry[1]
			except KeyError:
				mapping_data[kmer_entry[0]] = [0]*len(all_samples)
				mapping_data[kmer_entry[0]][idx] =  kmer_entry[1]
		del mapping_data_raw[sample] # gradually free up memory
#	print mapping_data
#	exit()
	print "total kmers:	", len(mapping_data.keys())
	min_shared = int(min_shared)
	kmer_list = [] # two lists replace dictionary as main data structure throughout
	mapped_list = []
	for kmer, mapped in mapping_data.items():
		if sum(mapped) >= min_shared: # kmers MUST occur in at least min_shared samples
			kmer_list.append( kmer )
			mapped_list.append( mapped )
	print "kmers with mapped reads passing count and shared thresholds:	", len( kmer_list )
	del mapping_data
	mapping_data = [ kmer_list, mapped_list ]		
	# get also a dictionary with the lsit indices of each kmer:
	mapping_data_idx_dict = {kmer:idx for idx, kmer in enumerate( kmer_list )  }
#	print mapping_data_idx_dict

	# make the stats, jacknifing over sample number numbers (100 replicate subsamples per jackknife-level)
	min_n = min( len(sexdict["male"]), len(sexdict["female"]) )

	# now the jackknife, multi-threaded!
	
	MT_return_dict = MT_resampling(min_stringency, min_n, sexdict_idx, mapping_data, n_resampling, CPUs, permute_sexes = False)
	
#	print MT_return_dict[3][1][1]
	
	# MT_return_dict structure:
	# keys: range(n_resampling)
	# values: [pres_abs_result, kmer_details]
	# pres_abs_result is a list of length range(1, min_n+1) "stringency levels" ; 
	#	for each stringency level, contains an element [male_specific, female_specific, male_total_loci, female_total_loci]
	# kmer_details is a dict; keys = confidence level i (pres-abs number of samples); values: [[male specific kmer IDs],[female specific kmer IDs]]
	
	# in the end, evaluate the resamplings by taking their mean of private RAD-loci per sex
	pres_abs_resampled_results = {}	
	cnt = 0
	for i in range(min_stringency, min_n+1):
		cnt += 1
		male_spec_resampled = [MT_return_dict[j][0][cnt-1][0] for j in range(n_resampling) ]
		female_spec_resampled = [MT_return_dict[j][0][cnt-1][1] for j in range(n_resampling) ] 
		total_kmers_m_resampled = [MT_return_dict[j][0][cnt-1][2] for j in range(n_resampling) ]
		total_kmers_f_resampled = [MT_return_dict[j][0][cnt-1][3] for j in range(n_resampling) ] 		

		male_specific = numpy.mean( male_spec_resampled )
		female_specific = numpy.mean( female_spec_resampled )
		male_specific_std = numpy.std( male_spec_resampled )
		female_specific_std = numpy.std( female_spec_resampled )
		male_specific_min = numpy.min( male_spec_resampled )
		female_specific_min = numpy.min( female_spec_resampled )
		total_kmers_m_mean = numpy.mean(total_kmers_m_resampled)
		total_kmers_m_std = numpy.std(total_kmers_m_resampled)
		total_kmers_f_mean = numpy.mean(total_kmers_f_resampled)
		total_kmers_f_std = numpy.std(total_kmers_f_resampled)


		pres_abs_resampled_results[i] = [male_specific, female_specific, male_specific_std, female_specific_std, male_specific_min, female_specific_min, total_kmers_m_mean, total_kmers_m_std, total_kmers_f_mean, total_kmers_f_std]
	
#	print pres_abs_resampled_results
		
	# in the end, evaluate the resampled loci by considering only those as truly specific that turn up as specific loci in 50% of re-/subsampling rounds:
	
	# clear prev. outputs files if present:
	with open("male_specific_candidates." + output_name + ".txt", "w") as OUTFILE:
		OUTFILE.write("subsample size per sex" + "\t" + "kmer_ID" + "\t" + "subsampling bootstrap support" + "\t" + "number of F with this kmer" + "\n")
	with open("female_specific_candidates." + output_name + ".txt", "w") as OUTFILE:
		OUTFILE.write("subsample size per sex" + "\t" + "kmer_ID" + "\t" + "subsampling bootstrap support" + "\t" + "number of M with this kmer" + "\n")
	
		
	consistently_specific_loci = {}
	print "qualitative evaluation of candidate sex-secific loci . . ."
	print "candidates passing bootstrap threshold:\nM F"
	for i in range(min_stringency, min_n+1):
#		print "get lists of loci from all resamplings"
		spec = [MT_return_dict[j][1][i][0] for j in range(n_resampling) ] # get lists of loci from all resamplings
#		print "flatten the 2dim list"
		spec_flat = [item for sublist in spec for item in sublist] # flatten the 2dim list
		a = set(spec_flat)		
		m_counts = count_item_occurence(spec_flat)
		good_male_specs = [kmer for kmer in a if m_counts[kmer]  >= min_support_to_report_kmers*n_resampling] # retain only those which occured n_resampling times
		
#		print m_counts
#		print "done males"
		
		spec = [MT_return_dict[j][1][i][1] for j in range(n_resampling) ] # get lists of loci from all resamplings
		spec_flat = [item for sublist in spec for item in sublist] # flatten the 2dim list
		a = set(spec_flat)
		f_counts = count_item_occurence(spec_flat)
		good_female_specs = [kmer for kmer in a if f_counts[kmer] >= min_support_to_report_kmers*n_resampling]  # retain only those which occured n_resampling times
		
		consistently_specific_loci[i] = [good_male_specs,good_female_specs]
	
		print len(good_male_specs), len(good_female_specs)
#		print mapping_data
#	print "got sex specific loci IDs that are consistent among 0.5 of subsampling rounds, outputting to file"		
#	print pres_abs_resampled_results
		with open("male_specific_candidates." + output_name + ".txt", "a") as OUTFILE:
			outlines = []
#			outlines.append( [ x+"\t"+str(i) for x in consistently_specific_loci[i][0] ] )
#			outlines.append( [ str(i) + "\t" +  x + "\t" + str((float(m_counts[x])/float(n_resampling))*100.0) for x in consistently_specific_loci[i][0] ] )		
			outlines.append( [ str(i) + "\t" +  x + "\t" + str((float(m_counts[x])/float(n_resampling))*100.0) + "\t" + str( sum([ mapping_data[1][ mapping_data_idx_dict[x]  ][idx] for idx in sexdict_idx["female"] ]) ) for x in good_male_specs ] )
			outlines = [ "\n".join(x) for x in outlines[:] ]
			OUTFILE.write( "\n".join(outlines) + "\n")
	
		with open("female_specific_candidates." + output_name + ".txt", "a") as OUTFILE:
			outlines = []
#			outlines.append( [ x+"\t"+str(i) for x in consistently_specific_loci[i][1] ] )
#			outlines.append( [ str(i) + "\t" +  x + "\t" + str((float(f_counts[x])/float(n_resampling))*100.0) for x in consistently_specific_loci[i][1] ] )
			outlines.append( [ str(i) + "\t" +  x + "\t" + str((float(f_counts[x])/float(n_resampling))*100.0) + "\t" + str( sum([ mapping_data[1][ mapping_data_idx_dict[x]  ][idx] for idx in sexdict_idx["male"] ]) ) for x in good_female_specs ] )
			outlines = [ "\n".join(x) for x in outlines[:] ]
			OUTFILE.write( "\n".join(outlines) + "\n")
	
	
	### now go for the null hypothesis that male and female are identical (permutation):
	print "testing null hypothesis that male and female are identical (permutation)"
	
	MT_return_dict = MT_resampling(min_stringency, min_n, sexdict_idx, mapping_data, n_resampling, CPUs, permute_sexes = True)
	pres_abs_resampled_results_perm = {}	
	cnt = 0
	for i in range(min_stringency, min_n+1):
		cnt += 1
		male_specific = [MT_return_dict[j][0][cnt-1][0] for j in range(n_resampling) ]
		female_specific = [MT_return_dict[j][0][cnt-1][1] for j in range(n_resampling) ]
		pres_abs_resampled_results_perm[i] = [male_specific, female_specific]

#	print pres_abs_resampled_results_perm
	
	# get p-value obs/permuted:
	# the p-value indicates the proportion of permuted sex-specific counts that are equal to or larger than the mean observed sex-specific count (the mean count from all SUBsampled observed data), i.e. the overlap of permuted vs. observed distributions.
	
	collected_results = {}
	for i in range(min_stringency, min_n+1):
		male_spec_obs = pres_abs_resampled_results[i][0]
		male_spec_obs_std = pres_abs_resampled_results[i][2]
		male_spec_obs_min = pres_abs_resampled_results[i][4]
		perm_greater_obs = len( [x for x in pres_abs_resampled_results_perm[i][0] if x >= male_spec_obs] )
		p_male = float(perm_greater_obs) / len(pres_abs_resampled_results_perm[i][0])
		
		female_spec_obs = pres_abs_resampled_results[i][1]
		female_spec_obs_std = pres_abs_resampled_results[i][3]
		female_spec_obs_min = pres_abs_resampled_results[i][5]
		perm_greater_obs = len( [x for x in pres_abs_resampled_results_perm[i][1] if x >= female_spec_obs] )
		p_female = float(perm_greater_obs) / len(pres_abs_resampled_results_perm[i][1])


		
		collected_results[i] = {"obs_means" : [male_spec_obs, female_spec_obs], "obs_std" : [male_spec_obs_std, female_spec_obs_std], "p_val" : [p_male, p_female], "total_kmers_stats" : pres_abs_resampled_results[i][6:] }
	
#	print collected_results
	
	## write to a file:
	with open("permutation_results." + output_name + ".txt", "w") as OUTFILE:
		outlines = [["n_samples_per_sex", "male_specific_mean", "male_specific_std", "p_val", "female_specific_mean", "female_specific_std", "p_val", "total_kmers_m_mean", "total_kmers_m_std", "total_kmers_f_mean", "total_kmers_f_std"]]
		for i in range(min_stringency, min_n+1):
			outlines.append( [str(x) for x in [i, collected_results[i]["obs_means"][0], collected_results[i]["obs_std"][0], collected_results[i]["p_val"][0], collected_results[i]["obs_means"][1], collected_results[i]["obs_std"][1], collected_results[i]["p_val"][1], collected_results[i]["total_kmers_stats"][0], collected_results[i]["total_kmers_stats"][1], collected_results[i]["total_kmers_stats"][2], collected_results[i]["total_kmers_stats"][3] ] ] )
		
		
		outlines = [ "\t".join(x) for x in outlines[:] ]
		OUTFILE.write( "\n".join(outlines)  + "\n")



def count_item_occurence(lst):
    
    res = collections.defaultdict(lambda: 0)
    for v in lst:
        res[v] += 1
    return res
	
def MT_resampling(min_stringency, min_n, sexdict_idx, mapping_data, n_resampling, CPUs, permute_sexes):
	
	print "resampling (multiprocessing) . . ."	
	results = {}
	
	pool = mp.Pool(CPUs) #use all available cores, otherwise specify the number you want as an argument
	for i in range(n_resampling): # there is one worker for each resampling round
		results[i] = pool.apply_async(pres_abs_MT, args=(min_stringency, min_n, sexdict_idx, mapping_data, permute_sexes))
	pool.close()
	pool.join()


	# Get process results from the output queue
	#print output
	results1 = {}
	for i, result in results.items():
		results1[i] = result.get()

	return results1

####

def get_pres_abs (mapping_data, males, females):
	
	# get the histogram of locus presence / absence: count for each sex
	# mapping_data is a dictionary; keys = kmers ; values = list of pres/abs (1/0) info for the samples; samples are represented by a fixed list index!
	# e.g. { '403848_L105': [0, 1, 0, 1] }
	
	pres_abs_data = {}
#	print mapping_data.keys()
	for i in range(len(mapping_data[1])):
		kmer = mapping_data[0][i]
		mappedlist = mapping_data[1][i]
		males_presence = sum( [ mappedlist[x] for x in males ] )	
		females_presence = sum( [ mappedlist[x] for x in females ] )
		pres_abs_data[kmer] = [males_presence, females_presence]
	
	return pres_abs_data


#####


def pres_abs_MT(min_stringency, min_n, sexdict_idx, mapping_data, permute_sexes):
	
#	print permute_sexes
	
	# random.seed() is necessary to ensure that random has a different seed in 
	# each thread; otherwise each thread will return the same random.choice! It uses system time to make the seed.
	random.seed()
	
	# pres_abs_result returns number (count) of "sex-specific" kmers for each stringency level i
	# for each stringency level, contains an element [male_specific, female_specific, male_total_loci, female_total_loci]
	pres_abs_result = []
	
	# a dictionary to hold the IDs of the reference kmers/loci/RADtags which are private to each sex
	# keys = confidence level i (pres-abs number of samples); values: [[male specific kmer IDs],[female specific kmer IDs]]
	kmer_details = {} 
	
	shared_among_all_samples = 0

	if permute_sexes == False:
		cnt = 0
		for i in range(min_stringency, min_n+1):
#			print i
			cnt += 1
			jackn_males = random.sample(sexdict_idx["male"], i)
			jackn_females = random.sample(sexdict_idx["female"], i)			
#			print jackn_males, jackn_females
		
			pres_abs_data = get_pres_abs (mapping_data, jackn_males, jackn_females)

			# [male_specific, female_specific, male_total_loci, female_total_loci]
			pres_abs_result.append([0,0,0,0])
			
			kmer_details[i] = [[],[]]
			
			# absence in exactly i samples is assured by the subsampling!
			# i.e. only i samples are scanned for pres-abs 	
			
			for loc in pres_abs_data.keys():
				if pres_abs_data[loc][0] == i:
					pres_abs_result[cnt-1][2] += 1 # counting overall presence: number of kmers present in all males (total number of kmers at this stringency)
					if pres_abs_data[loc][1] == 0:
						# if seen in exactly i males and absent in exactly i females				
						pres_abs_result[cnt-1][0] += 1
						kmer_details[i][0].append(loc)
						
				if pres_abs_data[loc][1] == i:
					pres_abs_result[cnt-1][3] += 1 # counting overall presence: number of kmers present in all females (total number of kmers at this stringency)
					if pres_abs_data[loc][0] == 0:
						# if seen in exactly i females and absent in exactly i males
						pres_abs_result[cnt-1][1] += 1
						kmer_details[i][1].append(loc)
				
				# counting overall presence: number of kmers present in all males and all females (total number of kmers at this stringency)			
							
	else:
		# the permutation option: kmer_details is left empty since meaningless; pres_abs_result returns number of "sex-specific" kmers for each stringency level i
		cnt = 0
		all_samples = sexdict_idx["female"] + sexdict_idx["male"]
		for i in range(min_stringency, min_n+1):
#			print i
			cnt += 1
			jackn_males = random.sample(sexdict_idx["male"], i)
			jackn_females = random.sample(sexdict_idx["female"], i)			
#			print jackn_males, jackn_females			
			pres_abs_data = get_pres_abs (mapping_data, jackn_males, jackn_females)

			pres_abs_result.append([0,0])

			for loc in pres_abs_data.keys():
				if pres_abs_data[loc][0] == i:
					if pres_abs_data[loc][1] == 0:				
						pres_abs_result[cnt-1][0] += 1
				if pres_abs_data[loc][1] == i:
					if pres_abs_data[loc][0] == 0:
						pres_abs_result[cnt-1][1] += 1
	
	return [pres_abs_result, kmer_details]
	

def get_pres_abs_per_kmer (dumpfile, min_count):

	count_data = []
	print dumpfile
	with open(dumpfile, "r") as F:	
		for line in F:
			if line.startswith(">"):
				cnt = int( line.rstrip("\n").lstrip(">"))
			else:
				tag = line.rstrip("\n")
				if cnt >= min_count:
					count_data.append( [tag, 1] )
				
	return count_data

		

######################## MAIN

args = get_commandline_arguments ()
	
check_congruence (args.sexlistfile, args.dump_dir, args.dump_suffix)

sexdict = read_sexlist (args.sexlistfile)

if int(args.min_stringency) > min( len(sexdict["male"]), len(sexdict["female"]) ):
	print "Error, min_stringency can not be greater than min number of samples in M or F"
	exit()

privacy_rarefaction_core(args.dump_dir, sexdict, int(args.n_resampling), int(args.CPUs), float(args.min_support_to_report_kmers), int(args.min_count), args.o, args.dump_suffix, args.min_shared, int(args.min_stringency))

print "Done!"
