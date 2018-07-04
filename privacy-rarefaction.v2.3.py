#!/usr/bin/python
# python 2.6 or 2.7 

# privacy-rarefaction.v2.2.py
# a python script that implements privacy-rarefaction to extract sex-specific RAD-tags from .bam alignment files
# Copyright (C) 2017, ETH Zurich, Mathias Scharmann
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 		
# If you use this code please cite:
#
# "Scharmann M, Grafe TU, Metali F, Widmer A. (2017) Sex-determination 
# and sex chromosomes are shared across the radiation of dioecious 
# Nepenthes pitcher plants. XXX"
# 	
# contact: mathias.scharmann[-at-]env.ethz.ch or msph52[-at-]gmail.com



import sys
import numpy
import time
import os
import argparse
import subprocess
import multiprocessing as mp
import collections
import random

"""
example usage:
python /privacy-rarefaction.v2.2.py --bam_dir ../ --CPUs 20 --sex_list fufulist.txt --n_resampling 20 --bam_suffix "-RG.bam "


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
	
	parser.add_argument("--bam_dir", required=True, help="path of directory with .BAM files", metavar="DIRECTORY")
	parser.add_argument("--bam_suffix", required=True, help="suffix of the .BAM files, e.g. -RG.bam; if suffix includes dash (-), please wrap it in quotes and add a terminal space character (known bug in argparse")
	parser.add_argument("--sex_list", required=True,
        dest="sexlistfile", type=extant_file,
		help="name and path of the sex_list file: 1st column barcode/sample name separated by tab from second column indicating the sex; males = 1 and females = 2", metavar="FILE")
	parser.add_argument("--CPUs", required=True, help="number of CPUs to use in multiprocessing/parallel parts of script", metavar="INT")
	parser.add_argument("--o", required=True, help="name for the output files", metavar="STRING")
	parser.add_argument("--n_resampling", nargs='?', help="number of resampled datasets to be drawn for jacknifing over sample size & number of permutations for sex vs. sample ID, default = 200", metavar="INT", default = "200")
	parser.add_argument("--min_support_to_report_loci", nargs='?', help="minimum boostrap support that a locus must reach to be reported in the qualitative results, default = 0.5", metavar="FLOAT", default = "0.5")
	parser.add_argument("--min_cov", nargs='?', help="minimum read coverage (depth) per contig to count as 'present', below contig will be calssified as 'absent', default = 1", metavar="INT", default = "1")
	
	args = parser.parse_args()
	
	linebreak_check(args.sexlistfile)
	
	return args

#######

def check_congruence (sexlistfile, bam_folder, bam_suffix):
	
	sexlistsamples = []
	with open(sexlistfile, "r") as INFILE:
		for line in INFILE:
			if len(line) > 1:
				fields = line.strip("\n").split("\t")
				sexlistsamples.append(fields[0])
	sexlistsamples = set(sexlistsamples)
	
	for sample in sexlistsamples:
		extant_file(bam_folder + sample + bam_suffix)
	print "all samples in sex_list also have bamfiles in the bam_dir, good to go!"
	
		


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

def privacy_rarefaction_core(bam_dir, sexdict, n_resampling, CPUs, bam_suffix, min_support_to_report_loci, min_cov, output_name):
	
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
	print "reading mapping information"
	mapping_data = {}
	for sample in all_samples:
		sample_mapped = get_pres_abs_per_contig (bam_dir + sample + bam_suffix, min_cov)
		for contig in sample_mapped.keys():
			try:
				mapping_data[contig].append(sample_mapped[contig])
			except KeyError:
				mapping_data[contig] = [ sample_mapped[contig] ]
	
	print "total contigs:	", len(mapping_data.keys())
	for contig, mapped in mapping_data.items():
		if sum(mapped) == 0:
			del mapping_data[contig]
	print "contigs with mapped reads passing coverage threshold:	", len(mapping_data.keys())
		
	# make the stats, jacknifing over sample number numbers (100 replicate subsamples per jackknife-level)
	min_n = min( len(sexdict["male"]), len(sexdict["female"]) )

	# now the jackknife, multi-threaded!
	
	MT_return_dict = MT_resampling(min_n, sexdict_idx, mapping_data, n_resampling, CPUs, permute_sexes = False)
	
	# MT_return_dict structure:
	# keys: range(n_resampling)
	# values: [pres_abs_result, contig_details]
	# pres_abs_result is a dict; keys = range(1, min_n+1) ; values: [count_male_specific, count_female_specific, count_total_loci_present] 
	# contig_details is a dict; keys = confidence level i (pres-abs number of samples); values: [[male specific contig IDs],[female specific contig IDs]]
	
	# in the end, evaluate the resamplings by taking their mean of private RAD-loci per sex
	pres_abs_resampled_results = {}	
	for i in range(1, min_n+1):
		male_spec_resampled = [MT_return_dict[j][0][i-1][0] for j in range(n_resampling) ]
		female_spec_resampled = [MT_return_dict[j][0][i-1][1] for j in range(n_resampling) ] 
		total_contigs_m_resampled = [MT_return_dict[j][0][i-1][2] for j in range(n_resampling) ]
		total_contigs_f_resampled = [MT_return_dict[j][0][i-1][3] for j in range(n_resampling) ] 		

		male_specific = numpy.mean( male_spec_resampled )
		female_specific = numpy.mean( female_spec_resampled )
		male_specific_std = numpy.std( male_spec_resampled )
		female_specific_std = numpy.std( female_spec_resampled )
		male_specific_min = numpy.min( male_spec_resampled )
		female_specific_min = numpy.min( female_spec_resampled )
		total_contigs_m_mean = numpy.mean(total_contigs_m_resampled)
		total_contigs_m_std = numpy.std(total_contigs_m_resampled)
		total_contigs_f_mean = numpy.mean(total_contigs_f_resampled)
		total_contigs_f_std = numpy.std(total_contigs_f_resampled)


		pres_abs_resampled_results[i] = [male_specific, female_specific, male_specific_std, female_specific_std, male_specific_min, female_specific_min, total_contigs_m_mean, total_contigs_m_std, total_contigs_f_mean, total_contigs_f_std]
	
#	print pres_abs_resampled_results
		
	# in the end, evaluate the resampled loci by considering only those as truly specific that turn up as specific loci in 50% of re-/subsampling rounds:
	
	# clear prev. outputs files if present:
	with open("male_specific_candidates." + output_name + ".txt", "w") as OUTFILE:
		OUTFILE.write("subsample size per sex" + "\t" + "contig_ID" + "\t" + "subsampling bootstrap support" + "\t" + "number of F mapping this contig" + "\n")
	with open("female_specific_candidates." + output_name + ".txt", "w") as OUTFILE:
		OUTFILE.write("subsample size per sex" + "\t" + "contig_ID" + "\t" + "subsampling bootstrap support" + "\t" + "number of M mapping this contig" + "\n")
	
		
	consistently_specific_loci = {}
	print "qualitative evaluation of candidate sex-secific loci . . ."
	print "candidates passing bootstrap threshold:\nM F"
	for i in range(1, min_n+1):
#		print "get lists of loci from all resamplings"
		spec = [MT_return_dict[j][1][i][0] for j in range(n_resampling) ] # get lists of loci from all resamplings
#		print "flatten the 2dim list"
		spec_flat = [item for sublist in spec for item in sublist] # flatten the 2dim list
		a = set(spec_flat)		
		m_counts = count_item_occurence(spec_flat)
		good_male_specs = [contig for contig in a if m_counts[contig]  >= min_support_to_report_loci*n_resampling] # retain only those which occured n_resampling times
		
#		print good_male_specs
#		print "done males"
		
		spec = [MT_return_dict[j][1][i][1] for j in range(n_resampling) ] # get lists of loci from all resamplings
		spec_flat = [item for sublist in spec for item in sublist] # flatten the 2dim list
		a = set(spec_flat)
		f_counts = count_item_occurence(spec_flat)
		good_female_specs = [contig for contig in a if f_counts[contig] >= min_support_to_report_loci*n_resampling]  # retain only those which occured n_resampling times
		
		consistently_specific_loci[i] = [good_male_specs,good_female_specs]
	
		print len(good_male_specs), len(good_female_specs)
#		print mapping_data
#	print "got sex specific loci IDs that are consistent among 0.5 of subsampling rounds, outputting to file"		
#	print pres_abs_resampled_results
		with open("male_specific_candidates." + output_name + ".txt", "a") as OUTFILE:
			outlines = []
#			outlines.append( [ x+"\t"+str(i) for x in consistently_specific_loci[i][0] ] )
#			outlines.append( [ str(i) + "\t" +  x + "\t" + str((float(m_counts[x])/float(n_resampling))*100.0) for x in consistently_specific_loci[i][0] ] )
			outlines.append( [ str(i) + "\t" +  x + "\t" + str((float(m_counts[x])/float(n_resampling))*100.0) + "\t" + str( sum([ mapping_data[x][idx] for idx in sexdict_idx["female"] ]) ) for x in good_male_specs ] )
			outlines = [ "\n".join(x) for x in outlines[:] ]
			OUTFILE.write( "\n".join(outlines) + "\n")
	
		with open("female_specific_candidates." + output_name + ".txt", "a") as OUTFILE:
			outlines = []
#			outlines.append( [ x+"\t"+str(i) for x in consistently_specific_loci[i][1] ] )
#			outlines.append( [ str(i) + "\t" +  x + "\t" + str((float(f_counts[x])/float(n_resampling))*100.0) for x in consistently_specific_loci[i][1] ] )
			outlines.append( [ str(i) + "\t" +  x + "\t" + str((float(f_counts[x])/float(n_resampling))*100.0) + "\t" + str( sum([ mapping_data[x][idx] for idx in sexdict_idx["male"] ]) ) for x in good_female_specs ] )
			outlines = [ "\n".join(x) for x in outlines[:] ]
			OUTFILE.write( "\n".join(outlines) + "\n")
	
	
	### now go for the null hypothesis that male and female are identical (permutation):
	print "testing null hypothesis that male and female are identical (permutation)"
	
	MT_return_dict = MT_resampling(min_n, sexdict_idx, mapping_data, n_resampling, CPUs, permute_sexes = True)
	pres_abs_resampled_results_perm = {}	
	for i in range(1, min_n+1):
		male_specific = [MT_return_dict[j][0][i-1][0] for j in range(n_resampling) ]
		female_specific = [MT_return_dict[j][0][i-1][1] for j in range(n_resampling) ]
		pres_abs_resampled_results_perm[i] = [male_specific, female_specific]

#	print pres_abs_resampled_results_perm
	
	# get p-value obs/permuted:
	# the p-value indicates the proportion of permuted sex-specific counts that are equal to or larger than the mean observed sex-specific count (the mean count from all SUBsampled observed data), i.e. the overlap of permuted vs. observed distributions.
	
	collected_results = {}
	for i in range(1, min_n+1):
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


		
		collected_results[i] = {"obs_means" : [male_spec_obs, female_spec_obs], "obs_std" : [male_spec_obs_std, female_spec_obs_std], "p_val" : [p_male, p_female], "total_contigs_stats" : pres_abs_resampled_results[i][6:] }
	
#	print collected_results
	
	## write to a file:
	with open("permutation_results." + output_name + ".txt", "w") as OUTFILE:
		outlines = [["n_samples_per_sex", "male_specific_mean", "male_specific_std", "p_val", "female_specific_mean", "female_specific_std", "p_val", "total_contigs_m_mean", "total_contigs_m_std", "total_contigs_f_mean", "total_contigs_f_std"]]
		for i in range(1, min_n+1):
			outlines.append( [str(x) for x in [i, collected_results[i]["obs_means"][0], collected_results[i]["obs_std"][0], collected_results[i]["p_val"][0], collected_results[i]["obs_means"][1], collected_results[i]["obs_std"][1], collected_results[i]["p_val"][1], collected_results[i]["total_contigs_stats"][0], collected_results[i]["total_contigs_stats"][1], collected_results[i]["total_contigs_stats"][2], collected_results[i]["total_contigs_stats"][3] ] ] )
		
		
		outlines = [ "\t".join(x) for x in outlines[:] ]
		OUTFILE.write( "\n".join(outlines)  + "\n")
	
	

def count_item_occurence(lst):
    
    res = collections.defaultdict(lambda: 0)
    for v in lst:
        res[v] += 1
    return res
	
def MT_resampling(min_n, sexdict_idx, mapping_data, n_resampling, CPUs, permute_sexes):
	
	print "resampling . . ."	
	results = {}
	
	pool = mp.Pool(CPUs) #use all available cores, otherwise specify the number you want as an argument
	for i in range(n_resampling): # there is one worker for each resampling round
		results[i] = pool.apply_async(pres_abs_MT, args=(min_n, sexdict_idx, mapping_data, permute_sexes))
	pool.close()
	pool.join()


	# Get process results from the output queue
	#print output
	results1 = {}
	for i, result in results.items():
		results1[i] = result.get()

	return results1

####


#####


def pres_abs_MT(min_n, sexdict_idx, mapping_data, permute_sexes):
	
#	print permute_sexes
	
	# random.seed() is necessary to ensure that random has a different seed in 
	# each thread; otherwise each thread will return the same random.choice! It uses system time to make the seed.
	random.seed()
	
	# pres_abs_result returns number (count) of "sex-specific" contigs for each stringency level i
	pres_abs_result = []
	
	# a dictionary to hold the IDs of the reference contigs/loci/RADtags which are private to each sex
	# keys = confidence level i (pres-abs number of samples); values: [[male specific contig IDs],[female specific contig IDs]]
	contig_details = {} 
	
	shared_among_all_samples = 0

	if permute_sexes == False:
		for i in range(1, min_n+1):
#			print i
			jackn_males = random.sample(sexdict_idx["male"], i)
			jackn_females = random.sample(sexdict_idx["female"], i)			
#			print jackn_males, jackn_females
		
			pres_abs_data = get_pres_abs (mapping_data, jackn_males, jackn_females)

			# [male_specific, female_specific, male_total_loci, female_total_loci]
			pres_abs_result.append([0,0,0,0])
			
			contig_details[i] = [[],[]]
			
			# absence in exactly i samples is assured by the subsampling!
			# i.e. only i samples are scanned for pres-abs 	
			
			for loc in pres_abs_data.keys():
				if pres_abs_data[loc][0] == i:
					pres_abs_result[i-1][2] += 1 # counting overall presence: number of contigs present in all males (total number of contigs at this stringency)
					if pres_abs_data[loc][1] == 0:
						# if seen in exactly i males and absent in exactly i females				
						pres_abs_result[i-1][0] += 1
						contig_details[i][0].append(loc)
						
				if pres_abs_data[loc][1] == i:
					pres_abs_result[i-1][3] += 1 # counting overall presence: number of contigs present in all females (total number of contigs at this stringency)
					if pres_abs_data[loc][0] == 0:
						# if seen in exactly i females and absent in exactly i males
						pres_abs_result[i-1][1] += 1
						contig_details[i][1].append(loc)
				
				# counting overall presence: number of contigs present in all males and all females (total number of contigs at this stringency)			
							
	else:
		# the permutation option: contig_details is left empty since meaningless; pres_abs_result returns number of "sex-specific" contigs for each stringency level i
		all_samples = sexdict_idx["female"] + sexdict_idx["male"]
		for i in range(1, min_n+1):
#			print i
			jackn_males = random.sample(sexdict_idx["male"], i)
			jackn_females = random.sample(sexdict_idx["female"], i)			
#			print jackn_males, jackn_females			
			pres_abs_data = get_pres_abs (mapping_data, jackn_males, jackn_females)

			pres_abs_result.append([0,0])

			for loc in pres_abs_data.keys():
				if pres_abs_data[loc][0] == i:
					if pres_abs_data[loc][1] == 0:				
						pres_abs_result[i-1][0] += 1
				if pres_abs_data[loc][1] == i:
					if pres_abs_data[loc][0] == 0:
						pres_abs_result[i-1][1] += 1
	
	return [pres_abs_result, contig_details]
	
def get_pres_abs (mapping_data, males, females):
	
	# get the histogram of locus presence / absence: count for each sex
	# mapping_data is a dictionary; keys = contigs ; values = list of pres/abs (1/0) info for the samples; samples are represented by a fixed list index!
	# e.g. { '403848_L105': [0, 1, 0, 1] }
	
	pres_abs_data = {}
#	print mapping_data.keys()
	for contig, mappedlist in mapping_data.items():
		males_presence = sum( [ mappedlist[x] for x in males ] )	
		females_presence = sum( [ mappedlist[x] for x in females ] )
		pres_abs_data[contig] = [males_presence, females_presence]
	
	return pres_abs_data
	

def get_pres_abs_per_contig (bamfile, min_cov):
	
	bash_command = "samtools idxstats {0}".format(bamfile)
	print bash_command
	
	p = subprocess.Popen(bash_command, shell=True, stdout=subprocess.PIPE) #, stderr=subprocess.STDOUT)
	stdout, stderr = p.communicate()
	stdout_dict = {}
#	print stdout
	for line in stdout.split("\n"):
#		line = p.stdout.readline()
#		print line
		if line == '' and p.poll() != None:
			break
		elif "*" in line: # discard last line
			continue
		else:			
			fields = line.strip("\n").split("\t")
			# 
			tag = fields[0]
			dp = int(fields[2])
			if dp >= min_cov:
				stdout_dict[tag] = 1
			else:
				stdout_dict[tag] = 0
	return stdout_dict

		

######################## MAIN

args = get_commandline_arguments ()
	
args.bam_suffix = args.bam_suffix.strip()

check_congruence (args.sexlistfile, args.bam_dir,args.bam_suffix)

sexdict = read_sexlist (args.sexlistfile)

privacy_rarefaction_core(args.bam_dir, sexdict, int(args.n_resampling), int(args.CPUs), args.bam_suffix, float(args.min_support_to_report_loci), int(args.min_cov), args.o)

print "Done!"
