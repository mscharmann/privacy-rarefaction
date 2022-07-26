#!/usr/bin/python
# python 2.6 or 2.7 

# 
# Mathias Scharmann
# 2020-03-23


import os
import argparse
import random

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
	
	parser.add_argument("--kmc_dir", required=True, help="path of directory with jellyfish dump files", metavar="DIRECTORY")
	parser.add_argument("--sex_list", required=True,
        dest="sexlistfile", type=extant_file,
		help="name and path of the sex_list file: 1st column barcode/sample name separated by tab from second column indicating the sex; males = 1 and females = 2", metavar="FILE")
	parser.add_argument("--CPUs", required=True, help="number of CPUs to use in multiprocessing/parallel parts of script", metavar="INT")
	parser.add_argument("--o", required=True, help="name for the output files", metavar="STRING")
	parser.add_argument("--n_resampling", nargs='?', help="number of resampled datasets to be drawn for jacknifing over sample size & number of permutations for sex vs. sample ID, default = 1", metavar="INT", default = "1")
	parser.add_argument("--max_sps", nargs='?', help="maximum number of samples per sex (sps) to explore", metavar="INT", required=True)
	args = parser.parse_args()
	
	linebreak_check(args.sexlistfile)
	
	return args

#######

def check_congruence (sexlistfile, dump_folder):
	
	sexlistsamples = []
	with open(sexlistfile, "r") as INFILE:
		for line in INFILE:
			if len(line) > 1:
				fields = line.strip("\n").split("\t")
				sexlistsamples.append(fields[0])
	sexlistsamples = set(sexlistsamples)
	
	for sample in sexlistsamples:
		extant_file(dump_folder + sample )
	print "all samples in sex_list also have KMC database files in the kmc_dir, good to go!"
	
		


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


def privacy_rarefaction_resample(kmc_dir, sexdict, n_resampling, CPUs, output_name, max_sps):
	
	n_resampling = int(n_resampling)
		
	for n in range(4,max_sps+1):
		for i in range(n_resampling):
			print "samples per sex = ", n, " iteration = ", i
			outfname = "sps_" + str(n) + "_iter_" + str(i) 
			
			# do not run again if output already exists
			existing_files = set( ["sps_" + x.strip("_f_spec.kmc_suf").strip("_m_spec.kmc_suf") for x in os.listdir(".") if x.endswith("_spec.kmc_suf")] )
			if not outfname in existing_files:
			
				jackn_males = random.sample(sexdict["male"], n)
				jackn_females = random.sample(sexdict["female"], n)			
			
				run_KMC_tools (outfname,CPUs,jackn_males,jackn_females,kmc_dir)
				



def run_KMC_tools (outfname,CPUs,jackn_males,jackn_females,kmc_dir):
	
	chunksize = 16
	
	if len(jackn_males) <= 0.5*chunksize:
	
		kmc_params_lines_both = ["INPUT:"]
		mspec_cmd1 = "m_spec = ("
		mspec_cmd2 = ")-("
		fspec_cmd1 = "f_spec = ("
		fspec_cmd2 = ")-("
			
		f_list = []
		m_list = []
			
		cnt = 0
		for s in jackn_males:
			cnt += 1
			kmc_params_lines_both.append("s" + str(cnt) + " = " + kmc_dir + "/" + s)
			m_list.append("s" + str(cnt))
		for s in jackn_females:
			cnt += 1
			kmc_params_lines_both.append("s" + str(cnt) + " = " + kmc_dir + "/" + s)
			f_list.append("s" + str(cnt))

		kmc_params_lines_m = list(kmc_params_lines_both)
		kmc_params_lines_m.append("OUTPUT:")
		kmc_params_lines_m.append( outfname + "_m_spec = (" + "*".join(m_list) + ")-(" + "+".join(f_list) + ")"  )
			
		kmc_params_lines_f = list(kmc_params_lines_both)	
		kmc_params_lines_f.append("OUTPUT:")
		kmc_params_lines_f.append( outfname + "_f_spec = (" + "*".join(f_list) + ")-(" + "+".join(m_list) + ")"  )
	
		with open("kmc_params_file_mspec", "w") as O:
			O.write("\n".join(kmc_params_lines_m) + "\n")
		with open("kmc_params_file_fspec", "w") as O:
			O.write("\n".join(kmc_params_lines_f) + "\n")
			
		cmd1 = "kmc_tools -t" + str(CPUs) + " complex kmc_params_file_mspec"
		cmd2 = "kmc_tools -t" + str(CPUs) + " complex kmc_params_file_fspec"
			
		os.system(cmd1)
		os.system(cmd2)
	
	# to keep RAM usage acceptable, we need to form set unions and set intersections ITERATIVELY!
	else:
		# split samples in chunks of chunksize; 4 => up to 64 samples per sex.

		# male union
		iterative_union (jackn_males, kmc_dir, chunksize, CPUs, "level1")	
		level1_outputs = list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("level1_union")]))

		iterative_union (level1_outputs, ".", chunksize, CPUs, "level2")
		level2_outputs = list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("level2_union")]))
			
		iterative_union (level2_outputs, ".", chunksize, CPUs, "level3")
		level3_outputs = list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("level3_union")]))

		os.system("rm level2* level1*")
		os.system("mv level3_union_1.kmc_pre males_union.kmc_pre ; mv level3_union_1.kmc_suf males_union.kmc_suf ")
			
		# female union
		iterative_union (jackn_females, kmc_dir, chunksize, CPUs, "level1")	
		level1_outputs = list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("level1_union")]))

		iterative_union (level1_outputs, ".", chunksize, CPUs, "level2")
		level2_outputs = list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("level2_union")]))

		iterative_union (level2_outputs, ".", chunksize, CPUs, "level3")
		level3_outputs = list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("level3_union")]))

		os.system("rm level2* level1*")
		os.system("mv level3_union_1.kmc_pre females_union.kmc_pre ; mv level3_union_1.kmc_suf females_union.kmc_suf ")

		# male intersect
		iterative_intersection (jackn_males, kmc_dir, chunksize, CPUs, "level1")	
		level1_outputs = list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("level1_intersection")]))

		iterative_intersection (level1_outputs, ".", chunksize, CPUs, "level2")
		level2_outputs = list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("level2_intersection")]))

		iterative_intersection (level2_outputs, ".", chunksize, CPUs, "level3")
		level3_outputs = list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("level3_intersection")]))

		os.system("rm level2* level1*")
		os.system("mv level3_intersection_1.kmc_pre males_intersection.kmc_pre ; mv level3_intersection_1.kmc_suf males_intersection.kmc_suf ")

		# female intersect
		iterative_intersection (jackn_females, kmc_dir, chunksize, CPUs, "level1")	
		level1_outputs = list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("level1_intersection")]))

		iterative_intersection (level1_outputs, ".", chunksize, CPUs, "level2")
		level2_outputs = list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("level2_intersection")]))

		iterative_intersection (level2_outputs, ".", chunksize, CPUs, "level3")
		level3_outputs = list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("level3_intersection")]))

		os.system("rm level2* level1*")
		os.system("mv level3_intersection_1.kmc_pre females_intersection.kmc_pre ; mv level3_intersection_1.kmc_suf females_intersection.kmc_suf ")
		
		
		os.system("kmc_tools simple males_intersection females_union kmers_subtract " + outfname + "_m_spec")
		os.system("kmc_tools simple females_intersection males_union kmers_subtract " + outfname + "_f_spec")
		os.system("rm males_intersection* females_intersection* males_union* females_union*")

def iterative_intersection (names_in, kmc_dir, chunksize, CPUs, prefix):
		
	chunks = list_to_chunks(names_in, chunksize)
	chunk_count = 0
	for c in chunks:
		chunk_count += 1
		KMC_tools_intersect (CPUs,c,prefix + "_intersection_"+str(chunk_count) ,kmc_dir)


def iterative_union (names_in, kmc_dir, chunksize, CPUs, prefix):
		
	chunks = list_to_chunks(names_in, chunksize)
	chunk_count = 0
	for c in chunks:
		chunk_count += 1
		KMC_tools_union(CPUs,c,prefix + "_union_"+str(chunk_count) ,kmc_dir)
	
def KMC_tools_union (CPUs,samplenames,outputname,kmc_dir):
	
	s_list = []
	kmc_params_lines = ["INPUT:"]
	cnt = 0
	for s in samplenames:
		cnt += 1
		kmc_params_lines.append("s" + str(cnt) + " = " + kmc_dir + "/" + s)
		s_list.append("s" + str(cnt))
	kmc_params_lines.append("OUTPUT:")
	kmc_params_lines.append( outputname + " = " + "+".join(s_list) )
	
	with open("kmc_params_file_union", "w") as O:
		O.write("\n".join(kmc_params_lines) + "\n")
	
	cmd1 = "kmc_tools -t" + str(CPUs) + " complex kmc_params_file_union"
	cmd2 = "rm kmc_params_file_union"
	os.system(cmd1)
	os.system(cmd2)


def KMC_tools_intersect (CPUs,samplenames,outputname,kmc_dir):
	
	s_list = []
	kmc_params_lines = ["INPUT:"]
	cnt = 0
	for s in samplenames:
		cnt += 1
		kmc_params_lines.append("s" + str(cnt) + " = " + kmc_dir + "/" + s)
		s_list.append("s" + str(cnt))
	kmc_params_lines.append("OUTPUT:")
	kmc_params_lines.append( outputname + " = " + "*".join(s_list) )
	
	with open("kmc_params_file_intersect", "w") as O:
		O.write("\n".join(kmc_params_lines) + "\n")
	
	cmd1 = "kmc_tools -t" + str(CPUs) + " complex kmc_params_file_intersect"
	cmd2 = "rm kmc_params_file_intersect"
	os.system(cmd1)
	os.system(cmd2)


def list_to_chunks (lst, n):
	
	return [lst[i:i + n] for i in range(0, len(lst), n)]


def privacy_rarefaction_evaluate(max_sps):
	
	with open("privrar_kmers.report.txt", "w") as O:
  		O.write("\t".join(["subsample_size_per_group","group1_specific_avg","group1_specific_std","group2_specific_avg","group2_specific_std"])+"\n")
	
	
	datasets =  list(set([x.split(".")[0] for x in os.listdir(".") if x.startswith("sps_")]))
	
	import numpy
	
				# do not run again if output already exists
	existing_files = set( [x for x in os.listdir(".") if x.startswith("cumulative")] )
	
	for i in range(4,max_sps+1):
		mspec = []
		fspec = []
		for d in datasets:
			if int(d.split("_")[1]) == i:
				if "m_spec" in d:
					if not "cumulative_sps_" + str(i) + "_m_spec.txt" in existing_files:
						cmd = "kmc_tools -hp transform " + d + " dump /dev/stdout |" + """ awk '{print ">\\n"$1}' > """ + d + ".fasta"
						print cmd
						os.system(cmd)

				if "f_spec" in d:
					if not "cumulative_sps_" + str(i) + "_f_spec.txt" in existing_files:
						cmd = "kmc_tools -hp transform " + d + " dump /dev/stdout |" + """ awk '{print ">\\n"$1}' > """ + d + ".fasta"
						print cmd
						os.system(cmd)
				
				cnt = int( os.popen("kmc_tools -hp transform " + d + " dump /dev/stdout | wc -l").read() )
				if d.endswith("m_spec"):
					mspec.append(cnt)
				else:
					fspec.append(cnt)
		
		cmd = "cat sps_" + str(i) + "_iter_*_f_spec.fasta > cumulative_sps_" + str(i) + "_f_spec.fasta ; rm sps_*f_spec.fasta"
		print cmd
		os.system(cmd)
		cmd = "cat sps_" + str(i) + "_iter_*_m_spec.fasta > cumulative_sps_" + str(i) + "_m_spec.fasta ; rm sps_*m_spec.fasta"
		print cmd
		os.system(cmd)
		
		cmd = "kmc -m2 -sm -k25 -fa -ci1 -t1 cumulative_sps_" + str(i) + "_f_spec.fasta FFF ./ ; kmc_tools -hp transform FFF dump cumulative_sps_" + str(i) + "_f_spec.txt ; rm cumulative_sps_" + str(i) + "_f_spec.fasta FFF.*"
		print cmd
		os.system(cmd)

		cmd = "kmc -m2 -sm -k25 -fa -ci1 -t1 cumulative_sps_" + str(i) + "_m_spec.fasta MMM ./ ; kmc_tools -hp transform MMM dump cumulative_sps_" + str(i) + "_m_spec.txt ; rm cumulative_sps_" + str(i) + "_m_spec.fasta MMM.*"
		print cmd
		os.system(cmd)
		
		
		
		print i, numpy.mean(mspec), numpy.std(mspec), numpy.mean(fspec), numpy.std(fspec) 			
		with open("privrar_kmers.report.txt", "a") as O:
  			O.write("\t".join([str(x) for x in [i, numpy.mean(mspec), numpy.std(mspec), numpy.mean(fspec), numpy.std(fspec)] ] )+"\n")
		
## dump to fasta for a meta-kmer counting!?

		
######################## MAIN

args = get_commandline_arguments ()

samples = list(set([x.rstrip(".kmc_pre") for x in os.listdir(args.kmc_dir) if x.endswith(".kmc_pre")]))

sexdict = read_sexlist (args.sexlistfile)

# How many samples per ses (sps) shall be sampled at most?
min_n = min( len(sexdict["male"]), len(sexdict["female"]) )
if int(args.max_sps) >= min_n:
	max_sps = min_n
else:
	max_sps = int(args.max_sps)

privacy_rarefaction_resample(args.kmc_dir, sexdict, args.n_resampling, args.CPUs, args.o, max_sps)

privacy_rarefaction_evaluate(max_sps)

