## naive_scoring.py

#!/usr/bin/python
# python 2.6 or 2.7 

# naive_scoring.py
# 
# Copyright (C) 2019, ETH Zurich, Mathias Scharmann
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




import os
import argparse
import subprocess

"""
example usage:
python /naive_scoring.py --bam_dir ../ --sex_list fufulist.txt --bam_suffix "-RG.bam "

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
	parser.add_argument("--o", required=True, help="name for the output files", metavar="STRING")
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

def naive_scoring_core(bam_dir, sexdict, bam_suffix, min_cov, output_name):
	
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
		
	# all M vs all F comparison: How many/which contigs never in F but at least 1 M, how many/which never in M but at least 1 F?
	male_specs = []
	female_specs = []
	for contig in mapping_data:
		males_pres = sum( [ mapping_data[contig][x] for x in sexdict_idx["male"] ] )
		females_pres = sum( [ mapping_data[contig][x] for x in sexdict_idx["female"] ] )
		if females_pres == 0 and males_pres != 0:
			male_specs.append(contig)
		elif males_pres == 0 and females_pres != 0:
			female_specs.append(contig)
	print "male-specific:		", len(male_specs)
	print "female-specific:	", len(female_specs)
	
	outlines = ["\t".join(["sex_specificity","contig_ID"])]
	for x in male_specs:
		outlines.append("male" + "\t" + x)
	for x in female_specs:
		outlines.append("female" + "\t" + x)

	with open(output_name + ".naive_scoring.txt", "w") as F:
		F.write("\n".join(outlines)+"\n")
			

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

naive_scoring_core(args.bam_dir, sexdict, args.bam_suffix, int(args.min_cov), args.o)

print "Done!"
