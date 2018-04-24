# privacy-rarefaction 
is a tool to scan population genomic datasets for sex-specific loci.

Privacy-rarefaction was developed to distinguish biological, population genetic polymorphisms from stochastic artefacts in presence-absence data. Common reduced-representation sequencing methods such as GBS, RAD-seq, ddRAD-seq, genome skimming, etc. are all partially stochastic in respect to which loci are sequenced in which individual, hence actual genomic presence-absence polymorphisms are confounded. When seeking loci private to one of two sets of individuals (e.g. females and males, resp. Y or W chromosomes), a practical dilemma arises. One may either compare few individuals per set and thus study many loci, but a high proportion of these will erroneously be classified as private (false positives). Alternatively, one may compare many individuals per set to increase the  confidence of classification, but this necessarily discards a large proportion of sequenced loci (possibly loosing all loci). By iterating from few to many individuals per set (stringency), privacy-rarefaction generates an objective and repeatable assessment of biological signal and stochasticity in presence-absence, including intuitive visualisation and significance tests. Privacy-rarefaction is model-free and uses permutations of the real data to generate specific null-hypotheses against which candidate sex-specific loci are evaluated. It helps to get the most out of your sequencing data without discarding any information. More about privacy-rarefaction can be found in the manuscript

Scharmann et al (). TBA

Please cite this manuscript if using privacy-rarefaction in your work.

## Requirements, Dependencies
Privacy-rarefaction is currently implemented as a command line script in python versions 2.6 and 2.7. It will run on desktop computers but it is advised to use a linux computer cluster for speed and simplicity. Mac and Windows operating systems are not tested. Required python modules:
- sys
- os
- argparse
- time
- numpy
- random
- subprocess
- multiprocessing

Privacy-rarefaction will call samtools idxstats (LINK), hence you need to have 'samtools' in the $PATH variable.

## Installation
No installation required. Download / copy the main script 'privacy-rarefaction.py' and the plotting function for R, or simply clone this repository.

## Inputs
- mapped sequencing reads in .BAM format, with bamfile index (.BAI)
- a sex list: a textfile, tab-delimited, one line per individual (=sample), 1st column sample name congruent with the .bam files except the suffix, 2nd column the sex: "1" for male and "2" for female.

### Comment:
Privacy-rarefaction is not limited to RAD-seq and similar data types, but it can be used with any genomic NGS reads that were mapped to any reference contigs. In the case of denovo assembled RAD-tags, the reference contigs are numerous and very short, similar to the read length. The many short, indipendent contigs largely ensure that truly sex-specific sequences (male-specific region of the Y; female-specific region of the W) are not eclipsed by phsyically linked but non-sex-specific sequences (pseudo-autosomal region, regions of sufficient XY homology resp. ZW homology). This may be an issue when larger contigs from whole genome assemblies are used in privacy-rarefaction: Contigs that contain both sex-specific AND non-sex-specific regions would not be detected as sex-specific because they map reads from both sexes. If you have a whole genome assembly I recommend to chop it into overlapping windows and map your population reads to these instead of the intact, large genome contigs. In that case you should retain mutliply (ambiguously) mapped reads in the .BAM files to reflect that the overlapping windows are partially redundant. Not that this strategy only makes sense when your reference genome was made from an individual of the heterogametic sex! 

RNA-seq reads are not adviseable for privacy-rarefaction because they indicate only sex-specific gene expression, but not necessarily sex-specific genome content. In contrast, denovo assembled transcriptome contigs from RNA-seq data by e.g. Trinity are very useful as reference contigs against which genomic reads can be mapped. Privacy-rarefaction was successfully used to identify the complete sequences of sex-specific genes by mapping ddRAD-seq reads against transcriptome contigs.


## Outputs
- quantitative results

- qualitative results


# Tutorial
Make sure that python 2.6 or 2.7 and samtools idxstats are available:
```
samtools idxstats
python --version 
```
Download or clone this repository to your machine. Then navigate to that directory on the commandline, and unpack the example data:
```
git clone https://github.com/mscharmann/privacy-rarefaction 
cd privacy-rarefaction
```
Now execute privacy-rarefaction to check that the relevant python modules are installed. This will also show you the command line arguments of the script:
```
python privacy-rarefaction.v2.2.py --help
```
If it complains about missing modules you need to install or load them first in your environment.
Else, you are ready to run the privacy rarefaction example, so unpack the example data:
```
tar -xzf example.tar.gz
```
In the directory 'example', you will find 20 small .BAM files. These are not real reads mapped to real contigs but for the purpose of this tutorial they will behave as-if, except running much faster. There is also a textfile 'sexlist.txt'
```
python privacy-rarefaction.v2.2.py --bam_dir ./example/ --bam_suffix .bam --sex_list ./example/sexlist.txt --CPUs 12 --o example_run

```
After finishing you can inspect the three output files:
- male_specific_candidates.example_run.txt
- female_specific_candidates.example_run.txt
- permutation_results.example_run.txt








## generating toy .BAM files for testing the script.

This python code was used to produce the toy .SAM files for 10 males and 10 females. Only 5,000 randomly selected loci out of total 10,000 loci are sequenced in each sample, and 1,000 loci only occur in the females (W-hemizygous):
```
import string
import random
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
	return ''.join(random.choice(chars) for _ in range(size))


len_loci = 10
n_loci_total = 10000
n_loci_private = 1000
n_loci_covered_per_sample = 5000
loci_dict = {}
for i in range(1,n_loci_total + 1):
	loci_dict["c_" + str(i)] = id_generator(len_loci, "ACTG")


l = ["c_" + str(i) for i in range(1,n_loci_total + 1)]
shared_loci = l[:(n_loci_total-n_loci_private)+1]
## i.e. the contig numbers 9001 - 10000 are female-specific!

genotypes_dict = {}
for i in range(1,11):
	id = "sample_" + str(i)
	seq_ids = random.sample(shared_loci, n_loci_covered_per_sample)
	genotypes_dict[id] = seq_ids


for i in range(11,21):
	id = "sample_" + str(i)
	seq_ids = random.sample(l, n_loci_covered_per_sample)
	genotypes_dict[id] = seq_ids


## now export this as mock-.SAM files:
for sample in genotypes_dict.keys():
	outlines = []
	readcount = 0
	for contig in l: ## ALL contigs MUST appear in the header, even if they didn't attract any reads!
		outlines.append("@SQ	SN:" + contig + "	LN:"+str(len_loci))
	for contig in genotypes_dict[sample]:
		readcount += 1
		QNAME = "r" + str(readcount)
		FLAG = "0"
		RNAME = contig
		POS = "1"
		MAPQ = "0"
		CIGAR = str(len_loci) + "M"
		MRNM_RNEXT = "="
		MPOS_PNEXT = "1"
		ISIZE_TLEN  = ""
		SEQ = loci_dict[contig]
		QUAL = "?"*len(SEQ)
		TAGs = ""
		outlines.append("\t".join([QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, MRNM_RNEXT, MPOS_PNEXT, ISIZE_TLEN, SEQ, QUAL, TAGs]))
	with open(sample + ".sam", "w") as OUTF:
		OUTF.write( "\n".join(outlines) + "\n")

```

After this, the .SAM files were converted to .BAM files and indexed using samtools:

```
for i in {1..20} ; do
	samtools view -Sb sample_${i}.sam > sample_${i}.bam
	samtools index sample_${i}.bam
	rm sample_${i}.sam
done

```


