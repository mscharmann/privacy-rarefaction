# privacy-rarefaction 
is a tool to scan population genomic datasets for sex-specific sequences (contigs of an assembly, or kmers if reference-free).

Privacy-rarefaction was developed to distinguish biological, population genetic polymorphisms from stochastic artefacts in presence-absence data. Common reduced-representation sequencing methods such as GBS, RAD-seq, ddRAD-seq, genome skimming, etc. are all partially stochastic in respect to which loci are sequenced in which individual, hence actual genomic presence-absence polymorphisms are confounded. When seeking loci private to one of two sets of individuals (e.g. females and males, resp. Y or W chromosomes), a practical dilemma arises. One may either compare few individuals per set and thus study many loci, but a high proportion of these will erroneously be classified as private (false positives). Alternatively, one may compare many individuals per set to increase the  confidence of classification, but this necessarily discards a large proportion of sequenced loci (possibly loosing all loci). By iterating from few to many individuals per set (stringency), privacy-rarefaction generates an objective and repeatable assessment of biological signal and stochasticity in presence-absence, including intuitive visualisation and significance tests. Privacy-rarefaction is model-free and uses permutations of the real data to generate specific null-hypotheses against which candidate sex-specific loci are evaluated. It helps to get the most out of your sequencing data without discarding any information. More about privacy-rarefaction can be found in the manuscript

Scharmann, M., Grafe, T. U., Metali, F. and Widmer, A. (2017). Sex-determination and sex chromosomes are shared across the radiation of dioecious Nepenthes pitcher plants. bioRxiv 240259.

Please cite this manuscript if using privacy-rarefaction in your work.

# NEW: KMER-based versions
for KMC version (```privacy-rarefaction.kmers.KMC.v1.py```), need additional dependency KMC. First run KMC count on each of the samples
```
samples=$( ls *.fq.gz | sed 's/.F.fq.gz//g' | sed 's/.R.fq.gz//g' | sort | uniq )
for sample in $samples ; do
echo $sample
echo ${sample}.F.fq.gz >input_file_names
kmc -m5 -sm -k25 -fq -ci2 -t12 @input_file_names ${sample} ./
rm input_file_names
done

```
Then run the python2.7 script, which will call KMC operations and produce a result table:
```
python2.7 privacy-rarefaction.kmers.KMC.v1.py --kmc_dir . --sex_list morph_phenotypes.txt --CPUs 12 --o FUFU --max_sps 15 --n_resampling 30
```


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

Privacy-rarefaction will call samtools idxstats (http://www.htslib.org/), hence you need to have 'samtools' in the $PATH variable.

## Installation
No installation required. Download / copy the main script 'privacy-rarefaction.py' and the plotting function for R, or simply clone this repository.

## Inputs
- mapped sequencing reads in .BAM format, with bamfile index (.BAI)
- a sex list: a textfile, tab-delimited, one line per individual (=sample), 1st column sample name congruent with the .bam files except the suffix, 2nd column the sex: "1" for male and "2" for female.

### Comment:
Privacy-rarefaction is not limited to RAD-seq and similar data types, but it can be used with any genomic NGS reads that were mapped to any reference contigs. In the case of denovo assembled RAD-tags, the reference contigs are numerous and very short, similar to the read length. The many short, indipendent contigs largely ensure that truly sex-specific sequences (male-specific region of the Y; female-specific region of the W) are not eclipsed by phsyically linked but non-sex-specific sequences (pseudo-autosomal region, regions of sufficient XY homology resp. ZW homology). This may be an issue when larger contigs from whole genome assemblies are used in privacy-rarefaction: Contigs that contain both sex-specific AND non-sex-specific regions would not be detected as sex-specific because they map reads from both sexes. If you have a whole genome assembly I recommend to chop it into overlapping windows and map your population reads to these instead of the intact, large genome contigs. In that case you should retain mutliply (ambiguously) mapped reads in the .BAM files to reflect that the overlapping windows are partially redundant. Note that this strategy only makes sense when your reference genome was made from an individual of the heterogametic sex! 

RNA-seq reads are not adviseable for privacy-rarefaction because they indicate only sex-specific gene expression, but not necessarily sex-specific genome content. In contrast, denovo assembled transcriptome contigs from RNA-seq data by e.g. Trinity are very useful as reference contigs against which genomic reads can be mapped. Privacy-rarefaction was successfully used to identify the complete sequences of sex-specific genes by mapping ddRAD-seq reads against transcriptome contigs.

### recent changes:
NEW in v 2.4: - parameter "--min_prop_contigs" to exclude samples with only few contigs covered; these can make the results worse / less stable. 
For example, if 9 out of 10 males are high-coverage but 1 has very few reads, it is extremely unlikely that any Y-specific contigs were sequenced in that low-coverage sample even though they may have been present in its genome. At 9 M versus 9 F, there is a small probability (10%) to sample only the 9 high-coverage males, but in 90% of resampled sets of males the problematic one will be included, causing the apparent absence of Y-contigs. Therefore, it can greatly increase the power of the method to exclude samples with unusually few contigs present (low coverage). 

## Outputs
- quantitative results: permutation_results.txt - to decide whether Y-hemizygous or W-hemizygous contigs exist

- qualitative results: female_specific_candidates.txt and male_specific_candidates.txt - to decide which contigs are the best candidates


# Tutorial
Make sure that python 2.6 or 2.7 and samtools idxstats are available:
```
samtools idxstats
python --version 
```
Download or clone this repository to your machine. Then navigate to that directory on the commandline:
```
git clone https://github.com/mscharmann/privacy-rarefaction 
cd privacy-rarefaction
```
Now execute privacy-rarefaction to check that the relevant python modules are installed. This will also show you the command line arguments of the script:
```
python privacy-rarefaction.v2.3.py --help
```
If it complains about missing modules you need to install or load them first in your environment.
Else, you are ready to run the privacy rarefaction example, so unpack the example data:
```
tar -xzf example.tar.gz
```
In the directory 'example', you will find 20 small .BAM files. These are not real reads mapped to real contigs but for the purpose of this tutorial they will behave as-if, except running much faster. There is also a textfile 'sexlist.txt'. You can now start the analysis like this:
```
python privacy-rarefaction.v2.3.py --bam_dir ./example/ --bam_suffix .bam --sex_list ./example/sexlist.txt --CPUs 12 --o example_run

```
After finishing you can inspect the three output files:
- permutation_results.example_run.txt

This is the quantitative result of privacy-rarefaction and is very intuitive when plotted with the included R script:
```
Rscript plot_privacy-rarefaction_curves.R permutation_results.example_run.txt
```
Open the resulting .pdf file and inspect the curves. The number of candidates for male- resp. female-specific candidates start out similar at low stringecies. In this example, only female-specific (W-hemizygous) contigs were simulated, so we know that any male-specific candidates are false-positives. As expected, the candidate counts start to diverge with increasing stringency, i.e. the male-specific candidates drop more steeply than the female-specific candidates. At stringency 4-5 the standard deviations do not overlap anymore. This indicates a significant quantitative difference in sex-specific candidates between the two mutually exclusive alternatives Y-hemizygous and W-hemizygous (there is also a p-value for this in the 'permutation_results.example_run.txt' file). However, there may still be some chance of false-positives among the female-specific candidates, as false-positive male-specific candidates are not yet eliminated at this stringency. But at stringecy 10, the count for male-specific candidates has dropped to zero (no data point on log-sale): here, the biologically plausible expectation of only one sex having sex-specific loci has been met, as well as the known truth (only female-specific loci were simulated). You have successfully distilled true presence-absence from noisy, confounded data.


The two further output files are the qualitative result of privacy-rarefaction:
- female_specific_candidates.example_run.txt
- male_specific_candidates.example_run.txt

They list female- resp. male-specific candidate loci identified by privacy-rarefaction, separately for each level of stringency (N males vs N females), and with the proportion of resampled male-female sets that the contig was sex-specific in ("bootstrap support"). In the toy example, only the contigs 9,000-10,000 were simulated as female-specific. You will notice that with adeaquate number of resampling rounds (>= 200), no false-positives will appear at the bottom of the list in 'female_specific_candidates.example_run.txt'. However, male-specific candidates (which are entirely false-positives in this example) never make it beyond stringency 3 or 4 and hence do not appear in 'male_specific_candidates.example_run.txt'. 

The quality, or confidence that a contig is actually sex-specific, increases with the stringency and the bootstrap support. Hence, if you were to design PCR primers for validation of these contigs, you should start with those at the very bottom of the lists first. These are the contigs that were most consistently sequenced from one of the sexes only.








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



# privacy-rarefaction with kmers
The resampling and bootstrapping idea behind privacy-rarefaction is NOT limited to reduced-representation data, but whole-genome data from a population can also be used. We provide a version here that scans for sex-specific kmers, as extracted by jellyfish from the read files.

Here is an example workflow:

## 1. extract kmers

Extract the kmers from each individual with jellyfish ( https://github.com/gmarcais/Jellyfish )

```
for sample in $( cut -f1 sex_list.txt ) ; do
echo $sample
ls ${sample}.*.fq.gz | xargs -n 1 echo gunzip -c > generators
jellyfish count -g generators -G 1 -m 93 -s 100M -t 10 -C -o ${sample}.mer_counts.jf
jellyfish dump ${sample}.mer_counts.jf > ${sample}.93mer_counts.jf.dump
done
```

## 2. run privacy-rarefaction

```
python privacy-rarefaction.kmers.v3.py --dump_dir ./ --sex_list sex_list.txt --CPUs 12 --o 93mers.mincount2.minshared6.minstringency8 --n_resampling 200 --min_support_to_report_kmers 0.2 --dump_suffix .93mer_counts.jf.dump --min_count 1 --min_shared 6 --min_stringency 7
```

