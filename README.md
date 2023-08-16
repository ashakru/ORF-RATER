# ORF-RATER
ORF-RATER (Open Reading Frame - Regression Algorithm for Translational Evaluation of Ribosome-protected footprints) comprises a series of scripts for coding sequence annotation based on ribosome profiling data.

The software was created at [Jonathan Weissman's lab at UCSF](http://weissmanlab.ucsf.edu/) and is described in [Fields, Rodriguez, et al., "A regression-based analysis of ribosome-profiling data reveals a conserved complexity to mammalian translation", *Molecular Cell* **60**, 816-827 (2015).](http://dx.doi.org/10.1016/j.molcel.2015.11.013)

Usage information can be found in the Detailed Protocol included in the paper's supplemental materials, or by running each script with the --help/-h flag.

Required packages include [numpy](http://www.numpy.org), [scipy](http://www.scipy.org), [pysam](https://github.com/pysam-developers/pysam), [biopython](http://www.biopython.org), [pandas](http://pandas.pydata.org/), [tables](http://www.pytables.org/), [scikit-learn](http://scikit-learn.org/), [pybedtools](https://pythonhosted.org/pybedtools/), and [plastid](https://pypi.python.org/pypi/plastid), all of which are available through [PyPI](https://pypi.python.org/pypi).

Some features require the [multiisotonic](https://github.com/alexfields/multiisotonic) package, which must be downloaded manually. Multiisotonic additionally requires [python-igraph](https://github.com/igraph/python-igraph).

Transcripts must be presented in UCSC's [BED12 format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). The most reliable method I've found to convert from GTF to BED12 involves first converting to [genePred format](https://genome.ucsc.edu/FAQ/FAQformat.html#format9), making use of UCSC's "gtfToGenePred" and "genePredToBed" scripts, which are available [here](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/). The full command is `gtfToGenePred INPUT_GTFFILE.gtf stdout | genePredToBed stdin OUTPUT_BEDFILE.bed`. Similarly, a BED file can be converted to a GTF using the command `bedToGenePred INPUT_BEDFILE.bed stdout | genePredToGtf file stdin OUTPUT_GTFFILE.gtf`.

Contact Alex Fields for further information or assistance.

---
### Detailed protocol
Addendum by @ashakru: pasted from supplementary material [Fields et al. 2015](http://dx.doi.org/10.1016/j.molcel.2015.11.013)

This protocol is intended to provide additional details to researchers intending to apply the ORF-RATER pipeline to a new dataset. Experimental aspects of ribosome profiling are beyond the scope of this protocol; instead, please refer to Ingolia et al. (2012). Nat. Protoc. 7, 1534-1550. Here, we assume that the reader has already acquired ribosome profiling sequencing reads from the sample of interest. What follows are step-by-step instructions for applying the ORF-RATER algorithm to those reads to identify translated ORFs.

In this protocol, familiarity with the UNIX command line and standard bioinformatics file formats is assumed. Software has not been tested in non-UNIX environment. Python must be installed, along with required packages numpy, scipy, pysam, biopython, pandas, scikit-learn, and plastid. Certain features require pybedtools.

#### Starting files
1.	The appropriate genome sequence, formatted as a FASTA file. For common research organisms, these are available from many online resources such as the UCSC Genome Browser (genome.ucsc.edu) or Ensembl (ensembl.org).
2.	The appropriate transcriptome annotation. Canonical transcriptomes are available from online sources. In some cases, it may be preferable to generate a sample-specific transcriptome directly from RNA sequencing data; many tools are available for this purpose. ORF-RATER requires that transcriptome annotations be structured in 12-column BED format; converters from other formats can be found online. Each transcript in the transcriptome should be assigned a unique identifier in its third column.
3.	CDS annotations. If a canonical transcriptome is being used, it will likely include these already, in which case no additional files are needed; otherwise, canonical CDSs can be downloaded from online sources. These should also be formatted as a 12-column BED file (or multiple such files), with the CDS extending from “thickstart” to “thickend”.
4.	Ribosome profiling reads aligned to the selected transcriptome. If multiple types of ribosome profiling data were obtained (e.g. following treatment with different translation inhibitors), separate alignments should be generated separately for each type. Multimapping reads should be reported in all (or as many as possible) positions to which they can align. Files should be structured in BAM format with a corresponding BAM index file. Each dataset may be spread across multiple files (no concatenation necessary).
5.	Optional: a lookup table assigning each transcript ID to a gene name.
6.	Optional: list of transcript IDs annotated as pseudogenes.

#### Procedure
For simplicity, it is recommended that the following scripts be executed in an empty folder. Input and output filenames can be set manually for all scripts, but it is much simpler to use the default filenames so that each script automatically finds the files generated by the others.
1.	Run the `prune_transcripts.py` script on the transcriptome. This will remove from further consideration any poorly behaved transcripts, such as those harboring too few reads. The aligned reads provided to this script should not be those treated with an initiation inhibitor; typically the most deeply sequenced CHX or no-drug dataset is appropriate. The values of `MINLEN` and `MAXLEN` should be set to encompass only the peak of the ribosome profiling read distribution (typically around 29–30 nt); setting too low of a value will overestimate the influence of multimapping and increase execution time. Reads of other lengths (and from other datasets if available) will be used in later scripts.
2.	Run the `make_tfams.py` script on the pruned transcriptome generated in the previous step. This script identifies transcripts that overlap on the same chromosome and strand, and assembles them into transcript families, which are the core groupings for later analyses. If available, a lookup file can be used to assign relevant names to each family of transcripts.
3.	Run `find_orfs_and_types.py` to identify ORFs in each transcript family, to determine which correspond to annotated CDSs, and to demark the ORF type (e.g. upstream, extension) for all other ORFs. At this step, the user can select which codons to consider as possible sites of translation initiation (e.g. ATG or NTG). By default, annotations in the provided transcriptome will be included; if coding sequences are not (or are incorrectly) specified in the source transcriptome, the --ignoreannotations flag should be set. Either way, additional CDSs may be provided using --extracdsbeds.
4.	For each dataset acquired with a different drug treatment:
a.	Run `psite_trimmed.py` to identify appropriate offsets for each read length in the dataset. If multiple datasets are to be analyzed, `SUBDIR` should be set so as to distinguish them (e.g. “CHX”, “HARR”). At this step, minimum and maximum read lengths are set for each dataset. (The user can delete lines in the output `offsets.txt` file corresponding to unwanted read lengths after the fact, if so desired.) It is recommended that offsets.txt be inspected to ensure that the offsets generated for each read length are reasonable (typically around 12 nt). It is also recommended that the optional TALLYFILE be generated and inspected to ensure that for each read length, a peak of density is found at the chosen offset.
b.	Run `regress_orfs.py` to perform the least-squares regressions for each ORF. The SUBDIR parameter should match the setting in the previous step. For initiation inhibitor-treated datasets, the `--startonly` flag should be specified. The metagene produced as part of this script is likely to be of interest for plotting.
5.	Run `rate_regression_output.py` to compile the regression outputs into a combined rating for each ORF.
6.	Run `make_orf_bed.py` to produce a BED file with the rated ORFs.
7.	If so desired, run `quantify_orfs.py` to quantify the expression of each ORF across a set of BAM files. BAM files collected under different drug treatments should be processed independently, as quantification depends on the relative degree of phasing in each.

#### Command line notes
Detailed information about the scripts and their arguments can be found by running each with the --help/-h flag. Sensible default values are provided for most arguments. Some scripts print warnings to stderr; all except `make_orf_bed.py` can be set to log progress to stdout through the `--verbose/-v` flag. Many can be parallelized using `--numproc/-p`. Here are sample command lines for an experiment with two experimental conditions (Harr and no-drug), run on 12 processors and with logging concatenated in log.out:

```
prune_transcripts.py --inbed source_transcriptome.bed --summarytable tid_removal_summary.txt --pseudogenes /refpath/pseudo_list.txt -p 12 -v /refpath/genome.fa /bampath/nodrug1.bam /bampath/nodrug2.bam > logfile.out
make_tfams.py -v -g tids_to_genenames.txt >> logfile.out
find_orfs_and_types.py /refpath/genome.fa --codons NTG --extracdsbeds /refpath/more_cds.bed /refpath/even_more_cds.bed –p 12 –v >> logfile.out
psite_trimmed.py /bampath/harr.bam --minrdlen 27 --maxrdlen 33 --subdir Harr --tallyfile tallies.txt --cdsbed source_transcriptome.bed -p 12 -v >> logfile.out
regress_orfs.py /bampath/harr.bam --subdir Harr --startonly -p 12 -v >> logfile.out
psite_trimmed.py /bampath/nodrug1.bam /bampath/nodrug2.bam --minrdlen 27 --maxrdlen 34 --subdir ND --tallyfile tallies.txt --cdsbed source_transcriptome.bed -p 12 -v >> logfile.out
regress_orfs.py /bampath/nodrug1.bam /bampath/nodrug2.bam --subdir ND --restrictbystarts Harr -p 12 -v >> logfile.out
rate_regression_output.py Harr ND -p 12 -v >> logfile.out
make_orf_bed.py
quantify_orfs.py /bampath/nodrug1.bam /bampath/nodrug2.bam --subdir ND –p 12 –v >> logfile.out
```

Although most of these operations must be performed in the specified order, psite_trimmed.py has no earlier dependencies, and can be run at any time prior to regress_orfs.py. The final annotation output (ratedorfs.bed by default) can be visualized in genome browsers such as IGV (broadinstitute.org/igv). The quantification output is saved as a pandas dataframe (quant.h5 by default) that can be read using the pandas.read_hdf() function.
