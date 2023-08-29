# scikit-ribo-ext
Extended and improved scikit-ribo (https://github.com/schatzlab/scikit-ribo) for accurate estimation and robust modelling of translation dynamics at codon resolution

[![Documentation Status](https://readthedocs.org/projects/scikit-ribo/badge/?version=latest)](http://scikit-ribo.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Please see original [paper](https://doi.org/10.1016/j.cels.2017.12.007), GitHub [repo](https://github.com/schatzlab/scikit-ribo/tree/master) and [documentation](https://scikit-ribo.readthedocs.io/en/latest/?badge=latest) for details on scikit-ribo method, functionality, and usage.

## Introduction

Scikit-ribo has two major modules:

Ribosome A-site location prediction using random forest with recursive feature selection
Translation efficiency inference using a codon-level generalized linear model with a ridge penalty
A complete analysis with scikit-ribo has two major procedures:

* The data pre-processing step to prepare the ORFs, codons for a genome: scikit-ribo-build.py
* The actual model training and fitting: scikit-ribo-run.py

Inputs:
* The alignment of Riboseq reads (bam)
* Gene-level quantification of RNA-seq reads (from either Salmon or Kallisto)
* A gene annotation file (gtf)
* A reference genome for the model organism of interest (fasta)

Outputs:
* Translation efficiency estimates for the genes
* Translation elongation rate/ribosome dwell time for 61 sense codons
* Ribosome profile plots for each gene
* Diagnostic plots of the models


## New features and concepts
* Translation efficiency estimates and dwell times for subsets of genes (`--genesubset`)
* Translation elongation rate/ribosome dwell time for 61 sense codons ribosome for A, P, or E site (`--site`)
* Drastically improved speed for model training and fitting by implementing [polars](https://www.pola.rs/) to replace pandas
* Enabled analysis of large genomes by building individual indices for each chromosome

## How to run scikit-ribo for multicellular organisms

#### To run scikit-ribo for a large multicellular genome (e.g., human), you need to consider each chromosome separately. This avoids memory errors during the indexing of the multicellular genomes. Briefly, we build the scikit-index for each chr; then we merge the files from every single chr into a single file to build the scikit-ribo-index. 

### Data processing and indexing
**Critical:** Please run the script `scikit-ribo-build.py` that is located inside the `scikit-ribo-ext` directory of your cloned GitHub repository (for example, `/path/to/repo/scikit-ribo-ext/scikit-ribo-ext/scikit-ribo-build.py`). This is crucial as inside your environment you also have an installed version of the original pipeline which can be run by simply calling `scikit-ribo-run.py` without specifying the path to this version. This will run the original pipeline and will be missing the new features and improvements listed above.

To run scikit-ribo, it's important to use an unmasked version of the genome (only upper-case letters), without soft-masked regions (i.e., repetitive and low-complexity sequences in lower-case letters). 

1. Clone this repository and create environment from supplied yml file
```bash
git clone https://github.com/nedialkova-lab/scikit-ribo-ext.git
conda env create -f scikit-ribo-ext.yml
```

2. Activate environment
```bash
conda activate scikit-ribo-ext
```

3. remove the UTR and put the start of the first exon equal to the start of CDS and the stop of the CDS equal to the stop of the last exon. For this, run the script:
```bash
python3 -W ignore removes_3UTR_5UTR_from_GTF.py input.gtf
```

This script calls awk to format the file correctly. Specifically, it adjusts the first exon of every transcript to start only at the coordinates of the start codon and not before, taking into account the gene strand. It finally removes transcript/UTR annotations (in line with scikit-ribo in yeast). The formatted gtf file is the input appended with *corrected_annotation.gtf. Moreover, it deletes all the transcripts having more than one annotated start/stop codon, or not having a start/stop codon.

**Note:** The .txt output output of this script is important for filter_kallisto_file.py (which filters the Kallisto file and writes it correctly for scikit-ribo-build.py).

**Note:** The gtf output of this script is important to build the Kallisto index. It's important to consider also UTR regions to quantify mRNA levels.

4. Split the genome fasta file to avoid MemoryError. For this, create a directory "chr_fasta" and run the following command lines inside:
  * `pip install pyfaidx`
  * Split the genome: `faidx -x GRCh38.p13.genome.fa` (replace this fasta with your genome)
  * Delete all the non-canonical chromosomes and chrM.

5. Split the gtf file based on the chromosomes. For this, create a directory "chr_gtf" and run the following command line inside:
```bash
for i in /PATH/TO/chr_fasta/*.fa; do out=$(basename $i .fa); grep -w $out /PATH/TO/OUTPUT/FROM/3/*corrected_annotation.gtf > '*corrected_annotation.'$out'.gtf'; done
```
Note: when you grep 'chr1' you end up getting chr11, chr12 etc., adding `-w` solve the problem.

6. Filter the Kallisto file (i.e., abundance.tsv), keeping only 'protein_coding' transcripts:
```bash
python filter_kallisto_file.py /PATH/TO/OUTPUT/FROM/3/*.txt /PATH/TO/gene_id_gene_name_transcript_id_association.txt /PATH/TO/abundance.tsv
```
**Note:** the output file (i.e., abundance_protein_coding.tsv) should contain the gene ID in the first column, not the transcript ID. Otherwise, the file should be formatted accordingly, generating an association file from the gtf file (field 9): 
```bash
awk '$3 == "CDS" { print substr($10, 2, length($10)-3), substr($12, 2, length($12)-3), substr($16, 2, length($16)-3)}' gencode.v41.selected_principal_isoforms_filtered_CDS.gtf > gene_id_gene_name_transcript_id_association.txt
```
7. Build the scikit-ribo index for each chromosome separately, running the following command line in in the directory containing the two python scripts scikit-ribo-build.py and scikit-ribo-run.py: 
```bash
for i in /PATH/TO/chr_fasta/*.fa; do out=$(basename $i .fa); python3 scikit-ribo-build.py -g /PATH/TO/chr_gtf/*corrected_annotation.$out.gtf -f '/PATH/TO/chr_fasta_human/'$out'.fa' -p mRNA -t /PATH/TO/abundance_protein_coding.tsv -o 'scikit_index_'$out; done
```
**Note:*** The index fasta fai for the chr should be created starting from the fasta file to avoid problems running scikit-ribo. If the fai is older, remove it and create it again running scikit-ribo. 

**Note:** While running scikit-ribo-build.py, some chromosomes could generate memory errors because of their size (e.g., chr1 and chr2). This could create a memory error during the creation of the index. To avoid this, look at the chromosomes generating the error. Try to run for all these chromosomes scikit-ribo-build.py separately, avoiding other jobs on the server to use all the RAM. Otherwise, separate manually each gtf file with the chromosome producing memory errors into smaller files. Dividing it in two should be enough. In this context, it's important to include the full annotation for each gene included in each smaller gtf (i.e., not breaking the annotation). To not get errors, the resulting gtf files should be 2.5 Mb at most.

**Note:** the chrM is not important. Also if we get errors, we can skip it. This is not included into scikit-ribo.

8. Merge all the files with the same name associated with each chromosome. For this, run the following command line for each file:
```bash
cat /PATH/TO/scikit_index_chr*/file_name > file_name_merged
```
Example: `cat /PATH/TO/scikit_index_chr*/mRNA.3utr.fasta > mRNA_merged.3utr.fasta`

**Note:** it could be useful to create a directory 'scikit_index_merged' to contain all the merged files. This improves the organization of the working directory.

**Note:** before doing step 8, delete chromosome M (chrM). 

9. Delete all the repetitive headers in the *_merged.pos_ranges.txt file to avoid `ValueError: You are trying to merge on int64 and object columns. If you wish to proceed you should use pd.concat` while running scikit-ribo. For this, remove the lines from the source file itself, using the -i option with sed command:
```bash
sed -i -e 1b -e '/^#/d' mRNA_merged.pos_ranges.txt
```

10. Delete all the repetitive headers in the *_merged.nt_table.txt file to avoid `ValueError: invalid literal for int() with base 10: 'pos_ranges'` while running scikit-ribo. For this, remove the lines from the source file itself, using the -i option with sed command:
```bash
sed -i -e 1b -e '/chrom/d' mRNA_merged.nt_table.txt
```

11. Delete all the repetitive headers in the *_merged.codons.df file to avoid `ValueError: You are trying to merge on object and float64 columns. If you wish to proceed you should use pd.concat` while running scikit-ribo. For this, remove the lines from the source file itself, using the -i option with sed command:
```bash
sed -i -e 1b -e '/chrom/d' mRNA_merged.codons.df
```
12. Delete all the repetitive headers in the *_merged.codons.bed file to avoid errors while running scikit-ribo. For this, remove the lines from the source file itself, using the -i option with sed command:
```bash
sed -i -e 1b -e '/chrom/d' mRNA_merged.codons.bed
```

### Running model fitting
**Critical:** Please run the script `scikit-ribo-run.py` that is located inside the `scikit-ribo-ext` directory of your cloned GitHub repository (for example, `/PATH/TO/REPO/scikit-ribo-ext/scikit-ribo-ext/scikit-ribo-run.py`). This is crucial as inside your environment you also have an installed version of the original pipeline which can be run by simply calling `scikit-ribo-run.py` without specifying the path to this version. This will run the original pipeline and will be missing the new features and improvements listed above.
```bash
python3 scikit-ribo-run.py -i /PATH/TO/your_RPF_bam_file.bam -f ./PATH/TO/scikit_index_merged -p WT_mRNA_merged -o scikit_results
```

The length of reads used (`-s` and `-l`), the ribosome site used for dwell time estimation (`--site`) and running scikit-ribo on a subset of genes (`--genesubset`) can all be customized for individual needs.

