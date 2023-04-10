# Differential-Expression-Analysis
A tutorial guiding the use through differential expression analysis on RNA-Seq data

## Input Data:
- Raw Reads from sequencing, typically in fastq or fastq.gz format
- Reference genome in fasta format with GTF and GFF files
  - Ex: https://www.ncbi.nlm.nih.gov/assembly/GCF_000002985.6/

## Workflow:

### 1. Trim raw reads with Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic): 

```
mkdir trimmed-reads

trimmomatic PE -threads 24 $forward $reverse\ 
    trimmed-reads/$paired_forward trimmed-reads/$unpaired_forward\
    trimmed-reads/$paired_reverse trimmed-reads/$unpaired_reverse\
    ILLUMINACLIP:$ADAPTERS:2:30:10:8:true\
    LEADING:3 TRAILING:3\
    SLIDINGWINDOW:4:10 MINLEN:36
```
This will perform trimming on one set of paired end reads.
Replace `$forward` with the R1 reads and `$reverse` with the R2 reads for each sample.
Replace `$paired_forward`, `$unpaired_forward`, `$paired_reverse`, and `$unpaired_reverse` with the desired names.
`$ADAPTERS` must be a path to a fasta file with the sequencing adapters used.


### 2. Build index for reference genome with hisat2-build (https://daehwankimlab.github.io/hisat2/manual/):

```
hisat2-build <reference_in> <ht2_base>

hisat2-build GCF_000002985.6_WBcel235_genomic.fna GCF_000002985.6_WBcel235_genomic.fna.idx
```

This will create several index files for your reference.
```
GCF_000002985.6_WBcel235_genomic.fna.idx.1.ht2
GCF_000002985.6_WBcel235_genomic.fna.idx.2.ht2
...
GCF_000002985.6_WBcel235_genomic.fna.idx.8.ht2
```

### 3. Map paired and trimmed reads to reference genome using hisat2:
The max intron length can be determined using the script `intron-length.awk`. This script is called with the gff as follows:
```
./intron-length.awk TYPE=CDS GCF_000002985.6_WBcel235_genomic.gff
```
The script will return the following info to standard out:
`MINIMUM_INTRON_LENGTH MAXIMUM_INTRON_LENGTH MAXIMUM_SUM_INTRON_LENGTHS`

```
hisat2 [options]* -x <hisat2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <hit>]

hisat2 -x GCF_000002985.6_WBcel235_genomic.fna.idx -1 $paired_forward -2 $paired_reverse -S $sam_name -p 24 --max-intronlen MAXIMUM_INTRON_LENGTH
```
This will map one set of paired and trimmed reads to the reference genome.
Replace `$paired_forward` with the trimmed paired R1 reads and `$paired_reverse` with the trimmed paired R2 reads for each sample.
`$sam_name` refers to the output name of the alignment file, this file will be in the SAM format (https://samtools.github.io/hts-specs/SAMv1.pdf)

### 4. Counting reads in features with htseq-count:
```
htseq-count [options] <alignment_files> <gff_file>

htseq-count -f sam $sam_name GCF_000002985.6_WBcel235_genomic.gtf -t gene > $sam_name.features.txt
```
This will count the reads for one sample per gene.
Replace `$sam_name` with the SAM file for each sample.

### 5. Build count matrix using bash commands:
The output `sample.sam.features.txt` will be a table with two columns. The first column will be gene names, and the second will be respective counts for the sample. Other aligners can be used as long as the output is tab delimited with the gene names in the first column and the counts in the second column.

```
$ less -S sample.sam.features.txt

ATP6	2
CELE_2L52.1	22
CELE_2L52.2	6
CELE_2RSSE.1	4
CELE_2RSSE.3	15
...
__no_feature	13501027
__ambiguous	70624
__too_low_aQual	370657
__not_aligned	206645
__alignment_not_unique	1111979
```

To combine these individual counts files into a single matrix we can use the bash script `htseq_count_to_counts_matrix.sh`:
```
./htseq_count_to_counts_matrix.sh sample1.sam.features.txt sample2.sam.features.txt ... sampleX.sam.features.txt
```

The resulting file gene.counts.matrix can be used in DESeq2. It should have one column with gene names, and one column per sample with counts.

View this file by typing `less -S gene.counts.matrix`

```
	sample1	sample2	sample3	sample4	sample5	sample6	sample7	sample8
ATP6	2	0	1	4	5	4	1	3
CELE_2L52.1	1	0	0	3	0	5	0	2
CELE_2L52.2	0	0	0	0	0	0	0	0
CELE_2RSSE.1	0	1	3	1	0	3	9	5
CELE_2RSSE.3	0	0	0	0	0	0	0	0
...
__no_feature	13501027	15563129	15994076	14471195	15563966	16712179	13239092	14673742
__ambiguous	70624	68283	80613	62025	68120	77720	60419	61401
__too_low_aQual	370657	420911	475688	374362	461920	435703	349582	390528
__not_aligned	206645	235034	224522	211690	235723	227077	192822	205528
__alignment_not_unique	1111979	911511	1157079	775521	1041590	847218	709790	715698
```

### 6.1. Running DESeq2 (Python Implementation):
The script `unfiltered-deseq.py` can be used to run DESeq2 in Python. This requires the pydeseq2 library from https://github.com/owkin/PyDESeq2.
For visualization, bioinfokit is used (https://github.com/reneshbedre/bioinfokit)
Both libraries can be installed using pip.

Once all libraries are installed, test the command using `./unfiltered-deseq.py -h`. This will bring up the help menu.
```
usage: unfiltered-deseq.py [-h] [--alpha ALPHA] [--no_cooks] [--no_independent] [--plot] [--save] [--abnormal ABNORMAL [ABNORMAL ...]] [--normal NORMAL [NORMAL ...]]
                           [--out OUT]
                           matrix

Performs DESeq2 in python

positional arguments:
  matrix                Read counts matrix generated by kallisto

options:
  -h, --help            show this help message and exit
  --alpha ALPHA, -p ALPHA
                        P-value and adjusted p-value significance threshold
  --no_cooks            whether to not filter and refit cooks outliers
  --no_independent      whether to not independent filter
  --plot                whether to plot final data
  --save                whether to store intermediates in pkl format
  --abnormal ABNORMAL [ABNORMAL ...]
                        sample prefix for abnormal samples
  --normal NORMAL [NORMAL ...]
                        sample prefix for normal samples
  --out OUT, -o OUT     output prefix for all files
```
Different options can be changed depending on desired alpha value, and whether cooks refitting, independent filtering, and plotting data are desired.
In addition you must specify the abnormal and normal sample prefixes.

If we want to run DESeq2 with alpha of 0.05, no cooks refitting, samples 1-4 as normal and samples 5-8 as abnormal, we can run the following:
```
./unfiltered-deseq.py gene.counts.matrix --alpha 0.05 --no_cooks --plot --normal sample1 sample2 sample3 sample4 --abnormal sample5 sample6 sample7 sample8 --out deseq
```
The results from the DESeq2 run will be stored in multiple files. These correspond to the raw results, the non NaN results, the LFC shrunk results,
and the isolated upregulated and downregulated gene csvs. By default, these upregulated and downregulated genes have an adjusted p-value < alpha
and log2FoldChange > 1 or < -1 respectively.

### 6.2. Running DESeq2 (R Alternative):
The script `deseq.r` can be used to run DESeq2 in R. The R version will produce PCA plots and heatmaps for the data.
The multiple test p-value adjustment method in R is different from the Python method, producing slightly different (less strict) results.
To run the script, the libraries DESeq2, pheatmap, RColorBrewer, apeglm, ggplot2, and data.table must be installed.
```
Rscript ./deseq.r <counts csv> <affected sample prefix> <alpha value>
```
The script uses positional arguments including the path to the `gene.counts.matrix` file, the affected sample condition prefix, and the alpha value. Multiple files are produced including `raw_results.csv`, `deseq_results_non_nan.csv`, `shrunk_deseq_results_non_nan.csv`, `upregulated.csv`, `downregulated.csv`, `shrunk_upregulated.csv`, and `shrunk_downregulated.csv`. The upregulated and downregulated results are filtered on the provided alpha value. In addition, a heatmap with sample cladogram is saved to `heatmap.pdf` while a PCA plot is saved to `PCA.png`.
