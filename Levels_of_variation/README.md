# Scripts for computing levels and patterns of genetic diversity and divergence

These scripts process VCFs to estimate levels and patterns of genetica variation (e.g. Watterson's theta, Pi and Tajima's D), as well as genetic divergence from allele frequencies (e.g. FST).
These scripts are essentially unchanged compared to this repository (although script names have been simplied here): https://github.com/andreaswallberg/Ecological-Genomics-Northern-Krill

In these examples, the "groups" specified in text files are:

```
Aut_plankt.txt = Autumn planktivorous
Spr_plankt.txt = Spring planktivorous
Nor_pisci.txt = Northern piscivorous (Slåttersill)
Sou_pisci.txt = Southern piscivorous
```

## vcf2diversity.pl

This script estimates Watterson's theta (population mutation rate), Pi (nucleotide diversity) and Tajima's D.

**Reuse of code from other projects:**

The code for Tajima's D originates from BioPerl and was originally written by Jason Stajich

https://metacpan.org/release/CJFIELDS/BioPerl-1.6.924/source/Bio/PopGen/Statistics.pm

The implementation here differs by reusing some already calculated variables.

### Usage example:

```
VCF=mydata.vcf.gz

./vcf2diversity.pl \
  --vcf $VCF \
  --output $VCF.diversity \
  --group \
    Aut_plankt=Aut_plankt.txt \
    Sou_pisci=Sou_pisci.txt \
    Nor_pisci=Nor_pisci.txt \
    Spr_plankt=Spr_plankt.txt \
  --window 100000000 \
  --coverage Ch_v2.0.2.fasta.mask.fasta \
  --region \
    1=any \
  --verbose
```

Additional options and other aspects of usage:

- The "--vcf" argument can take more than one VCF file.
- Samples should not be included in more than one group
- The "--window" argument can take any number of window sizes (at the cost of performance)
- If a "-genes" argument is provided together with a GFF/GTF file, the script will compute the average levels of diversity at different distances upstream/downstream of genes.
- The "groups" argument specifies the label of the group/population and points to a simple csv/txt file with the names of the samples for that group (one sample name per line).

Here, each group is defined by a simple list of individuals provided in separate files. A genome mask is given in fasta format that specifies accessible sites (e.g. 0=inaccessible; 1=accessible) and levels of variation are computed for the accessible sites ("any") in windows of an arbitrary size (100,000,000bp, i.e. one window per chromosome).

## vcf2allele_counts.pl

This script estimates SNP allele frequencies for populations by grouping together samples present in a VCF file.

The allele frequencies can be based on directly on the called genotypes ("GT") or on the genotype likelihoods ("GL") or on the genotype probablities ("GP").

The output is a tabular file that specifies the allele counts for every population at every SNP site. It generates two output files per "method" (GT, GL or GP), one containing the allele counts and the other containing the allele frequencies. The GL and GP methods incorporate some level of uncertainty in the genotyping.

The main advantage of this format is that it represents a condensed version of the variation data that is faster to parse downstream of this point compared to repeatedly parsing the VCF.

### Usage example:

```
VCF=mydata.vcf.gz

./vcf2allele_counts.pl \
  --input $VCF \
  --output ${VCF}.alleles \
  --gt \
  --groups \
    Aut_plankt=Aut_plankt.txt \
    Sou_pisci=Sou_pisci.txt \
    Nor_pisci=Nor_pisci.txt \
    Spr_plankt=Spr_plankt.txt
```

The "groups" argument specifies the label of the group/population and points to a simple csv/txt file with the names of the samples for that group (one sample name per line).

In depth description of the output format can be found here: https://github.com/andreaswallberg/Ecological-Genomics-Northern-Krill/tree/main/Population_genomics-SNP_processing

## allele_counts2fst.pl

This script computes FST per SNP (Weir-Cockerham) and/or per window (Reynolds). It sets FST values for SNPs <0 to 0.

### Usage example:

```
VCF=mydata.vcf.gz

./allele_counts2fst.pl \
  --group \
    Aut_plankt=Aut_plankt.txt \
    Sou_pisci=Sou_pisci.txt \
    Nor_pisci=Nor_pisci.txt \
    Spr_plankt=Spr_plankt.txt \
  --input $VCF.alleles.allele_counts.GT.csv \
  --window 10000
```

Per-SNP and Per-window output formats are described in detail here: https://github.com/andreaswallberg/Ecological-Genomics-Northern-Krill/tree/main/Population_genomics-Divergence

## fstwindows2counts.pl

A companion script to allele_counts2fst.pl. It scans the window-based FST output and counts the number of times a group has the minimal FST distance to another group.

```
VCF=mydata.vcf.gz

./fstwindows2counts.pl $VCF.alleles.allele_counts.GT.csv.fst.reynolds.windows.10000bp.csv \
  1> infile \
  2> $VCF.alleles.allele_counts.GT.csv.fst.reynolds.windows.10000bp.csv.windows.tsv
```
The output includes:

- An overall pairwise FST distance matrix between all groups in PHYLIP format ("infile").
- A list of all windows included in the analysis ("$VCF.alleles.allele_counts.GT.csv.fst.reynolds.windows.10000bp.csv.windows.tsv"), as windows with missing data are excluded.
- For each group, a "closest_group" tsv is generated with the counts for each chromosome / scaffold.
