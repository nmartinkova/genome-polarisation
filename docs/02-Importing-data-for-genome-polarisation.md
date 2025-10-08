

# Importing data for genome polarisation

*Natália Martínková and Stuart J. E. Baird*

## Conversion from VCF format
The package `diemr` provides a function `vcf2diem`, which converts diploid genotypes in a vcf file (optionally, the vcf file can be gzipped) or an `vcfR` object loaded into `R` to a file with format required by `diem`.  To start, load the `diemr` package or install it from CRAN if it is not yet available:


``` r
# Attempt to load a package, if the package was not available, install and load
if(!require("diemr", character.only = TRUE)){
    install.packages("diemr", dependencies = TRUE)
    library("diemr", character.only = TRUE)
}
```

Next, decide whether all markers in the vcf file should be polarised from one file, or whether the analysis will benefit from parallelisation. Parallelisation is available on UNIX-based platforms, and we advise to use it for large datasets. `diem` itself can also work in parallel on multiple data _compartments_. vcf files and _compartments_ can correspond to the same thing OR you can concatenate/split vcf2diem outputs to produce `diem` input _compartments_.


``` r
# Path to the vcf file
inputfile <- system.file("extdata", "myotis.vcf", package = "diemr")
# File name for the output
# If working directory does not have writing privileges, the filename must include a path to a suitable folder
outputfile <- "myotis"
# Convert the vcf file to two diem files
vcf2diem(SNP = inputfile, filename = outputfile, chunk = 2)
# Expecting to include 11 markers per diem file.
# If you expect more markers in the file, provide suitable chunk size.
# Done with chunk 1
```

The `vcf2diem` function processes calls for sites into the `diem` encoding for sites. The site order will be identical to that in the original vcf file (but see chapter [Sites not informative for genome polarisation]). 




## Principle of conversion to diem format

The `diemr` package uses a concise genome representation, differentiating homozygotes for the two most frequent alleles at each site, and their respective heterozygotes (mixed state of the two most common alleles). All other genotypes represent an unknown state. Specifically, the genotypes encoded as `0` represent homozygotes for one of the two most frequent alleles, 
`1` are heterozygous (mixed state) genotypes for the two most frequent alleles, 
`2` are homozygotes for the other allele, 
and `U` (symbol "_" is also allowed) is an unknown state that can be a missing genotype or a genotype that includes a third (or a fourth) allele.

### Invariant sites and indels

What logically follows from this description is that genome polarisation is meaningful only for variant sites (markers). Invariant sites will have **low support** and, though they will not change the `diem` outcome, they will slow down convergence and obscure the nature of change between taxa. During `diem` iterations their obscuring effect is removed by rescaling the hybrid indices onto $[0,1]$ 

Including invariant sites in genome polarisation is acceptable, and likely unavoidable. A proportion of truly invariant sites will be mis-diagnosed as variant due to sequencing errors.

At this time, genome polarisation uses variant markers where alleles differ in nucleotide substitutions. Insertions or deletions (indels) are not allowed to be classified among the most frequent alleles in the `vcf2diem` function. 



### Sites not informative for genome polarisation

Along side invariant sites, some variant sites are also uninformative for genome polarisation. Sites that include only missing data and heterozygous genotypes are unpolarisable. In general applications, one might wish to exclude sites that are only included in a `vcf` due to the presence of a single heterozygous individual. 

Uninformative sites will have support equal to 0, and their polarity will thus be equal to the null polarity. Including uninformative sites has similar consequences to including invariant sites. When a user wishes to retain such sites, we advice to then filter polarised markers with high support (or simply high diagnostic index, which is strongly correlated with support) for subsequent analyses and interpretation.


### Haploid data

Haploid data, including markers on sex chromosomes of the heterogametic sex, are usually shown as sequence alignments. To illustrate the principle of genome representation used in `diem`, let's have an 8bp-long alignment with five individuals. 


``` r
5 8
ind1   AACCTTGG
ind2   TACGATGG
ind3   TACC-TGT
ind4   CACGTTGG
ind5   AACCTTGT
```
We can see that the example alignment contains four variant markers at sites 1, 4, 5, and 8. Marker 5 contains a deletion in ind3, while other markers vary in different nucleotide substitutions. In marker 1, two individuals have an allele A, two individuals have an allele T and one has a C. This marker can be resolved in terms of two most frequent alleles A and T, where the individual ind4 has a third allele and is therefore considered as having an unknown genomic state. To convert the marker into `diem` format, we must make a decision on which allele will be encoded as 0. The choice is arbitrary. Good practice is to **keep a record of genotype encoding**. Let's decide that all alleles in ind1 will be encoded as 0. Marker 1 will then be `022_0`.

Contrary to DNA sequence alignments, diem format transposes the data. Rows in diem input files are markers and columns represent individual samples. The diem output for the example alignment will thus have four lines representing the four variant markers, and five columns representing the individuals. 


``` r
Marker1   022_0
Marker4   02020
Marker5   02_00
Marker8   00202
```

The diem format does not store information on marker identification or location. This metadata must be saved separately in a 'bed-like' file (with at minimum the scaffold and position for each site). `diem` uses a prefix `S` at line starts to ensure that all operating systems always read genotypes as strings without attempting to convert them into numbers. Encoding of markers 4 and 8 shown above is at risk of being interpreted as numeric. The final file in `diem` format will be:


``` r
S022_0
S02020
S02_00
S00202
```







### Diploid data

#### Heterozygotes in DNA sequence alignments

Diploid markers in DNA sequence alignment can be resolved in an analogous way to haploid markers. Heterozygotes in alignments can be identified according to their [IUPAC DNA ambiguity codes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2865858/table/T1/). For example, code W may represent an accepted heterozygote in Marker 1, where the most frequent alleles are A and T.


## Variant call format

Variant markers in genomic context will be most often identified using variant callers from sequence reads mapped onto a reference. Such data are stored or can be converted to variant call format (vcf). We implemented conversion of vcf files to diem format for the [vcf version 4.2](https://github.com/samtools/hts-specs/blob/master/VCFv4.2.pdf) in function `vcf2diem`. The function resolves indels and markers with more than two alleles to obtain biallelic genotypes for all markers included in the original vcf file.

1. Markers with the reference allele (column REF) containing insertions are set as unknown.
1. Markers with the alternative allele (column ALT) containing only indels are set as unknown.
1. Indels in the alternative alleles are set as unknown.
1. Overall allele counts in markers with more than two alleles representing substitutions are calculated and the two most frequent alleles are selected. Ties are (arbitrarily) resolved according to the allele order in the vcf file. 
1. Third and fourth most frequent substitutions are set as unknown. 

The conversion might result in some markers no longer being polymorphic. We advise checking how frequent invariant markers are after conversion.



``` r
# Import the first converted file for all individuals 
# and without changing marker polarity
myotis <- importPolarized("myotis-01.txt", 
                          changePolarity = rep(FALSE, 11), 
                          ChosenInds = 1:14)
# Check whether a marker includes heterozygotes 
# or that there is more than one genotype
apply(myotis, MARGIN = 2, 
      FUN = \(x) any(grepl("1", x)) | sum(levels(factor(x)) %in% c("0", "1", "2")) > 1) 
#    m1    m2    m3    m4    m5    m6    m7    m8    m9   m10   m11 
# FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
```

Starting from version 1.2.2, `vcf2diem` removes markers that include 

1. only missing data and heterozygotes, 
1. missing data, homozygotes for the most frequent allele and one heterozygote, 
1. missing data, homozygotes for one ALT allele and one heterozygote.

We include a list of the removed markers in a bed-like file ending with *-includedSites.txt* in the working directory. The file includes information from columns `CHROM`, `POS`, and `QUAL` for the respective markers, and records the (good practice) **record of genotype encoding**.

We include a list of the removed markers in a bed-like file ending with *-omittedSites.txt* in the working directory. The file includes information from columns `CHROM`, `POS`, and `QUAL` for the respective markers, and provides the reason, why the marker was omitted. **The file omittedLoci.txt must be checked when interpreting polarised markers** to ensure correct marker annotation.



## Polyploid data

The same principles apply to higher ploidy datasets. Preparing data in `diem` format for genome polarisation includes these steps:

1. Remove indels.
1. Identify two most frequent alleles for each site.
1. Encode homozygous and heterozygous (mixed) genotypes that include only the two most frequent alleles.

Note that `vcf2diem` only converts diploid datasets at this time, and should not be used to convert polyploid vcf files. 






