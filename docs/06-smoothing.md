# Smoothing polarised genotypes

Genome polarisation analyses are typically run on SNVs that are unevenly distributed along chromosomes. When genome coordinates are available, diagnostic markers can be mapped to their chromosomal positions, allowing quantitative analyses along the genome. Mapping also enables smoothing of the polarised signal, reducing the effect of uneven marker density and standing genetic variaton on the visualisation of ancestry tracts and hybrid index variation.

## Mapping diagnostic markers

The chromosomal coordinates of diagnostic markers are stored in the output file *-includedSites.txt*, created during VCF conversion by the function `vcf2diem`. This file lists all sites retained for polarisation and includes columns CHROM, POS, QUAL from the VCF file and allele0 and allele2 indicating alleles represented by the respective numbers when homozygous.

The positional information can be extracted directly from this file. The function `rank2map` links the ranked order of sites in the polarised genotype matrix to their physical coordinates.


``` r
# read the coordinates of included sites
inc <- readIncludedSites("my-sample-includedSites.txt")

# generate mapping of ranks to genome positions
map <- rank2map(inc$CHROM, inc$POS)

head(map)
```

Each entry corresponds to one marker in the polarised genotype matrix. For genomes composed of multiple scaffolds or chromosomes, smoothing should be performed separately for each element, and `rank2map` ensures that the windows are set up with respect to `CHROM` values.


## Smoothing across mapped distances

The `smoothPolarizedGenotypes` function applies kernel smoothing to the polarised genotype matrix, using physical distances between markers rather than their rank order. This approach preserves the broad genomic signal while suppressing local irregularities caused by uneven SNV spacing.


``` r
# smooth the polarised genotypes using a 250 kb window
genSmooth <- smoothPolarizedGenotypes("testdata-001.txt",
	includedSites = "my-sample-includedSites.txt",
	window = 2.5e5)
```

The smoothing window is expressed in base pairs. Larger windows produce a more continuous signal but can merge distinct ancestry blocks, while smaller windows follow local fluctuations more closely. The optimal setting depends on marker density and linkage along the genome.

> Smoothing is recommended only for mapped data. Applying it to unordered or scaffold-level data can distort the genomic signal.

## Visual comparison of raw and smoothed results

Smoothing does not change the overall polarity of sites but modifies their local continuity. The following example compares raw and smoothed genotypes.



The smoothed profile follows the same overall trajectory as the raw data but removes short-range noise caused by individual markers.
The resulting curve highlights broad transitions in ancestry and provides a clearer estimate of barrier width along chromosomes.

Interpretation and downstream use

Mapping and smoothing improve interpretability of diem output in three ways:

they stabilise HI profiles by reducing local fluctuations,

they allow comparison of barrier strength across chromosomes in consistent physical units, and

they enable integration with other genomic features, such as recombination rate or gene density.

When reporting smoothed results, specify the window size and clarify whether smoothing was applied per chromosome or across the whole genome.
Excessive smoothing can obscure narrow introgression tracts, while insufficient smoothing may fail to suppress sequencing artefacts.
