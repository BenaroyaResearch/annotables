annotables
================

[![DOI](https://zenodo.org/badge/3882/stephenturner/annotables.svg)](https://zenodo.org/badge/latestdoi/3882/stephenturner/annotables)

Provides tables for converting and annotating Ensembl Gene IDs.

## Installation

This is an [R](https://www.r-project.org) package.

### [Bioconductor](https://bioconductor.org) method

``` r
source("https://bioconductor.org/biocLite.R")
biocLite("stephenturner/annotables")
```

### [devtools](https://cran.r-project.org/package=devtools) method

``` r
install.packages("devtools")
devtools::install_github("stephenturner/annotables")
```

## Rationale

Many bioinformatics tasks require converting gene identifiers from one
convention to another, or annotating gene identifiers with gene symbol,
description, position, etc. Sure,
[biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
does this for you, but I got tired of remembering biomaRt syntax and
hammering Ensembl’s servers every time I needed to do this.

This package has basic annotation information from **Ensembl Genes 98**
for:

  - Human build 38 (`grch38`)
  - Human build 37 (`grch37`)
  - Mouse (`grcm38`)
  - Rat (`rnor6`)
  - Chicken (`galgal5`)
  - Worm (`wbcel235`)
  - Fly (`bdgp6`)
  - Macaque (`mmul801`)

Where each table contains:

  - `ensgene`: Ensembl gene ID
  - `entrez`: Entrez gene ID
  - `symbol`: Gene symbol
  - `chr`: Chromosome
  - `start`: Start
  - `end`: End
  - `strand`: Strand
  - `biotype`: Protein coding, pseudogene, mitochondrial tRNA, etc.
  - `description`: Full gene name/description

Additionally, there are `tx2gene` tables that link Ensembl gene IDs to
Ensembl transcript IDs.

## Usage

``` r
library(annotables)
```

Look at the human genes table (note the description column gets cut off
because the table becomes too wide to print nicely):

``` r
grch38
```

    ## # A tibble: 67,415 x 9
    ##    ensgene  entrez symbol  chr     start    end strand biotype description      
    ##    <chr>     <int> <chr>   <chr>   <int>  <int>  <int> <chr>   <chr>            
    ##  1 ENSG000…   7105 TSPAN6  X      1.01e8 1.01e8     -1 protei… tetraspanin 6 [S…
    ##  2 ENSG000…  64102 TNMD    X      1.01e8 1.01e8      1 protei… tenomodulin [Sou…
    ##  3 ENSG000…   8813 DPM1    20     5.09e7 5.10e7     -1 protei… dolichyl-phospha…
    ##  4 ENSG000…  57147 SCYL3   1      1.70e8 1.70e8     -1 protei… SCY1 like pseudo…
    ##  5 ENSG000…  55732 C1orf1… 1      1.70e8 1.70e8      1 protei… chromosome 1 ope…
    ##  6 ENSG000…   2268 FGR     1      2.76e7 2.76e7     -1 protei… FGR proto-oncoge…
    ##  7 ENSG000…   3075 CFH     1      1.97e8 1.97e8      1 protei… complement facto…
    ##  8 ENSG000…   2519 FUCA2   6      1.43e8 1.44e8     -1 protei… alpha-L-fucosida…
    ##  9 ENSG000…   2729 GCLC    6      5.35e7 5.36e7     -1 protei… glutamate-cystei…
    ## 10 ENSG000…   4800 NFYA    6      4.11e7 4.11e7      1 protei… nuclear transcri…
    ## # … with 67,405 more rows

Look at the human genes-to-transcripts table:

``` r
grch38_tx2gene
```

    ## # A tibble: 249,158 x 2
    ##    enstxp          ensgene        
    ##    <chr>           <chr>          
    ##  1 ENST00000373020 ENSG00000000003
    ##  2 ENST00000496771 ENSG00000000003
    ##  3 ENST00000494424 ENSG00000000003
    ##  4 ENST00000612152 ENSG00000000003
    ##  5 ENST00000614008 ENSG00000000003
    ##  6 ENST00000373031 ENSG00000000005
    ##  7 ENST00000485971 ENSG00000000005
    ##  8 ENST00000371588 ENSG00000000419
    ##  9 ENST00000466152 ENSG00000000419
    ## 10 ENST00000371582 ENSG00000000419
    ## # … with 249,148 more rows

Tables are saved in [tibble](http://tibble.tidyverse.org) format,
pipe-able with [dplyr](http://dplyr.tidyverse.org):

``` r
grch38 %>% 
    dplyr::filter(biotype == "protein_coding" & chr == "1") %>% 
    dplyr::select(ensgene, symbol, chr, start, end, description) %>% 
    head %>% 
    knitr::kable(.)
```

| ensgene         | symbol   | chr |     start |       end | description                                                                         |
| :-------------- | :------- | :-- | --------: | --------: | :---------------------------------------------------------------------------------- |
| ENSG00000000457 | SCYL3    | 1   | 169849631 | 169894267 | SCY1 like pseudokinase 3 \[Source:HGNC Symbol;Acc:HGNC:19285\]                      |
| ENSG00000000460 | C1orf112 | 1   | 169662007 | 169854080 | chromosome 1 open reading frame 112 \[Source:HGNC Symbol;Acc:HGNC:25565\]           |
| ENSG00000000938 | FGR      | 1   |  27612064 |  27635185 | FGR proto-oncogene, Src family tyrosine kinase \[Source:HGNC Symbol;Acc:HGNC:3697\] |
| ENSG00000000971 | CFH      | 1   | 196651878 | 196747504 | complement factor H \[Source:HGNC Symbol;Acc:HGNC:4883\]                            |
| ENSG00000001460 | STPG1    | 1   |  24356999 |  24416934 | sperm tail PG-rich repeat containing 1 \[Source:HGNC Symbol;Acc:HGNC:28070\]        |
| ENSG00000001461 | NIPAL3   | 1   |  24415802 |  24472976 | NIPA like domain containing 3 \[Source:HGNC Symbol;Acc:HGNC:25233\]                 |

Example with
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
results from the
[airway](https://bioconductor.org/packages/release/data/experiment/html/airway.html)
package, made tidy with
[biobroom](http://www.bioconductor.org/packages/devel/bioc/html/biobroom.html):

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:rlang':
    ## 
    ##     exprs

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## Loading required package: BiocParallel

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

``` r
library(airway)

data(airway)
airway <- DESeqDataSet(airway, design = ~cell + dex)
airway <- DESeq(airway)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
res <- results(airway)

# tidy results with biobroom
library(biobroom)
```

    ## Loading required package: broom

    ## Registered S3 methods overwritten by 'biobroom':
    ##   method      from 
    ##   glance.list broom
    ##   tidy.list   broom

``` r
res_tidy <- tidy.DESeqResults(res)
head(res_tidy)
```

    ## # A tibble: 6 x 7
    ##   gene            baseMean estimate stderror statistic   p.value p.adjusted
    ##   <chr>              <dbl>    <dbl>    <dbl>     <dbl>     <dbl>      <dbl>
    ## 1 ENSG00000000003  709.      0.381     0.101     3.79   0.000152    0.00128
    ## 2 ENSG00000000005    0      NA        NA        NA     NA          NA      
    ## 3 ENSG00000000419  520.     -0.207     0.112    -1.84   0.0653      0.197  
    ## 4 ENSG00000000457  237.     -0.0379    0.143    -0.264  0.792       0.911  
    ## 5 ENSG00000000460   57.9     0.0882    0.287     0.307  0.759       0.895  
    ## 6 ENSG00000000938    0.318   1.38      3.50      0.394  0.694      NA

``` r
res_tidy %>% 
    dplyr::arrange(p.adjusted) %>% 
    head(20) %>% 
    dplyr::inner_join(grch38, by = c("gene" = "ensgene")) %>% 
    dplyr::select(gene, estimate, p.adjusted, symbol, description) %>% 
    knitr::kable(.)
```

| gene            |   estimate | p.adjusted | symbol  | description                                                                                                           |
| :-------------- | ---------: | ---------: | :------ | :-------------------------------------------------------------------------------------------------------------------- |
| ENSG00000152583 | \-4.574919 |          0 | SPARCL1 | SPARC like 1 \[Source:HGNC Symbol;Acc:HGNC:11220\]                                                                    |
| ENSG00000165995 | \-3.291062 |          0 | CACNB2  | calcium voltage-gated channel auxiliary subunit beta 2 \[Source:HGNC Symbol;Acc:HGNC:1402\]                           |
| ENSG00000120129 | \-2.947810 |          0 | DUSP1   | dual specificity phosphatase 1 \[Source:HGNC Symbol;Acc:HGNC:3064\]                                                   |
| ENSG00000101347 | \-3.766995 |          0 | SAMHD1  | SAM and HD domain containing deoxynucleoside triphosphate triphosphohydrolase 1 \[Source:HGNC Symbol;Acc:HGNC:15925\] |
| ENSG00000189221 | \-3.353580 |          0 | MAOA    | monoamine oxidase A \[Source:HGNC Symbol;Acc:HGNC:6833\]                                                              |
| ENSG00000211445 | \-3.730403 |          0 | GPX3    | glutathione peroxidase 3 \[Source:HGNC Symbol;Acc:HGNC:4555\]                                                         |
| ENSG00000157214 | \-1.976773 |          0 | STEAP2  | STEAP2 metalloreductase \[Source:HGNC Symbol;Acc:HGNC:17885\]                                                         |
| ENSG00000162614 | \-2.035665 |          0 | NEXN    | nexilin F-actin binding protein \[Source:HGNC Symbol;Acc:HGNC:29557\]                                                 |
| ENSG00000125148 | \-2.210979 |          0 | MT2A    | metallothionein 2A \[Source:HGNC Symbol;Acc:HGNC:7406\]                                                               |
| ENSG00000154734 | \-2.345604 |          0 | ADAMTS1 | ADAM metallopeptidase with thrombospondin type 1 motif 1 \[Source:HGNC Symbol;Acc:HGNC:217\]                          |
| ENSG00000139132 | \-2.228903 |          0 | FGD4    | FYVE, RhoGEF and PH domain containing 4 \[Source:HGNC Symbol;Acc:HGNC:19125\]                                         |
| ENSG00000162493 | \-1.891217 |          0 | PDPN    | podoplanin \[Source:HGNC Symbol;Acc:HGNC:29602\]                                                                      |
| ENSG00000134243 | \-2.195712 |          0 | SORT1   | sortilin 1 \[Source:HGNC Symbol;Acc:HGNC:11186\]                                                                      |
| ENSG00000179094 | \-3.191750 |          0 | PER1    | period circadian regulator 1 \[Source:HGNC Symbol;Acc:HGNC:8845\]                                                     |
| ENSG00000162692 |   3.692661 |          0 | VCAM1   | vascular cell adhesion molecule 1 \[Source:HGNC Symbol;Acc:HGNC:12663\]                                               |
| ENSG00000163884 | \-4.459128 |          0 | KLF15   | Kruppel like factor 15 \[Source:HGNC Symbol;Acc:HGNC:14536\]                                                          |
| ENSG00000178695 |   2.528175 |          0 | KCTD12  | potassium channel tetramerization domain containing 12 \[Source:HGNC Symbol;Acc:HGNC:14678\]                          |
| ENSG00000198624 | \-2.918436 |          0 | CCDC69  | coiled-coil domain containing 69 \[Source:HGNC Symbol;Acc:HGNC:24487\]                                                |
| ENSG00000107562 |   1.911670 |          0 | CXCL12  | C-X-C motif chemokine ligand 12 \[Source:HGNC Symbol;Acc:HGNC:10672\]                                                 |
| ENSG00000148848 |   1.814543 |          0 | ADAM12  | ADAM metallopeptidase domain 12 \[Source:HGNC Symbol;Acc:HGNC:190\]                                                   |
