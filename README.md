<div align="center">
<h1 align="center">sgeQC</h1>
  <p align="center">A R package of quality control on SGE data</p>
</div>

## Table of Contents
<details open>
<summary><b>(click to expand or hide)</b></summary>

1. [Installation](#installation)
2. [Dependencies](#dependencies)
3. [Import Data](#import-data)
4. [Plasmid QC](#plasmid-qc)
    - [QC 1: ](#pqc1)
    - [QC 2: ](#pqc2)
5. [Screen QC](#screen-qc)
    - [QC 1: ](#sqc1)
    - [QC 2: ](#sqc2)

</details>

<!-- Installation-->
## Installation

Install from github
```sh
install.packages("devtools")

library(devtools)
install_github("wtsi-hgi/sgeQC")
```

Install from source file
```sh
install.packages("/path/of/sgeQC.tar.gz", type = "source")
```

<!-- Dependencies-->
## Dependencies

<!-- Import Data-->
## Import Data


```sh
library(sgeQC)
objA <- create_sge_object(file_libcount = "test/A.library_dependent_counts.tsv.gz",
                          file_allcount = "test/A.library_independent_counts.tsv.gz",
                          file_valiant_meta = "test/A.valiant_meta.csv.gz")
objA

objA@adapt5 <- "AGCAGCAGCTGACACAAAGT"
objA@adapt3 <- "CCTCCTCCCCACTCCCTG"
objA <- format_count(objA)
objA <- sge_stats(objA) 
objA <- sge_qc_stats(objA)
objA@sample <- "A"

show_stats(objA)
show_qc_stats(objA)

objB <- create_sge_object(file_libcount = "test/B.library_dependent_counts.tsv.gz",
                          file_allcount = "test/B.library_independent_counts.tsv.gz",
                          file_valiant_meta = "test/B.valiant_meta.csv.gz")
objB

objB@adapt5 <- "AGCTCTTCAGGTGGTTTTTGGA"
objB@adapt3 <- "TCTTACCAAGTTAAATGCATGATTCGT"
objB@sample <- "B"
objB <- format_count(objB)
objB <- sge_stats(objB) 
objB <- sge_qc_stats(objB)

priqc <- create_primaryqc_object(list(objA, objB))
```

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Plasmid QC -->
## Plasmid QC

<a id="pqc1"></a>
### QC 1: 

<p align="right">(<a href="#top">TOP</a>)</p>

<a id="pqc2"></a>
### QC 2: 

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Screen QC -->
## Screen QC

<a id="sqc1"></a>
### QC 1: 

<p align="right">(<a href="#top">TOP</a>)</p>

<a id="sqc2"></a>
### QC 2: 

<p align="right">(<a href="#top">TOP</a>)</p>

