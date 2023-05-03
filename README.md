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
myobj <- create_sge_object(file_libcount = "test/library_dependent_counts.tsv.gz",
                           file_allcount = "test/library_independent_counts.tsv.gz",
                           file_valiant_meta = "test/valiant_meta.csv.gz")
myobj

myobj@adapt5 <- "AGCAGCAGCTGACACAAAGT"
myobj@adapt3 <- "CCTCCTCCCCACTCCCTG"

myobj <- format_count(myobj)
myobj <- sge_stats(myobj) 
myobj <- sge_qc_stats(myobj) 

show_stats(myobj)
show_qc_stats(myobj)
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



