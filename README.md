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
    - [QC 1: ](#qc1)
    - [QC 2: ](#qc2)

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
```

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Plasmid QC -->
## Plasmid QC

<a id="qc1"></a>
### QC 1: 

<p align="right">(<a href="#top">TOP</a>)</p>



