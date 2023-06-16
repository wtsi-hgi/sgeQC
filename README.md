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

```sh
suppressMessages(library(data.table))
library(Ckmeans.1d.dp)
library(reshape2)
library(ggplot2)
library(gplots)
```

<!-- Import Data-->
## Import Data


```sh
library(sgeQC)

sge_objs <- import_sge_files("/path/to/input/directory", "sample_sheet.tsv")

samqc <- create_sampleqc_object(sge_objs)
samqc@samples_ref <- select_objects(sge_objs, c(2,5,8))
samqc <- run_sample_qc(samqc, "screen")

qcplot_readlens(samqc, "/path/to/out/plots")
qcplot_clusters(samqc, "screen", "/path/to/out/plots")
qcplot_stats(samqc, "/path/to/out/plots")
qcplot_position(samqc, "/path/to/out/plots")

qcout_bad_seqs(samqc, "/path/to/out/files")
qcout_sampleqc_stats(samqc, "/path/to/out/files")
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

