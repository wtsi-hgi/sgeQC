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
```R
install.packages("devtools")

library(devtools)
install_github("wtsi-hgi/sgeQC")
```

Install from source file
```R
install.packages("/path/of/sgeQC.tar.gz", type = "source")
```

<!-- Dependencies-->
## Dependencies

```R
install.packages("data.table")
install.packages("Ckmeans.1d.dp")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("gplots")
```

<!-- Import Data-->
## Import Data

### Import from a directory:
All the files are in the same directory including library dependent counts, library independent counts, valiant meta csv, vep annotation and the sample sheet.

#### Sample sheet format
| sample_name  | library_independent_count | library_dependent_count | valiant_meta | adapt5 | adapt3 | library_name | library_type|
| - | - | - | - | - | - | - | - |
| sample1 | s1.allcounts.tsv.gz | s1.libcounts.tsv.gz | meta.csv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | libA | screen |
| sample2 | s2.allcounts.tsv.gz | s2.libcounts.tsv.gz | meta.csv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | libA | screen |
| sample3 | s3.allcounts.tsv.gz | s3.libcounts.tsv.gz | meta.csv.gz | CTGACTGGCACCTCTTCCCCCAGGA | CCCCGACCCCTCCCCAGCGTGAATG | libA | screen |

* *please use the same headers in the example*
* *adapt5 and adapt3 are required if you don't provide the ref seq and pam seq*
* *library_name and library_type are not necessary*

```R
library(sgeQC)

sge_objs <- import_sge_files("/path/to/input/directory", "sample_sheet.tsv")
samqc <- create_sampleqc_object(sge_objs)
```

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Plasmid QC -->
## Plasmid QC

<a id="pqc1"></a>
### QC 1: 
```R
samqc <- run_sample_qc(samqc, "plasmid")

qcplot_readlens(samqc, "/path/to/out/plots")
qcplot_clusters(samqc, "plasmid", "/path/to/out/plots")
qcplot_stats(samqc, "/path/to/out/plots")
qcplot_position(samqc, "/path/to/out/plots")

qcout_bad_seqs(samqc, "/path/to/out/files")
qcout_sampleqc_stats(samqc, "/path/to/out/files")
```

<p align="right">(<a href="#top">TOP</a>)</p>

<a id="pqc2"></a>
### QC 2: 

<p align="right">(<a href="#top">TOP</a>)</p>

<!-- Screen QC -->
## Screen QC

<a id="sqc1"></a>
### QC 1: 

```R
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

<a id="sqc2"></a>
### QC 2: 

<p align="right">(<a href="#top">TOP</a>)</p>

