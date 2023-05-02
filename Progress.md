<div align="center">
<h1 align="center">Progress</h1>
</div>

Based on user stories from [the google sheet](https://docs.google.com/spreadsheets/d/15wh0OaS2DfawpNxuVeQ4yjGkzkevN43O27bUffV2bWw/edit#gid=0)

- [x] tested and passed
- [ ] not available yet

## Table of Contents
<details hide>
<summary><b>(click to expand or hide)</b></summary>

- [Ingestion](#ingestion)
  + [1. Ingestion of library-independent counts](#1-ingestion-of-library-independent-counts)
  + [2. Ingestion of library dependent counts](#2-ingestion-of-library-dependent-counts)
  + [3. Ingestion of library annotation](#3-ingestion-of-library-annotation)
- [Processing](#processing)
  + [1. Summary statistics per sample for library-independent counts](#1-summary-statistics-per-sample-for-library-independent-counts)
  + [2. Summary statistics per sample for library-dependent counts](#2-summary-statistics-per-sample-for-library-dependent-counts)
  + [3. Summary of unique read sequence length distribution](#3-summary-of-unique-read-sequence-length-distribution)
  + [4. Percentage of reads mapping to REF](#4-percentage-of-reads-mapping-to-ref)
  + [5. Percentage of reads mapping to PAM](#5-percentage-of-reads-mapping-to-pam)
  + [6. Percentage of reads mapping to library](#6-percentage-of-reads-mapping-to-library)
  + [7. Percentage of unmapped reads](#7-percentage-of-unmapped-reads)
  + [8. Percentage of variants missing](#8-percentage-of-variants-missing)
  + [9. Gini coefficient](#9-gini-coefficient)
  + [10. Percentage of adapters](#10-percentage-of-adapters)
  + [11. Percentage of paired reads with correct insert size](#11-percentage-of-paired-reads-with-correct-insert-size)
  + [12. Reads per sample](#12-reads-per-sample)
  + [13. Mapping statistic](#13-mapping-statistic)
  + [14. Beeswarm](#14-beeswarm)

</details>

## Ingestion
### 1. Ingestion of library-independent counts
**Essential:**
- [x] read count tsv file
- [x] returns a data frame 

**Preferable:**
- [x] can read in uncompressed and compressed input files
- [x] user can define column delimiter (e.g. tab, comma)
- [x] user can define any number of comment lines at the start which can be skipped
- [x] user can define whether the file has a header (or not)
- [x] column containing read sequence and counts is configurable
- [ ] allows sample (column) names which start with or are numeric values

### 2. Ingestion of library dependent counts
**Essential:**
- [x] read count tsv file
- [x] returns a data frame 

**Preferable:**
- [x] can read in uncompressed and compressed input files
- [x] user can define column delimiter (e.g. tab, comma)
- [x] user can define any number of comment lines at the start which can be skipped
- [x] user can define whether the file has a header (or not)
- [x] column containing read sequence and counts is configurable
- [ ] allows sample (column) names which start with or are numeric values

### 3. Ingestion of library annotation
**Essential:**
- [x] read meta csv file
- [x] returns a data frame 

**Preferable:**
- [x] user can define column delimiter (e.g. tab, comma)
- [x] user can define any number of comment lines at the start which can be skipped
- [x] user can define whether the file has a header (or not)
- [x] user can define whether to keep all columns or just a subset
- [x] user can define a column which contains the oligo identifiers and can check that identifiers have not been duplicated (i.e. unique ids)
- [ ] allows sample (column) names which start with or are numeric values

## Processing
### 1. Summary statistics per sample for library-independent counts
**Essential:**
- [x] Total number of quantified reads (sum of all unique read sequence counts)
- [x] Number of unique read sequences
- [x] Output as a dataframe (column: statistics, row: sample)

**Preferable:**
- [x] Maximum unique read sequence length
- [x] Minimum unique read sequence length
- [x] Mean unique read count
- [x] Median unique read count
- [x] Number of unique reads with count > user-defined value (library-independent low counts)
- [ ] User-defined output format: dataframe, list or JSON

### 2. Summary statistics per sample for library-dependent counts
**Essential:**
- [x] Total number of quantified reads (sum of all unique read sequence counts)
- [x] Number of unique read sequences
- [x] Output as a dataframe (column: statistics, row: sample)

**Preferable:**
- [x] Maximum unique read sequence length
- [x] Minimum unique read sequence length
- [x] Mean unique read count
- [x] Median unique read count
- [x] Number of unique reads with count > user-defined value (library-independent low counts)
- [ ] User-defined output format: dataframe, list or JSON

### 3. Summary of unique read sequence length distribution
```diff
+ need mulitqc report of read length
```
**Essential:**
- [ ] User can define a list of sequence length bins
- [ ] Function returns number of reads in each sequence length bin

### 4. Percentage of reads mapping to REF 
**Essential:**
- [x] percentage of REF reads
- [ ] sample ID or reference

**Preferable:**
- [ ] flag to exclude sample, based on predefined threshold

### 5. Percentage of reads mapping to PAM 
**Essential:**
- [x] percentage of PAM reads
- [ ] sample ID or reference

**Preferable:**
- [ ] flag to exclude sample, based on predefined threshold

### 6. Percentage of reads mapping to library
**Essential:**
- [x] percentage of Effective library reads
- [ ] sample ID or reference

**Preferable:**
- [ ] flag to exclude sample, based on predefined threshold

### 7. Percentage of unmapped reads
**Essential:**
- [x] percentage of unmappped reads
- [ ] sample ID or reference

**Preferable:**
- [ ] flag to exclude sample, based on predefined threshold

### 8. Percentage of variants missing
**Essential:**
- [ ] percentage of variants missing
- [ ] sample ID or reference

**Preferable:**
- [ ] flag to exclude sample, based on predefined threshold

### 9. Gini coefficient
**Essential:**
- [ ] gini-coefficient
- [ ] sample ID or reference

**Preferable:**
- [ ] flag to exclude sample, based on predefined threshold

### 10. Percentage of adapters
**Essential:**

**Preferable:**

### 11. Percentage of paired reads with correct insert size
**Essential:**

**Preferable:**

### 12. Reads per sample
**Essential:**
- [ ] number of reads
- [ ] sample ID or reference

**Preferable:**
- [ ] flag to exclude sample, based on predefined threshold

### 13. Mapping statistic
**Essential:**
- [ ] mapping score
- [ ] percentage of samples passing QC
- [ ] distribution of samples passing QC
- [ ] experiment ID or reference

**Preferable:**
- [ ] flag to exclude experiment, based on predefined criteria

### 14. Beeswarm
**Essential:**
- [ ] percentage of mutants with no functional impact
- [ ] experiment ID or reference

**Preferable:**
- [ ] flag to exclude experiment, based on predefined criteria
