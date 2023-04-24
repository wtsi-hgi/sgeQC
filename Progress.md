<div align="center">
<h1 align="center">Progress</h1>
</div>

Based on user stories from [the google sheet](https://docs.google.com/spreadsheets/d/15wh0OaS2DfawpNxuVeQ4yjGkzkevN43O27bUffV2bWw/edit#gid=0)

- [x] tested and passed
- [ ] not available yet

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
- [ ] column containing read sequence and counts is configurable
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
- [ ] column containing read sequence and counts is configurable
- [ ] allows sample (column) names which start with or are numeric values

### 3. Ingestion of library annotation
**Essential:**
- [x] read meta csv file
- [x] returns a data frame 

**Preferable:**
- [x] user can define column delimiter (e.g. tab, comma)
- [x] user can define any number of comment lines at the start which can be skipped
- [x] user can define whether the file has a header (or not)
- [ ] user can define whether to keep all columns or just a subset
- [ ] user can define a column which contains the oligo identifiers and can check that identifiers have not been duplicated (i.e. unique ids)
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