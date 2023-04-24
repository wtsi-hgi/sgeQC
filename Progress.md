## Progress
Based on user stories from https://docs.google.com/spreadsheets/d/15wh0OaS2DfawpNxuVeQ4yjGkzkevN43O27bUffV2bWw/edit#gid=0

### 1. reading files

#### Ingestion of library-independent counts
Essential:
- [x] read count tsv file
- [x] returns a data frame 

Preferable:
- [x] can read in uncompressed and compressed input files
- [x] user can define column delimiter (e.g. tab, comma)
- [x] user can define any number of comment lines at the start which can be skipped
- [x] user can define whether the file has a header (or not)
- [ ] column containing read sequence and counts is configurable
- [ ] allows sample (column) names which start with or are numeric values