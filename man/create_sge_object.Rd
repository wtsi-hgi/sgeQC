% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\name{create_sge_object}
\alias{create_sge_object}
\title{Create a new SGE object}
\usage{
create_sge_object(
  file_libcount,
  file_allcount,
  file_valiant_meta,
  file_vep_anno = NULL,
  file_libcount_hline = 3,
  file_allcount_hline = 3,
  file_valiant_meta_hline = 1,
  file_vep_anno_hline = 1,
  file_libcount_cols = vector(),
  file_allcount_cols = vector(),
  file_valiant_meta_cols = vector(),
  file_vep_anno_cols = vector()
)
}
\arguments{
\item{file_libcount}{QUANTS library-dependent count file, per sequence per count}

\item{file_allcount}{QUANTS library-independent count file, per sequence per count}

\item{file_valiant_meta}{VaLiAnT meta file}

\item{file_vep_anno}{vep annotation file}

\item{file_libcount_hline}{line number of header in library-dependent count file}

\item{file_allcount_hline}{line number of header in library-independent count file}

\item{file_valiant_meta_hline}{line number of header in VaLiAnT meta file}

\item{file_vep_anno_hline}{line number of header in vep annotation file}

\item{file_libcount_cols}{a vector of numbers of selected columns in library-dependent count file, default is none}

\item{file_allcount_cols}{a vector of numbers of selected columns in library-independent count file, default is none}

\item{file_valiant_meta_cols}{a vector of numbers of selected columns in VaLiAnT meta file, default is none}

\item{file_vep_anno_cols}{a vector of numbers of selected columns in vep annotation file, default is none}
}
\value{
An object of class SGE
}
\description{
Create a new SGE object
}
