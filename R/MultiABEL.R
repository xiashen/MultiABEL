#' MultiABEL: Multivariate Genome-Wide Association Analyses
#' 
#' Performing multivariate genome-wide association (MVGWA) analyses. 
#' The modules are compatible with existing *ABEL data formats. The GWA
#' analyses can be done on individual level data or on
#' single-trait GWA summary statistics only. 
#' 
#' For converting data from other formats, see
#' 
#' \code{{convert.snp.illumina}} (Illumina/Affymetrix-like format). This is 
#' our preferred converting function, very extensively tested. Other conversion 
#' functions include: 
#' \code{{convert.snp.text}} (conversion from human-readable GenABEL format),
#' \code{{convert.snp.ped}} (Linkage, Merlin, Mach, and similar files),
#' \code{{convert.snp.mach}} (Mach-format),
#' \code{{convert.snp.tped}} (from PLINK TPED format),
#' \code{{convert.snp.affymetrix}} (BRML-style files).
#' 
#' For converting of GenABEL's data to other formats, see
#' \code{{export.merlin}} (MERLIN and MACH formats), 
#' \code{{export.impute}} (IMPUTE, SNPTEST and CHIAMO formats),
#' \code{{export.plink}} (PLINK format, also exports phenotypic data). 
#' 
#' To load the data, see \code{{load.gwaa.data}}.
#' 
#' For conversion to DatABEL format (used by ProbABEL and some other 
#' GenABEL suite packages), see 
#' \code{{impute2databel}}, 
#' \code{{impute2mach}}, 
#' \code{{mach2databel}}. 
#' 
#' For data managment and manipulations see
#' \code{{merge.gwaa.data}},
#' \code{{merge.snp.data}},
#' \code{{gwaa.data-class}},
#' \code{{snp.data-class}},
#' \code{{snp.names}},
#' \code{{snp.subset}}.
#' 
#' @author Xia Shen
#' 
#' @references 
#' If you use the MultiABEL package in your analysis, please cite the
#' papers in \code{citation("MultiABEL")}.
#' 
#' @seealso \code{GenABEL}, \code{DatABEL}
#'
#'
#' @name MultiABEL
#' @docType package
#' @title Multivariate GWAS in R
#' @aliases MultiABEL multiabel MultiABEL-package
#' @keywords package
#'
NULL