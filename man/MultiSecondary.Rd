% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiSecondary.R
\name{MultiSecondary}
\alias{MultiSecondary}
\alias{MultiSecondary,}
\alias{multi.secondary}
\title{Conditional multivariate association analysis using summary statistics}
\usage{
MultiSecondary(gwa.region, LD.ref, snp.ref, R.ref, p.threshold = 5e-08,
  tol = 0.8, traits, v_y = NULL, T2.return = FALSE)
}
\arguments{
\item{gwa.region}{GWAS summary statistics, includes A1, A2 and three columns for each trait: beta, se and N}

\item{LD.ref}{Regional LD matrix including SNPs in gwa.region}

\item{snp.ref}{The reference alleles of SNPs in the reference LD correlation matrix. The names of the vector 
should be SNP names in reference sample}

\item{R.ref}{Shrinkage phenotypic correlation matrix, achieved from \code{load.summary()}}

\item{p.threshold}{P-value threshold in conditional analysis}

\item{tol}{Tolerance for multicollinearity}

\item{traits}{Traits to be analyzed}

\item{v_y}{The variance of the traits}

\item{T2.return}{Returning conditional T2 statistic or not}
}
\value{
The function returns a list with elements of \code{T2.sele}: The conditional test statistic of the selected variants. 
It will be provided if \code{T2.return = TRUE}; \code{p.sele}: The conditional p-value of the selected variants;
\code{b_joint.sele}: The conditional effect size of the selected variants; \code{se_b_joint.sele}: The conditional 
standard error of the selected variants.
}
\description{
This function is developed to implement cMVA based on multivariate results
}
\examples{
\dontrun{
data(example.MultiSecondary)
##### 474 snps around rs905938 #####

## Six-traits cMVA ##
traits <- c("HEIGHT", "BMI", "HIP", "WC", "WHR", "WEIGHT")
MultiSecondary(gwa.region = example.gwas, LD.ref = example.LD, 
				  snp.ref = example.snp.ref, R.ref = example.R.ref,
               p.threshold = 5e-8, tol = 0.9, traits = traits, T2.return = TRUE)

## Three-traits cMVA ##
traits <- c("HEIGHT", "BMI", "HIP")
MultiSecondary(gwa.region = example.gwas, LD.ref = example.LD, 
               snp.ref = example.snp.ref, R.ref = example.R.ref,
               p.threshold = 5e-4, tol = 0.9, traits = traits, T2.return = TRUE)
}
}
\references{
Zheng Ning, Yakov Tsepilov, Sodbo Zh. Sharapov, Alexander K. Grishenko, Masoud Shirali, Peter K. Joshi,
James F. Wilson, Yudi Pawitan, Chris S. Haley, Yurii S. Aulchenko, Xia Shen (2018).
Multivariate discovery, replication, and interpretation of pleiotropic loci using summary association statistics. \emph{Submitted}.
}
\seealso{
\code{MultiSummary}
}
\author{
Zheng Ning, Xia Shen
}
\keyword{analysis}
\keyword{conditional}
\keyword{meta-analysis,}
\keyword{multivariate,}
