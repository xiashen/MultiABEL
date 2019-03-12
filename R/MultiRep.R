#' Replication analysis of multivariate genome-wide association signal
#' 
#' The function performs replication analysis of multivariate GWA signals.
#' 
#' @param training.pheno An (optional) matrix or data frame contains the phenotype data for the discovery
#' sample, preferrably adjusted for fixed effects and population structure before multivariate GWA analysis.
#' @param training.phenofile An (optional) plain text file contains phenotypes for the discovery sample. 
#' If this is provided, it will serve as \code{training.pheno}.
#' @param test.pheno An (optional) matrix or data frame contains the phenotype data for the replication
#' sample, preferrably adjusted for fixed effects and population structure.
#' @param test.phenofile An (optional) plain text file contains phenotypes of the replication sample. 
#' If this is provided, it will serve as \code{test.pheno}.
#' @param pheno.names A vector (length > 1) giving the column names of the phenotypes to be analyzed. 
#' @param training.geno A matrix or data.frame that contains the discovery sample genotype dosages 
#' of the variants to replicate.
#' @param test.geno A matrix or data.frame that contains the replication sample genotype dosages 
#' of the variants to replicate. This object should have the same column names and order 
#' as \code{training.geno}.
#' 
#' @note Either \code{.pheno} or \code{.phenofile} has to be provided.
#' If both are provided, only \code{phenofile} will be used. Individual IDs 
#' in \code{.pheno} or \code{.phenofile} and \code{.geno} have to match!
#' 
#' @return The function returns a list of 3 matrices. \code{$replication} contains the estimate of 
#' variant effect on the corresponding compound phenotype (\code{beta_c}), standard error (\code{s.e.}), 
#' replication P-value (\code{P}), and proportion of phenotypic variance explained (\code{R-squared}).
#' \code{$training.coef} contains the estimated coefficients in the discovery sample of each phenotype 
#' for each variant to construct the compound phenotype. \code{$test.coef} contains similar coefficients 
#' as in \code{$training.coef} but estimated in the replication sample, but these are just for the record,
#' NOT used in the replication procedure.
#' 
#' @author Xia Shen
#' 
#' @references 
#' Shen X, Klaric L, Sharapov S, Mangino M, Ning Z, Wu D, 
#' Trbojevic-Akmacic I, Pucic-Bakovic M, Rudan I, Polasek O, 
#' Hayward C, Spector TD, Wilson JF, Lauc G, Aulchenko YS (2017): 
#' Multivariate discovery and replication of five novel loci 
#' associated with Immunoglobulin G N-glycosylation. 
#' \emph{Nature Communications}, 8, 447; doi: 10.1038/s41467-017-00453-3.
#' 
#' @seealso 
#' \code{\link{Multivariate}}
#' 
#' @examples 
#' \dontrun{
#' ## loading example discovery sample gwaa.data in GenABEL
#' data(ge03d2)
#' 
#' ## running multivariate GWAS for 3 traits: height, weight, bmi
#' res <- Multivariate(gwaa.data = ge03d2, trait.cols = c(5, 6, 8), 
#'                     covariate.cols = c(2, 3))
#' 
#' ## extracting 5 significant variants
#' (top <- res[order(res[,'P.F']),][2:6,])
#' snps <- rownames(top)
#' training.geno <- as.double(gtdata(ge03d2)[,snps])
#' 
#' ## loading example test sample gwaa.data in GenABEL
#' data(ge03d2c)
#' 
#' ## extracting genotypes of the 5 variants
#' test.geno <- as.double(gtdata(ge03d2c)[,snps])
#' 
#' ## try replication
#' rep <- MultiRep(training.pheno = phdata(ge03d2), test.pheno = phdata(ge03d2c), 
#'                 pheno.names = c('height', 'weight', 'bmi'),
#'                    training.geno = training.geno, test.geno = test.geno)
#' }
#' @aliases MultiRep, multirep
#' @keywords multivariate, replication
#' 
`MultiRep` <- function(training.pheno = NULL, training.phenofile = NULL, 
        test.pheno = NULL, test.phenofile = NULL, pheno.names = NULL,
        training.geno, test.geno) {
    cat('loading data ...')
    if (is.null(training.phenofile) & is.null(training.pheno)) {
        stop('no phenotypes provided for the training set!')
    }
    if (is.null(test.phenofile) & is.null(test.pheno)) {
        stop('no phenotypes provided for the test set!')
    }
    if (!is.null(training.phenofile)) {
        training.pheno <- read.table(training.phenofile, header = TRUE)
    } 
    if (!is.null(test.phenofile)) {
        test.pheno <- read.table(test.phenofile, header = TRUE)
    }
    if (is.null(pheno.names)) {
        idx <- which(colnames(training.pheno) %in% colnames(test.pheno))
        pheno.names <- colnames(training.pheno)[idx]
    }
    if (!all(pheno.names %in% colnames(training.pheno)) |
        !all(pheno.names %in% colnames(test.pheno))) stop('phenotype(s) missing!')
    if (ncol(training.geno) != ncol(test.geno) |
        nrow(training.geno) != nrow(training.pheno) |
        nrow(test.geno) != nrow(test.pheno)) stop('dimensions do not match!')
    cat(' OK\n')
    cat('initializing analysis ...')
    Y0 <- as.matrix(training.pheno[,pheno.names])
    Y1 <- as.matrix(test.pheno[,pheno.names])
    training.coef <- test.coef <- matrix(NA, length(pheno.names), ncol(training.geno))
    dimnames(training.coef) <- dimnames(test.coef) <- list(pheno.names, colnames(training.geno))
    reptab <- matrix(NA, ncol(training.geno), 4)
    dimnames(reptab) <- list(colnames(training.geno), c('beta_c', 's.e.', 'P', 'R-squared'))
    cat(' OK\n')
    if (ncol(training.geno) != 1) {
        cat('replication analysis ...\n')
    } else {
		cat('replication analysis ...')
	}
    for (k in 1:ncol(training.geno)) {
        m0 <- lm(training.geno[,k] ~ Y0)
        training.coef[,k] <- m0$coef[-1]
        m1 <- lm(test.geno[,k] ~ Y1)
        test.coef[,k] <- m1$coef[-1]
        cy <- Y1 %*% training.coef[,k]
        mrep <- lm(cy ~ test.geno[,k])
        reptab[k,] <- c(summary(mrep)$coef[2,c(1:2, 4)], summary(mrep)$r.squared)
        if (ncol(training.geno) != 1) {
			progress(k/ncol(training.geno)*100)
		}
    }
    if (ncol(training.geno) == 1) {
		cat(' OK\n')
	} else {
		cat('\n')
	}
    return(list(replication = reptab, training.coef = training.coef, test.coef = test.coef))
}
