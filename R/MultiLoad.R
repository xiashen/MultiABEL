#' Load individual-level data for multivariate GWA analysis
#' 
#' The function imports GenABEL (gwaa.data class) or DatABEL (.fv*) data formats
#' to perform multivariate test for each genetic variant.
#' 
#' @param gwaa.data An (optional) object of \code{{gwaa.data-class}}.
#' @param phenofile An (optional) plain text file contains phenotypic outcomes and covariates.
#' @param genofile An (optional) object of \code{{databel-class}} containing genotype data.
#' @param trait.cols A vector (length > 1) giving the column indices of the phenotypes to be analyzed. 
#' @param covariate.cols An (optional) vector giving the column indices of the covariates to be included.
#' @param cuts An integer telling how many pieces the genotype data matrix will be splitted for analyze.
#' The value can be set depending on the memory size. The smaller the value is, potentially the faster
#' the analysis will be.
#' @param impute An (optional) logical argument telling whether missing genotypes should be imputed.
#' @param gaussianize An (optional) logical argument telling whether the phenotypes should be gaussianized
#' via inverse-Gaussian transformation.
#' @param ... not used.
#' 
#' @note Either \code{gwaa.data} (for GenABEL data format) or the combination of 
#' \code{phenofile} and \code{genofile} (for DatABEL data format) has to be provided.
#' If all are provided, only \code{phenofile} and \code{genofile} will be used. When using
#' DatABEL format input, individual IDs in \code{phenofile} and \code{genofile} have to match!
#' 
#' @return The function returns a list of cleaned statistics for subsequent, with class \code{multi.loaded}. 
#' 
#' @author Xia Shen
#' 
#' @references 
#' Xia Shen, ..., Jim Wilson, Gordan Lauc, Yurii Aulchenko (2015).
#' Multi-omic-variate analysis identified novel loci associated with 
#' compound N-Glycosylation of human Immunoglobulin G. \emph{Submitted}.
#' 
#' @seealso 
#' \code{\link{Multivariate}}
#' 
#' @examples 
#' \dontrun{
#' ## loading example gwaa.data in GenABEL
#' require(GenABEL)
#' data(ge03d2ex.clean)
#' 
#' ## running multivariate GWAS for 3 traits: height, weight, bmi
#' loaded <- MultiLoad(gwaa.data = ge03d2ex.clean, trait.cols = c(5, 6, 8), 
#'                     covariate.cols = c(2, 3))
#' 
#' ## converting the same dataset into DatABEL format files
#' require(DatABEL)
#' write.table(phdata(ge03d2ex.clean), 'pheno.txt', col.names = TRUE, row.names = TRUE, 
#'             quote = FALSE, sep = '\t')
#' geno <- as.double(ge03d2ex.clean)
#' matrix2databel(geno, 'geno')
#' 
#' ## running the multivariate GWAS again
#' loaded <- MultiLoad(phenofile = 'pheno.txt', genofile = 'geno', trait.cols = c(5, 6, 8), 
#'                     covariate.cols = c(2, 3))
#' }
#' @aliases MultiLoad, multiload
#' @keywords multivariate, multiload, load
#' 
`MultiLoad` <- function(gwaa.data = NULL, phenofile = NULL, genofile = NULL,  
                           trait.cols, covariate.cols = NULL, cuts = 20, 
						   impute = TRUE, gaussianize = TRUE, ...) {
	set.seed(911)
	#if (!require(GenABEL) | !require(DatABEL)) {
	#	stop('GenABEL and DatABEL packages required!')
	#}
    cat('loading data ...')
    if (length(trait.cols) == 1) {
		stop('select multiple traits to analyze!')
	}
    if (!is.null(phenofile) & !is.null(genofile)) {
        pheno <- read.table(phenofile, header = TRUE)
        geno <- DatABEL::databel(genofile)
        if (nrow(pheno) != nrow(geno)) {
			stop('sizes of phenofile and genofile do not match!')
		}
        GenABEL <- FALSE
        cat(' OK\n')
    } else if (!is.null(gwaa.data)) {
        if (!is.null(phenofile)) {
			pheno <- read.table(phenofile, header = TRUE)
		} else {
			pheno <- GenABEL::phdata(gwaa.data)
		}
        GenABEL <- TRUE
        cat(' OK\n')
    } else {
        stop('insufficient data input!')
    }
    cat('preparing data ...')
    Y <- as.matrix(pheno[,trait.cols]) 
	okidx <- which(!is.na(rowSums(Y)))
	Y <- Y[okidx,]
	n <- length(okidx)
	if (gaussianize) {
		for (i in 1:ncol(Y)) {
			Y[,i] <- zscore(Y[,i])
		}
	}
    X <- matrix(1, length(okidx), 1)
    if (!is.null(covariate.cols)) {
		X <- cbind(X, as.matrix(pheno[okidx,covariate.cols]))
	}
    if (!GenABEL) {
		nsnp <- ncol(geno)
		snp <- colnames(geno)
	} else {
		nsnp <- ncol(gwaa.data@gtdata)
		snp <- gwaa.data@gtdata@snpnames
	}
	bs <- matrix(NA, nsnp, ncol(Y)*2)
	cvg <- numeric(nsnp)
    if (cuts != 1) {
        cutpoints <- round(seq(1, nsnp, length = cuts + 1)) 
        starts <- cutpoints[1:cuts]
        ends <- c(cutpoints[2:cuts] - 1, cutpoints[cuts + 1])
    } else {
        starts <- 1
        ends <- nsnp
    }
    if (cuts != 1) {
        cat('\n')
    }
    for (k in 1:cuts) {
        idx <- starts[k]:ends[k]
        if (!GenABEL) {
			g <- DatABEL::databel2matrix(geno[,idx])[okidx,]
		} else {
			g <- as.double(gwaa.data@gtdata[,idx])[okidx,]
		}
		if (any(is.na(g))) {
			if (impute) {
				for (j in which(is.na(colSums(g)))) {
					naidx <- which(is.na(g[,j]))
					g[naidx,j] <- sample(g[-naidx,j], length(naidx), replace = TRUE)
				}
			}
		}
		cvg[idx] <- colVars(g)
		U3 = crossprod(X, g)
		U4 = solve(crossprod(X), U3)
		Str = g - X %*% U4
		Str2 = colSums(Str ^ 2)
		for (i in 1:ncol(Y)) {
			y <- Y[,i]
			U1 = crossprod(X, y)
			U2 = solve(crossprod(X), U1)
			ytr = y - X %*% U2
			b = as.vector(crossprod(ytr, Str) / Str2)
			## calculate residual error
			sig = (sum(ytr ^ 2) - b ^ 2 * Str2) / (n - k - 2)
			## calculate standard error for beta
			err = sqrt(sig * (1 / Str2))
			bs[idx,i*2 - 1] <- b
			bs[idx,i*2] <- err
		}
        if (cuts != 1) {
			progress(k/cuts*100)
		}
    }
    if (cuts == 1) {
		cat(' OK\n')
	} else {
		cat('\n')
	}
	R <- cor(Y)
	vy <- colSums((t(t(Y) - colMeans(Y)))^2)/(nrow(Y) - 1)
	rownames(bs) <- snp
    res <- list(gwa = bs, cor.pheno = R, var.pheno = vy, cvg = cvg, snp = snp, nid = nrow(Y), nsnp = nsnp)
	class(res) <- 'multi.loaded'
    return(res)
}
