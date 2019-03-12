#' Multivariate genome-wide association scan using summary statistics
#' 
#' This function performs multivariate GWA analysis using meta-GWAS summary statistics
#' 
#' @param x A data object of class \code{multi.summary} loaded by the function \code{load.summary}.
#' @param index A numeric vector that gives the indices of the traits to be analyzed jointly.
#' @param type A string gives the type of analysis. Default is \code{"outbred"}, referring to 
#' general outbred populations, following Hardy-Weinberg equilibrium. \code{"inbred"} refers to 
#' inbred populations, where no heterzygotes exists, namely, allele frequency = genotype frequency.
#' \code{"precise"} refers to precise test statistics, especially when the individual-level data 
#' are available, for which the argument \code{vars} has to be given. \code{"direct"} refers to 
#' test statistics directly constructed from the T-statistics in univariate GWAS, this provides a
#' scale-invariant test most similar to the direct MANOVA, but may be less powerful in some scenarios.
#' @param vars A numeric vector gives the variance of the genotypes at each SNP, e.g. coded as 0, 1 and 2.
#' Only used when \code{type = "precise"}.
#' @param high.dim Are the phenotypes high-dimensional or not? This is particularly important when the ratio
#' of the number of individuals (n) to the number of phenotypes being analyzed (p) is not big enough, e.g when
#' analyzing a big number of omics phenotypes in a small cohort. Default = \code{FALSE}.
#' 
#' @return The function returns a data frame containing the multi-trait GWAS results, where the row names are
#' the variants names. The column names are: variant name (\code{Marker}), allele frequency (\code{Freq}),
#' the effective sample size of the multiple traits (\code{N}), effect on the phenotype score (\code{Beta.S}, see reference),
#' standard error (\code{SE}), p-value (\code{P}), and the rest the coefficients to construct the phenotype score
#' (see reference).
#' 
#' @author Xia Shen
#' 
#' @references 
#' Zheng Ning, Yakov Tsepilov, Peter K. Joshi,
#' James F. Wilson, Yudi Pawitan, Chris S. Haley, 
#' Yurii S. Aulchenko, Xia Shen (2018).
#' Pleiotropic meta-analysis for genomic studies: discovery, replication, 
#' and interpretation. \emph{Submitted}.
#' 
#' @seealso 
#' \code{load.summary}
#' 
#' @examples 
#' \dontrun{
#' ## download the six example files from:
#' ## https://www.dropbox.com/sh/hhta45cewvvea2s/AADfj4OXlbroToZAwIii2Buha?dl=0
#' ## the summary statistics from Randall et al. (2013) PLoS Genet
#' ## for males only
#' ## bmi: body mass index
#' ## hip: hip circumference
#' ## wc: waist circumference
#' ## whr: waist-hip ratio
#' 
#' ## load the prepared set of independent SNPs
#' indep.snps <- as.character(read.table('indep.snps')$V1)
#' 
#' ## load summary statistics of the six traits
#' stats.male <- load.summary(files = c('bmi.txt', 'height.txt', 
#'                            'weight.txt', 'hip.txt', 'wc.txt', 
#'                            'whr.txt'), indep.snps = indep.snps)
#' 
#' ## perform multi-trait meta-GWAS
#' result <- MultiSummary(stats.male)
#' head(result)
#' }
#' @aliases MultiSummary, multi.summary
#' @keywords multivariate, meta-analysis
#' 
`MultiSummary` <- function(x, index = NULL,  type = 'direct', vars = NULL, high.dim = FALSE) {
	if (class(x) != 'multi.summary') {
		stop('Wrong data type!')
	}
	cat('Multi-trait genome scan ... ')
	if (!is.null(index) & length(index) < nrow(x$cor.pheno)) {
		m <- nrow(x$cor.pheno)
		x$cor.pheno <- x$cor.pheno[index,index]
		x$var.pheno <- x$var.pheno[index]
		idx <- which(!(1:m %in% index))
		x$gwa <- x$gwa[,-c(idx*2 - 1, idx*2)]
	}
	m <- nrow(x$cor.pheno)
	f <- x$gwa[,2*m + 1]
	n <- x$gwa[,2*m + 2]
	k <- nrow(x$gwa)
	betamat <- x$gwa[,seq(1, 2*m, 2)]
	semat <- x$gwa[,seq(2, 2*m, 2)]
	tmat <- betamat/semat
	bad <- which(rowSums(is.na(tmat)) > 0)
	if (length(bad) > 0) {
		betamat <- betamat[-bad,]
		semat <- semat[-bad,]
		tmat <- tmat[-bad,]
		x$gwa <- x$gwa[-bad,]
		k <- nrow(x$gwa)
	}
	dimnames(betamat) <- list(NULL, NULL)
	res <- data.frame(marker = rownames(x$gwa), freq = f, n = n,
		stringsAsFactors = FALSE)
	rownames(res) <- rownames(x$gwa)
	if (type != 'direct') {
		bnames <- paste0('beta.cond.', rownames(x$cor.pheno))
		snames <- paste0('se.cond.', rownames(x$cor.pheno))
		pnames <- paste0('p.cond.', rownames(x$cor.pheno))
		cns <- rep(NA, m*3)
		cns[seq(1, m*3, 3)] <- bnames
		cns[seq(2, m*3, 3)] <- snames
		cns[seq(3, m*3, 3)] <- pnames
		res3 <- matrix(NA, k, m*3)
		colnames(res3) <- cns
		bnames <- paste0('est.coef.', rownames(x$cor.pheno))
		snames <- paste0('se.coef.', rownames(x$cor.pheno))
		pnames <- paste0('p.coef.', rownames(x$cor.pheno))
		cns[seq(1, m*3, 3)] <- bnames
		cns[seq(2, m*3, 3)] <- snames
		cns[seq(3, m*3, 3)] <- pnames
		resc <- matrix(NA, k, m*3)
		colnames(resc) <- cns
	}
	if (type == 'outbred') {
		scan1 <- .Fortran('MultiSummaryLoopDirect', k = as.integer(k), m = as.integer(m), nn = as.numeric(n), 
				tmat = tmat, invR = solve(x$cor.pheno), 
				pil = numeric(k), fstat = numeric(k), PACKAGE = "MultiABEL")
		fs <- scan1$pil/m/(n - scan1$pil)*(n - m - 1)
		mf <- mean(fs)
		n.eff <- mf*2/(mf - 1) + m + 1 # effective N: fit a F-distribution to the observed fstat using moment estimator
		res$n <- n 
		# res$effective.n <- round(n.eff, digits = 2)
		if (!high.dim) {
			# if (n.eff < 5*m) {
			#	warning('Effective sample size < 5 times number of phenotypes. Consider argument high.dim = TRUE.')
			# }
			# pv <- pf(scan1$fstat, m, n - m - 1, lower.tail = FALSE) 
			# modified to the following line 2019-02-08
			pv <- pf(fs, m, n - m - 1, lower.tail = FALSE)
		} else {
			pv <- pf(fs, m, n - m - 1, lower.tail = FALSE)	
			# scan1 <- .Fortran('MultiSummaryLoopDirectFstat', k = as.integer(k), m = as.integer(m), nn = as.numeric(n), 
			#		betamat = betamat, invsemat = 1/semat, invse = diag(m), invR = solve(x$cor.pheno), invV = diag(m), 
			#		pil = numeric(k), fstat = numeric(k), PACKAGE = "MultiABEL")
		}
		scan2 <- .Fortran('MultiSummaryLoop', k = as.integer(k), m = as.integer(m), nn = as.numeric(n), 
				f = as.numeric(f), betamat = betamat, R = x$cor.pheno, invR = solve(x$cor.pheno), D = matrix(0, m, m),
				sdY = diag(sqrt(x$var.pheno)), invsdY = diag(1/sqrt(x$var.pheno)), sY = sqrt(x$var.pheno), 
				b = betamat, s = betamat, pil = numeric(k), coef = betamat, ss = betamat, PACKAGE = "MultiABEL")
		res2 <- data.frame(p = pv, beta.score = scan2$pil)
		res2$se.score <- res2$beta.score/sqrt(qchisq(res2$p, 1, lower.tail = FALSE))
	    res3[,seq(1, m*3, 3)] <- scan2$b
		res3[,seq(2, m*3, 3)] <- scan2$s
		res3[,seq(3, m*3, 3)] <- pchisq(scan2$b**2/scan2$s**2, 1, lower.tail = FALSE)
		res <- cbind(res, res2, res3)
		resc[,seq(1, m*3, 3)] <- scan2$coef
		resc[,seq(2, m*3, 3)] <- sqrt(scan2$ss)
		resc[,seq(3, m*3, 3)] <- pchisq(scan2$coef**2/scan2$ss, 1, lower.tail = FALSE)
	} else if (type == 'inbred') {
		scan1 <- .Fortran('MultiSummaryLoopDirect', k = as.integer(k), m = as.integer(m), nn = as.numeric(n), 
				tmat = tmat, invR = solve(x$cor.pheno), 
				pil = numeric(k), fstat = numeric(k), PACKAGE = "MultiABEL")
		fs <- scan1$pil/m/(n - scan1$pil)*(n - m - 1)
		mf <- mean(fs)
		n.eff <- mf*2/(mf - 1) + m + 1 # effective N: fit a F-distribution to the observed fstat using moment estimator
		res$n <- n 
		# res$effective.n <- round(n.eff, digits = 2)
		if (!high.dim) {
			# if (n.eff < 5*m) {
			#	warning('Effective sample size < 5 times number of phenotypes. Consider argument high.dim = TRUE.')
			# }
			# pv <- pf(scan1$fstat, m, n - m - 1, lower.tail = FALSE) 
			# modified to the following line 2019-02-08
			pv <- pf(fs, m, n - m - 1, lower.tail = FALSE)
		} else {
			pv <- pf(fs, m, n.eff - m - 1, lower.tail = FALSE)	
			# scan1 <- .Fortran('MultiSummaryLoopDirectFstat', k = as.integer(k), m = as.integer(m), nn = as.numeric(n), 
			#		betamat = betamat, invsemat = 1/semat, invse = diag(m), invR = solve(x$cor.pheno), invV = diag(m), 
			#		pil = numeric(k), fstat = numeric(k), PACKAGE = "MultiABEL")
		}
		scan2 <- .Fortran('MultiSummaryLoopInbred', k = as.integer(k), m = as.integer(m), nn = as.numeric(n), 
				f = as.numeric(f), betamat = betamat, R = x$cor.pheno, invR = solve(x$cor.pheno), D = matrix(0, m, m),
				sdY = diag(sqrt(x$var.pheno)), invsdY = diag(1/sqrt(x$var.pheno)), sY = sqrt(x$var.pheno), 
				b = betamat, s = betamat, pil = numeric(k), coef = betamat, ss = betamat, PACKAGE = "MultiABEL")
		res2 <- data.frame(p = pv, beta.score = scan2$pil)
		res2$se.score <- res2$beta.score/sqrt(qchisq(res2$p, 1, lower.tail = FALSE))
		res3[,seq(1, m*3, 3)] <- scan2$b
		res3[,seq(2, m*3, 3)] <- scan2$s
		res3[,seq(3, m*3, 3)] <- pchisq(scan2$b**2/scan2$s**2, 1, lower.tail = FALSE)
		res <- cbind(res, res2, res3)
		resc[,seq(1, m*3, 3)] <- scan2$coef
		resc[,seq(2, m*3, 3)] <- sqrt(scan2$ss)
		resc[,seq(3, m*3, 3)] <- pchisq(scan2$coef**2/scan2$ss, 1, lower.tail = FALSE)
	} else if (type == 'precise' & !is.null(vars)) {
		scan1 <- .Fortran('MultiSummaryLoopDirect', k = as.integer(k), m = as.integer(m), nn = as.numeric(n), 
				tmat = tmat, invR = solve(x$cor.pheno), 
				pil = numeric(k), fstat = numeric(k), PACKAGE = "MultiABEL")
		fs <- scan1$pil/m/(n - scan1$pil)*(n - m - 1)
		mf <- mean(fs)
		n.eff <- mf*2/(mf - 1) + m + 1 # effective N: fit a F-distribution to the observed fstat using moment estimator
		res$n <- n 
		# res$effective.n <- round(n.eff, digits = 2)
		if (!high.dim) {
			# if (n.eff < 5*m) {
			#	warning('Effective sample size < 5 times number of phenotypes. Consider argument high.dim = TRUE.')
			# }
			# pv <- pf(scan1$fstat, m, n - m - 1, lower.tail = FALSE) 
			# modified to the following line 2019-02-08
			pv <- pf(fs, m, n - m - 1, lower.tail = FALSE)
		} else {
			pv <- pf(fs, m, n.eff - m - 1, lower.tail = FALSE)	
			# scan1 <- .Fortran('MultiSummaryLoopDirectFstat', k = as.integer(k), m = as.integer(m), nn = as.numeric(n), 
			#		betamat = betamat, invsemat = 1/semat, invse = diag(m), invR = solve(x$cor.pheno), invV = diag(m), 
			#		pil = numeric(k), fstat = numeric(k), PACKAGE = "MultiABEL")
		}
		scan2 <- .Fortran('MultiSummaryLoopPrecise', k = as.integer(k), m = as.integer(m), nn = as.numeric(n), 
				varg = vars, betamat = betamat, R = x$cor.pheno, invR = solve(x$cor.pheno), D = matrix(0, m, m),
				sdY = diag(sqrt(x$var.pheno)), invsdY = diag(1/sqrt(x$var.pheno)), sY = sqrt(x$var.pheno), 
				b = betamat, s = betamat, pil = numeric(k), coef = betamat, ss = betamat, PACKAGE = "MultiABEL")
		res2 <- data.frame(p = pv, beta.score = scan2$pil)
		res2$se.score <- res2$beta.score/sqrt(qchisq(res2$p, 1, lower.tail = FALSE))
		res3[,seq(1, m*3, 3)] <- scan2$b
		res3[,seq(2, m*3, 3)] <- scan2$s
		res3[,seq(3, m*3, 3)] <- pchisq(scan2$b**2/scan2$s**2, 1, lower.tail = FALSE)
		res <- cbind(res, res2, res3)
		resc[,seq(1, m*3, 3)] <- scan2$coef
		resc[,seq(2, m*3, 3)] <- sqrt(scan2$ss)
		resc[,seq(3, m*3, 3)] <- pchisq(scan2$coef**2/scan2$ss, 1, lower.tail = FALSE)
	} else if (type == 'direct') {
		scan1 <- .Fortran('MultiSummaryLoopDirect', k = as.integer(k), m = as.integer(m), nn = as.numeric(n), 
				tmat = tmat, invR = solve(x$cor.pheno), 
				pil = numeric(k), fstat = numeric(k), PACKAGE = "MultiABEL")
		fs <- scan1$pil/m/(n - scan1$pil)*(n - m - 1)
		mf <- mean(fs)
		n.eff <- mf*2/(mf - 1) + m + 1 # effective N: fit a F-distribution to the observed fstat using moment estimator
		res$n <- n 
		# res$effective.n <- round(n.eff, digits = 2)
		if (!high.dim) {
			# if (n.eff < 5*m) {
			#	warning('Effective sample size < 5 times number of phenotypes. Consider argument high.dim = TRUE.')
			# }
			# pv <- pf(scan1$fstat, m, n - m - 1, lower.tail = FALSE) 
			# modified to the following line 2019-02-08
			pv <- pf(fs, m, n - m - 1, lower.tail = FALSE)
		} else {
			pv <- pf(fs, m, n.eff - m - 1, lower.tail = FALSE)	
			# scan1 <- .Fortran('MultiSummaryLoopDirectFstat', k = as.integer(k), m = as.integer(m), nn = as.numeric(n), 
			#		betamat = betamat, invsemat = 1/semat, invse = diag(m), invR = solve(x$cor.pheno), invV = diag(m), 
			#		pil = numeric(k), fstat = numeric(k), PACKAGE = "MultiABEL")
		}
		res$p <- pv
	}else {
		stop('Wrong type of analysis!')
	}
	cat('Done.\n')
	if (type %in% c('outbred', 'inbred', 'precise')) {
		colnames(scan2$coef) <- rownames(x$cor.pheno)
		rownames(scan2$coef) <- rownames(res)
		return(list(scan = res, coef = resc))
	} else {
		return(list(scan = res, coef = NULL))
	}
}
