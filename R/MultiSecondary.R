#' Conditional multivariate association analysis using summary statistics
#' 
#' This function is developed to implement cMVA based on multivariate results
#' 
#' @param gwa.region GWAS summary statistics, includes A1, A2 and three columns for each trait: beta, se and N
#' @param LD.ref Regional LD matrix including SNPs in gwa.region
#' @param snp.ref The reference alleles of SNPs in the reference LD correlation matrix. The names of the vector 
#' should be SNP names in reference sample
#' @param R.ref Shrinkage phenotypic correlation matrix, achieved from \code{MultiSummary()}
#' @param p.threshold P-value threshold in conditional analysis
#' @param tol Tolerance for multicollinearity
#' @param traits Traits to be analyzed
#' @param v_y The variance of the traits
#' @param T2.return Returning conditional T2 statistic or not
#' 
#' @return The function returns a list with elements of \code{T2.sele}: The conditional test statistic of the selected variants. 
#' It will be provided if \code{T2.return = TRUE}; \code{p.sele}: The conditional p-value of the selected variants;
#' \code{b_joint.sele}: The conditional effect size of the selected variants; \code{se_b_joint.sele}: The conditional 
#' standard error of the selected variants.
#' 
#' @author Zheng Ning, Xia Shen
#' 
#' @references 
#' Zheng Ning, Yakov Tsepilov, Sodbo Zh. Sharapov, Alexander K. Grishenko, Masoud Shirali, Peter K. Joshi,
#' James F. Wilson, Yudi Pawitan, Chris S. Haley, Yurii S. Aulchenko, Xia Shen (2018).
#' Multivariate discovery, replication, and interpretation of pleiotropic loci using summary association statistics. \emph{Submitted}.
#' 
#' @seealso 
#' \code{MultiSummary}
#' 
#' @examples 
#' \dontrun{
#' data(example.MultiSecondary)
#' ##### 474 snps around rs905938 #####
#' 
#' ## Six-traits cMVA ##
#' traits <- c("HEIGHT", "BMI", "HIP", "WC", "WHR", "WEIGHT")
#' MultiSecondary(gwa.region = example.gwas, LD.ref = example.LD, 
#' 				  snp.ref = example.snp.ref, R.ref = example.R.ref,
#'                p.threshold = 5e-8, tol = 0.9, traits = traits, T2.return = TRUE)
#' 
#' ## Three-traits cMVA ##
#' traits <- c("HEIGHT", "BMI", "HIP")
#' MultiSecondary(gwa.region = example.gwas, LD.ref = example.LD, 
#'                snp.ref = example.snp.ref, R.ref = example.R.ref,
#'                p.threshold = 5e-4, tol = 0.9, traits = traits, T2.return = TRUE)
#' }
#' @aliases MultiSecondary, multi.secondary
#' @keywords multivariate, meta-analysis, conditional analysis
#' 

`MultiSecondary` <- function(gwa.region, LD.ref, snp.ref, R.ref, p.threshold = 5e-8, tol = 0.8, traits, v_y = NULL, T2.return = FALSE){
  if(is.null(v_y) == TRUE){
    v_y <- rep(1,length(traits))
  }
  adj <- function(x,y){
    (-abs(x-y)+x+y)/2/x/y
  }
  
  snp.names <- intersect(rownames(LD.ref), rownames(gwa.region))
  idx.diff.ref <- which(gwa.region[snp.names,"A2"] != snp.ref[snp.names])
  diff.adj.vec <- 1 - 2*(gwa.region[snp.names,"A2"] != snp.ref[snp.names])
  names(diff.adj.vec) <- snp.names
  
  k <- length(traits)
  N.vec <- as.matrix(gwa.region[snp.names,paste0(traits, ".N")])
  for(i in 1:k){
    assign(paste0(traits[i],".N.adj"), outer(N.vec[,i],N.vec[,i],FUN = "adj"))
  }
  
  R.ref <- R.ref[traits, traits]
  
  snp.sele <- c()
  snp.left <- snp.names
  m <- length(snp.left)
  
  e.value <- svd(R.ref)$d
  e.vector <- svd(R.ref)$u
  
  b_uni <- as.matrix(gwa.region[snp.names, paste0(traits,".beta")])
  if(length(idx.diff.ref) > 0){
    b_uni[idx.diff.ref] <- -b_uni[idx.diff.ref]
  }
  
  se_uni <- as.matrix(gwa.region[snp.names, paste0(traits,".se")])
  
  t_uni <- b_uni / se_uni
  vr <- t_uni / sqrt(N.vec)
  vl <- sqrt(N.vec*se_uni^2)
  var_X <- 1 / (N.vec*se_uni^2)
  T2.sele <- p.sele <- b_joint.sele <- se_b_joint.sele <- NULL
  
  while(m > 0){
    t_left_mat <- matrix(0,m,k)
    for(j in 1:m){
      cor_X_inv <- solve(LD.ref[c(snp.sele, snp.left[j]), c(snp.sele, snp.left[j])])
      cor_X_next <- cor_X_inv[length(snp.sele) + 1,]
      se_b_joint_next <- numeric(k)
      for(i in 1:k){
        sd_X_next <- sqrt(var_X[c(snp.sele, snp.left[j]),i])  ## variance of X_c(snp.sele, snp.left[j]) for trait i
        cov_X_inv <- diag(1/sd_X_next, nrow = length(sd_X_next)) %*% cor_X_inv %*% diag(1/sd_X_next, nrow = length(sd_X_next))
        cov_X_adj <- (diag(sd_X_next, nrow = length(sd_X_next)) %*% 
                        LD.ref[c(snp.sele, snp.left[j]), c(snp.sele, snp.left[j])] %*% 
                        diag(sd_X_next, nrow = length(sd_X_next))) * 
          get(paste0(traits[i],".N.adj"))[c(snp.sele, snp.left[j]), c(snp.sele, snp.left[j])]
        se_b_joint_next[i] <- sqrt((v_y[i]*cov_X_inv %*% cov_X_adj %*% cov_X_inv)[length(snp.sele) + 1,length(snp.sele) + 1])
      }
      b_joint_next <- vl[snp.left[j],]*crossprod(cor_X_next, vr[c(snp.sele, snp.left[j]),])
      t_left_mat[j,] <- b_joint_next / se_b_joint_next    
    }
    # t_summary <- rowSums( (gwa.region[snp.left, paste0(traits,".t")] %*% e.vector %*% diag(sqrt(1 / e.value)))^2)
    t_summary <- rowSums((t_left_mat %*% e.vector %*% diag(sqrt(1 / e.value)))^2)
    p_sum_ref <- sapply(t_summary, pchisq, df = nrow(R.ref), lower.tail = F)
    if(min(p_sum_ref) < p.threshold){
      snp.sele.this.time <- snp.left[which.min(p_sum_ref)]
      snp.sele <- c(snp.sele, snp.sele.this.time)
      
      snp.r2 <- (LD.ref[snp.left,snp.sele.this.time])^2
      snp.ld <- snp.left[which(snp.r2 > tol)]
      snp.left <- setdiff(snp.left, c(snp.sele, snp.ld))
      
      while(length(snp.sele) > 0){
        cor_X_inv.sele <- solve(LD.ref[snp.sele, snp.sele])
        b_joint.sele <- se_b_joint.sele <- matrix(0, length(snp.sele), k) # each column is joint se for a trait
        for(i in 1:k){
          sd_X.sele <- sqrt(var_X[snp.sele,i])
          cov_X_inv.sele <- diag(1/sd_X.sele, nrow = length(sd_X.sele)) %*% cor_X_inv.sele %*% diag(1/sd_X.sele, nrow = length(sd_X.sele))
          cov_X_adj.sele <- (diag(sd_X.sele, nrow = length(sd_X.sele)) %*%
                               LD.ref[snp.sele, snp.sele] %*%
                               diag(sd_X.sele, nrow = length(sd_X.sele))) *
            get(paste0(traits[i],".N.adj"))[snp.sele, snp.sele]
          se_b_joint.sele[,i] <- sqrt(diag(v_y[i]*cov_X_inv.sele %*% cov_X_adj.sele %*% cov_X_inv.sele))
          b_joint.sele[,i] <- vl[snp.sele,i] * crossprod(cor_X_inv.sele, vr[snp.sele,i])
        }
        t_joint.sele <- b_joint.sele / se_b_joint.sele
        hotelling.t.sele <- rowSums((t_joint.sele %*% e.vector %*% diag(sqrt(1 / e.value)))^2)
        p.sele <- sapply(hotelling.t.sele, pchisq, df = nrow(R.ref), lower.tail = F)
        if(max(p.sele) > p.threshold){
          snp.sele <- snp.sele[ - which.max(p.sele)]
        }
        if(max(p.sele) < p.threshold){
          cor_X_inv.sele <- solve(LD.ref[snp.sele, snp.sele])
          
          snp.r2 <- rowSums((LD.ref[snp.left,snp.sele] %*% cor_X_inv.sele) * LD.ref[snp.left,snp.sele])
          snp.ld <- snp.left[which(snp.r2 > tol)]
          snp.left <- setdiff(snp.left, snp.ld)
          break
        }
      }
    }else{
      if(length(snp.sele) > 0){
        names(p.sele) <- snp.sele
        rownames(b_joint.sele) <- rownames(se_b_joint.sele) <- snp.sele
        colnames(b_joint.sele) <- colnames(se_b_joint.sele) <- traits
        b_joint.sele <- diag(diff.adj.vec[snp.sele], nrow = length(snp.sele)) %*% b_joint.sele
        rownames(b_joint.sele) <- snp.sele
      }
      if(T2.return == T){
        if(length(snp.sele) > 0){
          T2.sele <- qchisq(p.sele, df = length(traits), lower.tail = F)
        }
        return(list(T2.sele = T2.sele, p.sele = p.sele, b_joint.sele = b_joint.sele, se_b_joint.sele = se_b_joint.sele))
      }
    }
    m <- length(snp.left)
  }
  
  if(length(snp.sele) > 0){
    names(p.sele) <- snp.sele
    rownames(b_joint.sele) <- rownames(se_b_joint.sele) <- snp.sele
    colnames(b_joint.sele) <- colnames(se_b_joint.sele) <- traits
    b_joint.sele <- diag(diff.adj.vec[snp.sele], nrow = length(snp.sele)) %*% b_joint.sele
    rownames(b_joint.sele) <- snp.sele
  }
  if(T2.return == T){
    if(length(snp.sele) > 0){
      T2.sele <- qchisq(p.sele, df = length(traits), lower.tail = F)
    }
    return(list(T2.sele = T2.sele, p.sele = p.sele, b_joint.sele = b_joint.sele, se_b_joint.sele = se_b_joint.sele))
  }
  
  return(list(p.sele = p.sele, b_joint.sele = b_joint.sele, se_b_joint.sele = se_b_joint.sele))
}
