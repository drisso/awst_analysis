get_CPM <- function (gene_per_sample_matrix) 
{
  tmp <- 10^6/apply(gene_per_sample_matrix, 2, sum)
  invisible(t(gene_per_sample_matrix) * tmp)
}
###
get_TPM <- function(gene_per_sample_matrix, gene_length) {
  ans <- 1e3 * gene_per_sample_matrix/gene_length
  #head(ans[1,] * gene_legth[1])
  tmp <- 1e6/apply(ans, 2, sum, na.rm = TRUE)
  ans <- tmp * t(ans)
  invisible(ans)
}
###
get_FPKM <- function (gene_per_sample_matrix, gene_length) {
  tmp <- 1e+06/apply(gene_per_sample_matrix, 2, sum, na.rm = TRUE)
  ans <- 1000 * gene_per_sample_matrix/gene_length
  ans <- tmp * t(ans)
  invisible(ans)
}
###
get_Hart2013 <- function(sample_per_gene_TPM_matrix) {
  ans <- sample_per_gene_TPM_matrix
  for(i in 1:nrow(sample_per_gene_TPM_matrix)) {
    wd <- log2(1 + sample_per_gene_TPM_matrix[i,])
    wwhich <- which(wd > 0)
    d <- density(wd[wwhich], n = 1000)
    ccenter <- (d$x[d$x > 1])[which.max(d$y[d$x > 1])]
    those_greater_than_center <- wd[which(wd > ccenter)]
    U <- mean(those_greater_than_center)
    ssd <- (U - ccenter) * sqrt(pi/2)
    ans[i,] <- (wd - ccenter)/ssd
  }
  invisible(ans)
}
###
calinski <- function (hhc, dist = NULL, gMax = round(1 + 3.3 * log(length(hhc$order), 10))) {
  if (is.null(dist)) {
    dist <- cophenetic(hhc)
    attr(dist, "method") <- "cophenetic"
  } else if(attr(dist, "method") != "euclidean") dist <- cophenetic(hhc)
  
  dist <- as.matrix(dist)^2
  A <- -dist/2
  A_bar <- apply(A, 1, mean)
  totalSum <- sum(diag(A) - 2 * A_bar + mean(A))
  n <- length(hhc$order)
  ans <- rep(0, gMax)
  for (g in 2:gMax) {
    cclust <- cutree(hhc, k = g)
    withinSum <- 0
    for (k in 1:g) {
      if (sum(cclust == k) == 1) 
        next
      A <- as.matrix(-dist/2)[cclust == k, cclust == k]
      A_bar <- apply(A, 1, mean)
      withinSum <- withinSum + sum(diag(A) - 2 * A_bar + 
                                     mean(A))
    }
    betweenSum <- totalSum - withinSum
    betweenSum <- betweenSum/(g - 1)
    withinSum <- withinSum/(n - g)
    ans[g] <- betweenSum/withinSum
  }
  class(ans) <- "calinski"
  return(ans)
}
###
add_Calinski_curve <- function(height = 0.6, from = 5, to = 70, shift = 0.25) {
  obj <- aCalinski/max(aCalinski) * height
  G <- length(obj)
  ccol <- rep("black", G)
  for (g in 2:(G - 1)) {
    check <- obj[g - 1] < obj[g] & obj[g + 1] < obj[g]
    if (check) ccol[g] <- "red"
  }
  xx <- floor(seq(from = from, to = to, length.out = length(aCalinski)))
  nums <- paste(1:G); nums[1] <- ""
  obj <- obj + shift
  text(xx, obj, nums, col = ccol)
  lines(xx, obj, col = "gray30", lty = "longdash")
}
vCramer <- function (x, y) 
{
  ans <- chisq.test(x, y)
  ans$vCramer <- ans$statistic/sum(ans$observed)
  ans$vCramer <- sqrt(ans$vCramer/(min(dim(ans$observed) - 
                                         1)))
  class(ans) <- "vCramer"
  ans
}
print.vCramer <- function(obj) {
  cat(obj$method)
  str <- paste0("\nX-stat = ", round(obj$statistic, 4), " (vCramer = ", round(100*obj$vCramer, 2), "%), ")
  str <- paste0(str, "df = ", obj$parameter, ", p-value = ", format(obj$p.value, digits = 6))
  cat(str)
}

