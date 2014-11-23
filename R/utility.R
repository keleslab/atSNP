.structure_diff <- function(pval_mat) {
  id <- apply(pval_mat[, c(2, 4)], 1, which.min)
  return(cbind(pval_mat[, c(1, 3)][cbind(seq_along(id), id)],
               pval_mat[, c(2, 4)][cbind(seq_along(id), id)]))
}

.structure <- function(pval_mat) {
  id1 <- apply(pval_mat[, c(2, 4)], 1, which.min)
  id2 <- apply(pval_mat[, c(6, 8)], 1, which.min)
  return(cbind(
               pval_mat[, c(1, 3)][cbind(seq_along(id1), id1)],
               pval_mat[, c(5, 7)][cbind(seq_along(id2), id2)],
               pval_mat[, c(2, 4)][cbind(seq_along(id1), id1)],
               pval_mat[, c(6, 8)][cbind(seq_along(id2), id2)])
         )
}

CheckSameLength <- function(x) {
  if(length(x) == 1) {
    return(TRUE)
  }
  return(var(unlist(sapply(x, length))) == 0)
}