bellmann_ford_product <- function(U, B, Tree) {
  stopifnot(all.equal(dim(U), dim(B)))
  stopifnot(nrow(U) == ncol(U))

  n <- nrow(U)

  # V = U x B
  V <- diag(Inf, nrow = n)
  for (i in 1:n) {
    for (j in 1:n) {
      sums <- U[i,] + B[,j]
      m <- min(sums)
      min_index <- which(sums == m)[1]

      if (is.finite(m) && (i != j) && m < U[i, j]) {
        Tree[i, j] <- min_index
      }

      V[i, j] <- m
    }
  }

  list(V, Tree)
}
