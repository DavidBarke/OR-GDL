#' @export
bellmann_ford <- function(B) {
  Tree <- init_tree(B)

  B_list <- structure(
    list(
      list(
        B,
        Tree
      )
    ),
    class = "bellmann_ford"
  )

  while (TRUE) {
    B_list[[length(B_list) + 1]] <- bellmann_ford_product(
      U = B_list[[length(B_list)]][[1]],
      B = B,
      Tree = B_list[[length(B_list)]][[2]]
    )

    if (
      all(
        B_list[[length(B_list)]][[1]] == B_list[[length(B_list) - 1]][[1]]
      )
    ) break()

    if (any(diag(B_list[[length(B_list)]][[1]]) < 0)) {
      print(B_list)
      stop("Negative cycle detected!", call. = FALSE)
    }

    if (length(B_list) == nrow(B)) break()
  }

  B_list
}

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

init_tree <- function(B) {
  Tree <- B
  for (i in 1:nrow(B)) {
    for (j in 1:ncol(B)) {
      val <- Tree[i, j]
      if (val == 0 || !is.finite(val)) {
        Tree[i, j] <- -1
      } else {
        Tree[i, j] <- i
      }
    }
  }
  Tree
}

#' @export
print.bellmann_ford <- function(x, ...) {
  cat("Bellmann-Ford-Algorithm\n")
  for (i in seq_along(x)) {
    cat("Step", i, "\n")
    cat("U:\n")
    print(x[[i]][[1]])
    cat("Tree:\n")
    print(x[[i]][[2]])
    cat("\n")
  }
}
