#' @export
simplex <- function(m) {
  m <- label_simplex(m)

  do_simplex(m, step = 1)
}

label_simplex <- function(m) {
  x_col <- seq_len(ncol(m) - 1)
  x_row <- seq_len(nrow(m) - 1)
  colnames(m) <- c(paste0("x", x_col), "b")
  rownames(m) <- c(paste0("x", rev(rev(x_col)[x_row])), "z")
  m
}

do_simplex <- function(m, step) {
  z <- z_vec(m)
  b <- b_vec(m)

  # Simplex finished
  if (all(z >= 0)) {
    cat("Simplex finished\n")
    return(m)
  }

  cat("Step:", step, "\n")
  print(m)

  ## Select maximum negative z as pivot column (take first when tied)
  pc_index <- which(z == min(z[z < 0]))[1]
  pc <- p_col(m, pc_index)

  ## Select pivot row
  nn <- which(b / pc >= 0)
  if (!length(nn)) stop("Simplex failed.")
  quot <- b[nn] / pc[nn]
  pr_nn_index <- which(quot == min(quot))
  pr_index <- nn[pr_nn_index]

  pv <- m[pr_index, pc_index]

  cat("Pivot cell: m[", pr_index, ", ", pc_index, "]: ", pv, "\n\n", sep = "")

  # Standardise pivot row
  m[pr_index, ] <- m[pr_index, ] / pv

  # Generate basis vector in pivot col
  other_row_indices <- setdiff(1:nrow(m), pr_index)
  for (i in other_row_indices) {
    m[i, ] <- m[i, ] - m[i, pc_index] * m[pr_index, ]
  }

  # Adjust basis vector names
  rownames(m)[pr_index] <- colnames(m)[pc_index]

  do_simplex(m, step + 1)
}

z_vec <- function(m) {
  m[nrow(m), -ncol(m)]
}

b_vec <- function(m) {
  m[-nrow(m), ncol(m)]
}

p_col <- function(m, i) {
  m[-nrow(m),i]
}
