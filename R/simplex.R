#' @export
simplex <- function(m) {
  m <- label_simplex(m)

  m <- dual_simplex(m, step = 1)

  cat("Primal Simplex\n")
  primal_simplex(m, step = 1)
}

label_simplex <- function(m) {
  x_col <- seq_len(ncol(m) - 1)
  x_row <- seq_len(nrow(m) - 1)
  colnames(m) <- c(paste0("x", x_col), "b")
  rownames(m) <- c(paste0("x", rev(rev(x_col)[x_row])), "z")
  m
}

primal_simplex <- function(m, step) {
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
  # b_j / A_ij >= 0
  nn <- which(b / pc >= 0)
  if (!length(nn)) stop("Primal simplex: no solution.")
  quot <- b[nn] / pc[nn]
  # Take first row index when min is tied
  pr_nn_index <- which(quot == min(quot))[1]
  pr_index <- nn[pr_nn_index]

  # Pivot cell value
  pv <- m[pr_index, pc_index]
  cat("Pivot cell: m[", pr_index, ", ", pc_index, "]: ", pv, "\n\n", sep = "")

  m <- simplex_step(m, pr_index, pc_index)

  primal_simplex(m, step + 1)
}

dual_simplex <- function(m, step) {
  z <- z_vec(m)
  b <- b_vec(m)

  if (all(b >= 0)) return(m)

  # Only cat once
  if (step == 1) cat("Dual simplex\n")

  ## Select first minimum b as pivot row
  pr_index <- which(b == min(b))[1]
  pr <- p_row(m, pr_index)

  ## select pivot col
  neg <- which(pr < 0)
  if (!length(neg)) stop("Dual simplex: no solution.")
  quot <- z[neg] / pr[neg]
  pc_index <- which(quot == max(quot))[1]

  cat("Step:", step, "\n")
  print(m)

  # Pivot cell value
  pv <- m[pr_index, pc_index]
  cat("Pivot cell: m[", pr_index, ", ", pc_index, "]: ", pv, "\n\n", sep = "")

  m <- simplex_step(m, pr_index, pc_index)

  dual_simplex(m, step + 1)
}

simplex_step <- function(m, i, j) {
  ## Standardise pivot row
  m[i, ] <- m[i, ] / m[i, j]

  ## Generate basis vector in pivot col
  row_indices <- setdiff(1:nrow(m), i)
  for (k in row_indices) {
    # Subtract standardised pivot row from all other rows m[k, j] / 1
    # times
    m[k, ] <- m[k, ] - m[k, j] * m[i, ]
  }

  ## Adjust basis vector names
  rownames(m)[i] <- colnames(m)[j]

  m
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

p_row <- function(m, i) {
  m[i, -ncol(m)]
}
