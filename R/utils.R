# function to extract correlation from a var_comp output from a brms model
extract_corr <- function(x, idx) {
  x <- x$cor[idx, 1, idx]
  tibble(
    value = c(x),
    pos_y = c(row(x)),
    pos_x = c(col(x)),
    name1 = rep(rownames(x), times = ncol(x)),
    name2 = rep(colnames(x), each = nrow(x))
  )
}
