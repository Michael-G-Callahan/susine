#' Compute X %*% b using precomputed scaling
#'
#' @param X n x p (unstandardized) matrix with attributes `"scaled:center"`,
#'   `"scaled:scale"` and `"d"` set (see `initialize_X`).
#' @param b Length-p coefficient vector.
#'
#' @return Length-n numeric vector equal to X %*% b on the original scale.
#' @importFrom Matrix t
#' @importFrom Matrix tcrossprod
compute_Xb = function (X, b) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")

  # Scale Xb (X is an ordinary sparse/dense matrix).
  scaled.Xb = tcrossprod(X,t(b/csd))

  # Center Xb.
  Xb = scaled.Xb - sum(cm*b/csd)
  return(as.numeric(Xb))
}

#' Compute t(X) %*% y using precomputed scaling
#'
#' @param X n x p (unstandardized) matrix with attributes as in `compute_Xb`.
#' @param y Length-n response vector.
#'
#' @return Length-p numeric vector equal to t(X) %*% y on the standardized scale.
#' @importFrom Matrix t
#' @importFrom Matrix crossprod
compute_Xty = function (X, y) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")
  ytX = crossprod(y,X)

  # Scale Xty (X is an ordinary sparse/dense matrix).
  scaled.Xty = t(ytX/csd)

  # Center Xty.
  centered.scaled.Xty = scaled.Xty - cm/csd * sum(y)
  return(as.numeric(centered.scaled.Xty))
}

#' Compute M %*% t(X) using precomputed scaling
#'
#' @param M L x p matrix (e.g., stacked `b_hat` rows) to multiply by t(X).
#' @param X n x p (unstandardized) matrix with standardization attributes.
#'
#' @return L x n numeric matrix equal to M %*% t(X) on the original scale.
#' @importFrom Matrix t
compute_MXt = function (M, X) {
  cm = attr(X,"scaled:center")
  csd = attr(X,"scaled:scale")

  if (!is.null(attr(X,"matrix.type")))

    # When X is a trend filtering matrix.
    return(as.matrix(t(apply(M,1,function(b) compute_Xb(X,b)))))
  else

    # When X is an ordinary sparse/dense matrix.
    return(as.matrix(t(X %*% (t(M)/csd)) - drop(M %*% (cm/csd))))

  # This should be the same as
  #
  #   t(apply(M, 1, function(b) compute_Xb(X, b))))
  #
  # as well as
  #
  #   M %*% (t(X)/csd) - drop(tcrossprod(M,t(cm/csd)))
  #
  # but should be more memory-efficient.
}

#' Compute column stats and quadratic terms
#'
#' @description Computes column means and sds of X, and d = rowSums(Y^2),
#'   where Y is the centered/scaled version of X as requested.
#'
#' @param X n x p matrix (dense or sparse).
#' @param center Logical; center columns.
#' @param scale Logical; scale columns to unit variance (sd = 1 when variance = 0).
#'
#' @return List with `cm` (length p), `csd` (length p), and `d` (length p).
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
compute_colstats = function (X, center = TRUE, scale = TRUE) {
  n = nrow(X)
  p = ncol(X)

  # X is an ordinary dense or sparse matrix. Set sd = 1 when the
  # column has variance 0.
  if (center)
    cm = colMeans(X,na.rm = TRUE)
  else
    cm = rep(0,p)
  if (scale) {
    csd = compute_colSds(X)
    csd[csd == 0] = 1
  } else
    csd = rep(1,p)

  # These two lines of code should give the same result as
  #
  #   Y = (t(X) - cm)/csd
  #   d = rowSums(Y^2)
  #
  # for all four combinations of "center" and "scale", but do so
  # without having to modify X, or create copies of X in memory. In
  # particular the first line should be equivalent to colSums(X^2).
  d = n*colMeans(X)^2 + (n-1)*compute_colSds(X)^2
  d = (d - n*cm^2)/csd^2

  return(list(cm = cm,csd = csd,d = d))
}

#' Compute column standard deviations (dense or sparse)
#'
#' @description Equivalent to `matrixStats::colSds` for dense matrices, but
#'   also supports sparse matrices.
#'
#' @param X n x p matrix (dense or sparse).
#'
#' @return Length-p vector of column standard deviations.
#' @importFrom matrixStats colSds
#' @importFrom Matrix summary
compute_colSds = function(X) {
  return(colSds(X))
}
