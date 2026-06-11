## Common routines for non-convex isotonic projection-and-pooling simulations
## The scripts use a finite product partial order and a shared quartic component.

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Please install it.", call. = FALSE)
  }
}

require_pkg("Matrix")
require_pkg("data.table")

make_index <- function(K, L, Tlev) {
  array(seq_len(K * L * Tlev), dim = c(K, L, Tlev))
}

make_product_order_edges <- function(K, L, Tlev) {
  idx <- make_index(K, L, Tlev)
  edges <- matrix(NA_integer_, nrow = 0, ncol = 2)
  for (tt in seq_len(Tlev)) {
    if (K >= 2L) {
      for (k in seq_len(K - 1L)) for (m in seq_len(L)) {
        edges <- rbind(edges, c(idx[k, m, tt], idx[k + 1L, m, tt]))
      }
    }
    if (L >= 2L) {
      for (k in seq_len(K)) for (m in seq_len(L - 1L)) {
        edges <- rbind(edges, c(idx[k, m, tt], idx[k, m + 1L, tt]))
      }
    }
  }
  if (Tlev >= 2L) {
    for (tt in seq_len(Tlev - 1L)) for (k in seq_len(K)) for (m in seq_len(L)) {
      edges <- rbind(edges, c(idx[k, m, tt], idx[k, m, tt + 1L]))
    }
  }
  edges
}

make_constraint_matrix <- function(edges, n_nodes) {
  n_edges <- nrow(edges)
  ii <- rep(seq_len(n_edges), each = 2L)
  jj <- as.vector(t(edges))
  xx <- rep(c(1, -1), times = n_edges)
  Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n_edges, n_nodes))
}

solve_iso_gurobi <- function(g, weights, edges, threads = 1L) {
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    stop("gurobi package is not available.", call. = FALSE)
  }
  n <- length(g)
  A <- make_constraint_matrix(edges, n)
  model <- list(
    modelsense = "min",
    Q = Matrix::Diagonal(n = n, x = 2 * weights),
    obj = -2 * weights * g,
    A = A,
    sense = rep("<", nrow(edges)),
    rhs = rep(0, nrow(edges)),
    lb = rep(-Inf, n),
    ub = rep( Inf, n)
  )
  params <- list(OutputFlag = 0, Threads = threads)
  fit <- gurobi::gurobi(model, params)
  if (!fit$status %in% c("OPTIMAL", "SUBOPTIMAL")) {
    stop("Gurobi failed in isotonic projection. Status: ", fit$status)
  }
  as.numeric(fit$x)
}

solve_iso_osqp <- function(g, weights, edges) {
  if (!requireNamespace("osqp", quietly = TRUE)) {
    stop("Neither gurobi nor osqp is available. Install osqp or gurobi.", call. = FALSE)
  }
  n <- length(g)
  P <- Matrix::Diagonal(n = n, x = 2 * weights)
  q <- -2 * weights * g
  A <- make_constraint_matrix(edges, n)
  l <- rep(-Inf, nrow(A))
  u <- rep(0, nrow(A))
  settings <- if ("osqpSettings" %in% getNamespaceExports("osqp")) {
    osqp::osqpSettings(verbose = FALSE, eps_abs = 1e-8, eps_rel = 1e-8, max_iter = 200000)
  } else {
    list(verbose = FALSE, eps_abs = 1e-8, eps_rel = 1e-8, max_iter = 200000)
  }
  fit <- osqp::osqp(P = P, q = q, A = A, l = l, u = u, pars = settings)
  res <- fit$Solve()
  if (is.null(res$x)) stop("OSQP failed to return a solution.")
  as.numeric(res$x)
}

solve_isotonic_projection <- function(g, weights, edges, solver = c("auto", "gurobi", "osqp"), threads = 1L) {
  solver <- match.arg(solver)
  if (solver == "gurobi") return(solve_iso_gurobi(g, weights, edges, threads = threads))
  if (solver == "osqp") return(solve_iso_osqp(g, weights, edges))
  if (requireNamespace("gurobi", quietly = TRUE)) {
    return(solve_iso_gurobi(g, weights, edges, threads = threads))
  }
  solve_iso_osqp(g, weights, edges)
}

quartic_transform_one <- function(u, lambda, tol = 1e-10) {
  ## Maximize u*x - 0.5*x^2 + lambda*x^4 over [0,1].
  ## We use the maximal selection if the scalar maximizer is not unique.
  if (abs(lambda) < .Machine$double.eps) return(min(max(u, 0), 1))
  roots <- polyroot(c(u, -1, 0, 4 * lambda))
  real_roots <- Re(roots[abs(Im(roots)) < 1e-8])
  cand <- c(0, 1, real_roots[real_roots >= -1e-10 & real_roots <= 1 + 1e-10])
  cand <- unique(round(pmin(pmax(cand, 0), 1), 14))
  vals <- u * cand - 0.5 * cand^2 + lambda * cand^4
  max_val <- max(vals)
  max(cand[vals >= max_val - tol])
}

quartic_transform <- function(u, lambda) {
  vapply(u, quartic_transform_one, numeric(1), lambda = lambda)
}

make_true_surface <- function(K, L, alphas, k0, m0) {
  Tlev <- length(alphas)
  if (Tlev != 2L) stop("This DGP assumes two quantile levels.")
  x <- array(0, dim = c(K, L, Tlev))
  d_sigma <- c(0.50, 0.90)
  d_ell   <- c(0.35, 0.60)
  d_int   <- c(0.40, 0.70)
  for (tt in seq_len(Tlev)) {
    for (k in seq_len(K)) for (m in seq_len(L)) {
      x[k, m, tt] <- 0.35 * k / K + 0.25 * m / L +
        d_sigma[tt] * as.numeric(k >= k0) +
        d_ell[tt]   * as.numeric(m >= m0) +
        d_int[tt]   * as.numeric(k >= k0 && m >= m0)
    }
  }
  x
}

calibrate_t_location_scale <- function(x_star, alphas, nu) {
  K <- dim(x_star)[1]; L <- dim(x_star)[2]
  p1 <- 1 - alphas[1]
  p2 <- 1 - alphas[2]
  z1 <- stats::qt(p1, df = nu)
  z2 <- stats::qt(p2, df = nu)
  mu <- sigma <- matrix(0, K, L)
  for (k in seq_len(K)) for (m in seq_len(L)) {
    target1 <- -x_star[k, m, 1]
    target2 <- -x_star[k, m, 2]
    sigma[k, m] <- (target1 - target2) / (z1 - z2)
    mu[k, m] <- target1 - sigma[k, m] * z1
  }
  list(mu = mu, sigma = sigma)
}

make_pi_matrix <- function(K, L, k0, m0, pi_base, pi_near, pi_deep) {
  pi <- matrix(pi_base, K, L)
  for (k in seq_len(K)) for (m in seq_len(L)) {
    near_boundary <- (k >= k0 - 1L && m >= m0 - 1L)
    deep_stress <- (k >= k0 && m >= m0)
    pi[k, m] <- pi[k, m] + pi_near * near_boundary + pi_deep * deep_stress
  }
  pmin(pi, 0.25)
}

generate_returns <- function(mu, sigma, N, nu, perturb = FALSE, pi_mat = NULL, shock_size = 0) {
  K <- nrow(mu); L <- ncol(mu)
  R <- array(0, dim = c(K, L, N))
  for (k in seq_len(K)) for (m in seq_len(L)) {
    r <- mu[k, m] + sigma[k, m] * stats::rt(N, df = nu)
    if (perturb) {
      xi <- stats::rbinom(N, size = 1L, prob = pi_mat[k, m])
      r <- r - shock_size * xi
    }
    R[k, m, ] <- r
  }
  R
}

make_empirical_loss_quantiles <- function(R, alphas) {
  K <- dim(R)[1]; L <- dim(R)[2]; Tlev <- length(alphas)
  q <- array(0, dim = c(K, L, Tlev))
  Loss <- -R
  for (k in seq_len(K)) for (m in seq_len(L)) {
    for (tt in seq_len(Tlev)) {
      q[k, m, tt] <- as.numeric(stats::quantile(Loss[k, m, ], probs = alphas[tt], type = 8, names = FALSE))
    }
  }
  if (Tlev >= 2L) {
    for (tt in 2:Tlev) q[, , tt] <- pmax(q[, , tt], q[, , tt - 1L])
  }
  q
}

array_to_vector <- function(x) as.numeric(x)
vector_to_array <- function(v, K, L, Tlev) array(v, dim = c(K, L, Tlev))

pinball_loss <- function(loss, q, alpha) {
  e <- loss - q
  (alpha - as.numeric(e < 0)) * e
}

compute_quantile_loss <- function(R_test, x_hat, alpha, tt, stress_only = FALSE, k0 = NULL, m0 = NULL) {
  K <- dim(R_test)[1]; L <- dim(R_test)[2]
  vals <- numeric(0)
  for (k in seq_len(K)) for (m in seq_len(L)) {
    if (stress_only && !(k >= k0 && m >= m0)) next
    loss <- -R_test[k, m, ]
    vals <- c(vals, pinball_loss(loss, x_hat[k, m, tt], alpha))
  }
  mean(vals)
}

jump_components <- function(x, k0, m0, tt) {
  K <- dim(x)[1]; L <- dim(x)[2]
  left_cols <- seq_len(max(m0 - 1L, 1L))
  low_rows <- seq_len(max(k0 - 1L, 1L))
  high_rows <- k0:K
  high_cols <- m0:L
  base <- mean(x[low_rows, left_cols, tt])
  vol <- mean(x[high_rows, left_cols, tt]) - base
  lev <- mean(x[low_rows, high_cols, tt]) - base
  inter <- mean(x[high_rows, high_cols, tt]) - base - vol - lev
  c(vol = vol, drawdown = lev, interaction = inter)
}

jump_error <- function(x_hat, x_true, k0, m0, tt) {
  mean(abs(jump_components(x_hat, k0, m0, tt) - jump_components(x_true, k0, m0, tt)))
}

evaluate_fit <- function(x_hat, x_true, R_test, alphas, k0, m0, method, lambda) {
  K <- dim(x_true)[1]; L <- dim(x_true)[2]
  stress <- array(FALSE, dim = dim(x_true))
  stress[k0:K, m0:L, ] <- TRUE
  e2 <- (x_hat - x_true)^2
  out <- data.table::data.table(
    method = method,
    lambda = lambda,
    grid_mse = mean(e2),
    stress_mse = mean(e2[stress]),
    jump_err_095 = jump_error(x_hat, x_true, k0, m0, 1L),
    jump_err_099 = jump_error(x_hat, x_true, k0, m0, 2L)
  )
  for (tt in seq_along(alphas)) {
    nm <- sprintf("ql_%03d", as.integer(round(1000 * alphas[tt])))
    out[[nm]] <- compute_quantile_loss(R_test, x_hat, alphas[tt], tt)
  }
  out
}

validation_stress_mse <- function(x_hat, q_val, k0, m0) {
  K <- dim(x_hat)[1]; L <- dim(x_hat)[2]
  mean((x_hat[k0:K, m0:L, ] - q_val[k0:K, m0:L, ])^2)
}

select_lambda <- function(g_proj_scaled, x_upper, q_val, K, L, Tlev, lambdas, k0, m0,
                          boundary_tol = 1e-8, tol_frac = 0.005, min_gain_frac = 0.005) {
  cand <- list()
  for (lam in lambdas) {
    z <- quartic_transform(g_proj_scaled, lam)
    x <- vector_to_array(z * x_upper, K, L, Tlev)
    boundary <- any(z <= boundary_tol | z >= 1 - boundary_tol)
    score <- validation_stress_mse(x, q_val, k0, m0)
    cand[[length(cand) + 1L]] <- data.table::data.table(lambda = lam, score = score, boundary = boundary)
  }
  tab <- data.table::rbindlist(cand)
  convex_score <- tab[lambda == 0, score][1]
  eligible <- tab[lambda == 0 | !boundary]
  best <- min(eligible$score, na.rm = TRUE)
  selected <- eligible[score <= best * (1 + tol_frac)]
  selected_lambda <- min(selected$lambda)
  selected_score <- tab[lambda == selected_lambda, score][1]
  if (selected_lambda > 0 && (convex_score - selected_score) / convex_score < min_gain_frac) {
    selected_lambda <- 0
  }
  list(lambda = selected_lambda, table = tab)
}

simulate_quartic_replication <- function(rep_id, cfg) {
  set.seed(cfg$base_seed + rep_id)
  K <- cfg$K; L <- cfg$L; alphas <- cfg$alphas; Tlev <- length(alphas)
  x_true <- make_true_surface(K, L, alphas, cfg$k0, cfg$m0)
  par <- calibrate_t_location_scale(x_true, alphas, cfg$nu_t)
  pi_mat <- make_pi_matrix(K, L, cfg$k0, cfg$m0, cfg$pi_base, cfg$pi_near, cfg$pi_deep)
  R_train <- generate_returns(par$mu, par$sigma, cfg$n_train, cfg$nu_t, TRUE, pi_mat, cfg$shock_size)
  R_val   <- generate_returns(par$mu, par$sigma, cfg$n_val,   cfg$nu_t, FALSE)
  R_test  <- generate_returns(par$mu, par$sigma, cfg$n_test,  cfg$nu_t, FALSE)
  q_train <- make_empirical_loss_quantiles(R_train, alphas)
  q_val <- make_empirical_loss_quantiles(R_val, alphas)

  x_upper <- cfg$x_upper_mult * max(x_true)
  g_scaled <- array_to_vector(q_train / x_upper)
  weights <- rep(cfg$n_train, length(g_scaled))
  edges <- make_product_order_edges(K, L, Tlev)
  g_proj <- solve_isotonic_projection(g_scaled, weights, edges,
                                      solver = cfg$solver, threads = cfg$threads_per_worker)

  sel <- select_lambda(g_proj, x_upper, q_val, K, L, Tlev, cfg$lambda_grid,
                       cfg$k0, cfg$m0, cfg$boundary_tol, cfg$validation_tol_frac,
                       cfg$min_validation_gain_frac)

  methods <- list(
    convex = 0,
    quartic_fixed = cfg$fixed_lambda,
    quartic_selected = sel$lambda
  )
  out <- list()
  for (nm in names(methods)) {
    lam <- methods[[nm]]
    z <- quartic_transform(g_proj, lam)
    x_hat <- vector_to_array(z * x_upper, K, L, Tlev)
    out[[nm]] <- evaluate_fit(x_hat, x_true, R_test, alphas, cfg$k0, cfg$m0, nm, lam)
  }
  ans <- data.table::rbindlist(out, use.names = TRUE)
  ans[, rep := rep_id]
  ans[, selected_lambda := sel$lambda]
  ans
}

summarize_simulation <- function(dt) {
  metric_cols <- setdiff(names(dt), c("method", "lambda", "rep", "selected_lambda"))
  out <- dt[, lapply(.SD, function(x) c(mean = mean(x), sd = stats::sd(x))), by = method, .SDcols = metric_cols]
  out
}

gain_summary <- function(dt, base_method = "convex") {
  metric_cols <- setdiff(names(dt), c("method", "lambda", "rep", "selected_lambda"))
  wide <- data.table::dcast(dt, rep ~ method, value.var = metric_cols)
  rows <- list(); rr <- 1L
  for (metric in metric_cols) {
    base_col <- paste(metric, base_method, sep = "_")
    for (method in setdiff(unique(dt$method), base_method)) {
      comp_col <- paste(metric, method, sep = "_")
      gain <- 100 * (wide[[base_col]] - wide[[comp_col]]) / wide[[base_col]]
      rows[[rr]] <- data.table::data.table(
        method = method,
        metric = metric,
        mean_gain_pct = mean(gain, na.rm = TRUE),
        se_gain_pct = stats::sd(gain, na.rm = TRUE) / sqrt(sum(is.finite(gain)))
      )
      rr <- rr + 1L
    }
  }
  data.table::rbindlist(rows)
}
