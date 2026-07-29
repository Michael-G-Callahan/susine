# mini_sim_prior_methods.R  —  v2
# Expanded factorial benchmark: real M1-stratified genotype matrices + varied
# p*, annotation quality (r2), and noise fraction.
# Run from this directory: Rscript mini_sim_prior_methods.R
#
# Factorial: 3 X sources x 3 p_star x 2 annotation_r2 x 2 noise x 1 seed = 36 runs/method
#
# X sources:
#   N3           — SuSiE_N3_X  (574 x 1001, standard benchmark matrix)
#   low_m1_KLK11 — KLK11 locus  (600 x  997, M1 = 0.084, simple LD)
#   high_m1_AMT  — AMT locus    (600 x  974, M1 = 0.229, complex LD)

# Self-locate: set working directory to this script's own directory so that
# relative paths (../../, ../../../) resolve correctly regardless of how Rscript
# is invoked (absolute path vs. from within the directory).
local({
  args <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", args, value = TRUE)
  if (length(file_flag) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_flag[1])))
    setwd(script_dir)
  }
})

devtools::load_all("../../", quiet = TRUE)

# ---------------------------------------------------------------------------
# 1. Load X matrices
# ---------------------------------------------------------------------------
cat("Loading X matrices...\n")

# 1a. SuSiE_N3_X (bundled in test_susine package data)
n3_rda_path <- "../../../test_susine/data/SuSiE_N3_X.rda"
if (!file.exists(n3_rda_path)) stop("Cannot find SuSiE_N3_X.rda at: ", n3_rda_path)
load(n3_rda_path)   # loads object SuSiE_N3_X into environment

# 1b. Real genotype matrices (M1-stratified, scenario_1, closest-1000 SNPs, n=600)
data_root <- "../../../test_susine/data/sampled_simulated_genotypes/scenario_1/"

load_tsv_gz <- function(path) {
  if (!file.exists(path)) stop("Matrix file not found: ", path)
  mat <- readr::read_tsv(
    path,
    col_types      = readr::cols(.default = "d", participant_ID = "c"),
    show_col_types = FALSE,
    progress       = FALSE
  )
  m <- as.matrix(mat[, -1])   # drop participant_ID column
  storage.mode(m) <- "double"
  m
}

X_low_m1 <- load_tsv_gz(paste0(
  data_root, "KLK11_ENST00000453757_3_chr19_bp51529870_closest1000_n600.tsv.gz"))
X_high_m1 <- load_tsv_gz(paste0(
  data_root, "AMT_ENST00000495436_1_chr3_bp49459884_closest1000_n600.tsv.gz"))

x_list <- list(
  N3           = as.matrix(SuSiE_N3_X),   # ensure matrix, not data.frame
  low_m1_KLK11 = X_low_m1,
  high_m1_AMT  = X_high_m1
)
cat(sprintf("  N3:           %d x %d\n",   nrow(SuSiE_N3_X),  ncol(SuSiE_N3_X)))
cat(sprintf("  low_m1_KLK11: %d x %d  (M1 = 0.084)\n", nrow(X_low_m1),  ncol(X_low_m1)))
cat(sprintf("  high_m1_AMT:  %d x %d  (M1 = 0.229)\n", nrow(X_high_m1), ncol(X_high_m1)))

# ---------------------------------------------------------------------------
# 2. Simulation helpers (no test_susine dependency)
# ---------------------------------------------------------------------------

# Annotation vector: cor(a[causal], beta[causal]) ~= sqrt(annotation_r2).
# Non-causal SNPs get a = 0 (no spurious signal).
make_annotation <- function(beta, annotation_r2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  causal <- which(beta != 0)
  a <- rep(0, length(beta))
  if (length(causal) == 0 || annotation_r2 <= 0) return(a)
  b_c <- beta[causal]
  # abs() for single causal SNP where sd() would return NA
  b_sd <- if (length(b_c) > 1) sd(b_c) else abs(b_c)
  if (b_sd == 0) b_sd <- 1
  noise_sd  <- b_sd * sqrt((1 - annotation_r2) / annotation_r2)
  a[causal] <- b_c + rnorm(length(b_c), 0, noise_sd)
  a
}

# Phenotype: noise_fraction = proportion of total variance attributed to noise.
make_phenotype <- function(X, beta, noise_fraction, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Xb         <- drop(X %*% beta)
  signal_var <- mean(Xb^2)
  if (signal_var == 0 || noise_fraction >= 1) {
    sigma2 <- 1
  } else {
    sigma2 <- signal_var * noise_fraction / (1 - noise_fraction)
  }
  y <- Xb + rnorm(nrow(X), 0, sqrt(sigma2))
  y - mean(y)
}

# AUPRC via trapezoidal integration
compute_auprc <- function(pips, is_causal) {
  ord          <- order(pips, decreasing = TRUE)
  truth_sorted <- is_causal[ord]
  tp           <- cumsum(truth_sorted)
  fp           <- cumsum(!truth_sorted)
  precision    <- tp / (tp + fp)
  recall       <- tp / max(sum(is_causal), 1L)
  n_pts        <- length(recall)
  if (n_pts < 2) return(NA_real_)
  sum(diff(recall) * (precision[-1] + precision[-n_pts]) / 2)
}

# ---------------------------------------------------------------------------
# 3. Experimental grid
# ---------------------------------------------------------------------------
p_star_grid        <- c(1, 3, 10)
annotation_r2_grid <- c(0.2, 0.8)
noise_grid         <- c(0.5, 0.8)
n_seeds            <- 1
L                  <- 5

methods <- c("none", "var", "scale_var", "scale")

annotation_methods <- c("scale_var", "scale")

total_runs <- length(x_list) * length(p_star_grid) *
              length(annotation_r2_grid) * length(noise_grid) * n_seeds
cat(sprintf(
  "\nGrid: %d X sources x %d p_star x %d r2 x %d noise x %d seeds = %d runs/method\n",
  length(x_list), length(p_star_grid), length(annotation_r2_grid),
  length(noise_grid), n_seeds, total_runs))
cat(sprintf("Total fits: %d runs x %d methods = %d\n\n",
            total_runs, length(methods), total_runs * length(methods)))

# ---------------------------------------------------------------------------
# 4. Result storage
# ---------------------------------------------------------------------------
run_results <- data.frame(
  run_id        = integer(),
  x_source      = character(),
  p_star        = integer(),
  annotation_r2 = numeric(),
  noise         = numeric(),
  seed          = integer(),
  method        = character(),
  auprc         = numeric(),
  final_elbo    = numeric(),
  n_iter        = integer(),
  error         = character(),
  stringsAsFactors = FALSE
)

# mu_eff_l = sum(b_hat[l, ]) — expected total contribution of effect l to fitted values
effect_results <- data.frame(
  run_id        = integer(),
  x_source      = character(),
  p_star        = integer(),
  annotation_r2 = numeric(),
  noise         = numeric(),
  seed          = integer(),
  method        = character(),
  effect        = integer(),
  c_l           = numeric(),
  sigma_0_2_l   = numeric(),
  mu_eff_l      = numeric(),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 5. Main loop
# ---------------------------------------------------------------------------
run_id <- 0
for (x_name in names(x_list)) {
  X_full <- x_list[[x_name]]
  p_full <- ncol(X_full)

  for (p_star in p_star_grid) {
    for (r2 in annotation_r2_grid) {
      for (noise_frac in noise_grid) {
        for (s in seq_len(n_seeds)) {
          run_id    <- run_id + 1
          base_seed <- run_id * 17 + s

          # Generate data
          set.seed(base_seed)
          causal       <- sample(p_full, min(p_star, p_full))
          beta         <- rep(0, p_full)
          beta[causal] <- rnorm(p_star, 0, 0.8)
          is_causal    <- (beta != 0)

          y <- make_phenotype(X_full, beta, noise_frac, seed = base_seed + 1)
          a <- make_annotation(beta, r2,        seed = base_seed + 2)

          for (m in methods) {
            mu_0_input <- if (m %in% annotation_methods) a else NULL

            fit <- tryCatch(
              susine(
                L                   = L,
                X                   = X_full,
                y                   = y,
                mu_0                = mu_0_input,
                sigma_0_2           = 0.2,
                prior_update_method = m,
                ig_alpha            = 2,
                ig_beta             = 0.1,
                max_iter            = 50,
                tol                 = 1e-4
              ),
              error = function(e) list(error = conditionMessage(e))
            )

            if (!is.null(fit$error)) {
              run_results <- rbind(run_results, data.frame(
                run_id = run_id, x_source = x_name, p_star = p_star,
                annotation_r2 = r2, noise = noise_frac, seed = s, method = m,
                auprc = NA, final_elbo = NA, n_iter = NA,
                error = fit$error, stringsAsFactors = FALSE
              ))
            } else {
              run_results <- rbind(run_results, data.frame(
                run_id        = run_id,
                x_source      = x_name,
                p_star        = p_star,
                annotation_r2 = r2,
                noise         = noise_frac,
                seed          = s,
                method        = m,
                auprc         = compute_auprc(fit$model_fit$PIPs, is_causal),
                final_elbo    = tail(fit$model_fit$elbo, 1),
                n_iter        = length(fit$model_fit$elbo) - 1L,
                error         = "",
                stringsAsFactors = FALSE
              ))

              # Per-effect storage
              # After finalize_effect_fits(): alpha and b_hat are L x p matrices
              # After finalize_priors(): sigma_0_2 is L x p, mu_0_scale_factor is length L
              alpha_mat <- fit$effect_fits$alpha         # L x p
              b_hat_mat <- fit$effect_fits$b_hat         # L x p
              sigma_mat <- fit$priors$sigma_0_2          # L x p
              c_vec     <- fit$priors$mu_0_scale_factor  # length L

              for (l in seq_len(L)) {
                mu_eff_l    <- sum(b_hat_mat[l, ])
                c_l         <- if (!is.null(c_vec) && length(c_vec) >= l) c_vec[l] else 1
                sigma_0_2_l <- sigma_mat[l, 1]

                effect_results <- rbind(effect_results, data.frame(
                  run_id        = run_id,
                  x_source      = x_name,
                  p_star        = p_star,
                  annotation_r2 = r2,
                  noise         = noise_frac,
                  seed          = s,
                  method        = m,
                  effect        = l,
                  c_l           = c_l,
                  sigma_0_2_l   = sigma_0_2_l,
                  mu_eff_l      = mu_eff_l,
                  stringsAsFactors = FALSE
                ))
              }
            }
          }
        }
      }
    }
  }
  cat(sprintf("  Completed X source: %-16s  (run_id = %d)\n", x_name, run_id))
}

# ===========================================================================
# 6. Summary output
# ===========================================================================
ok  <- run_results[run_results$error == "", ]
err <- run_results[run_results$error != "", ]

cat("\n")
cat("=======================================================================\n")
cat(sprintf("SECTION A -- Overall AUPRC per method (all %d runs)\n", total_runs))
cat("=======================================================================\n")
cat(sprintf("%-22s  %3s  %6s  %6s  %8s  %5s\n",
            "method", "n", "AUPRC", "sd", "ELBO", "iter"))
cat(strrep("-", 60), "\n")
for (m in methods) {
  sub <- ok[ok$method == m, ]
  cat(sprintf("%-22s  %3d  %6.3f  %6.3f  %8.1f  %5.1f\n",
              m, nrow(sub),
              mean(sub$auprc,      na.rm = TRUE),
              sd(sub$auprc,        na.rm = TRUE),
              mean(sub$final_elbo, na.rm = TRUE),
              mean(sub$n_iter,     na.rm = TRUE)))
}

cat("\n")
cat("=======================================================================\n")
cat("SECTION B -- AUPRC by method x scenario axis\n")
cat("=======================================================================\n")

print_axis_table <- function(col, vals, label, width = 7) {
  cat(sprintf("\n  --- By %s ---\n", label))
  hdr <- sprintf("%-22s", "method")
  for (v in vals) hdr <- paste0(hdr, sprintf(paste0("  %", width, "s"), as.character(v)))
  cat(hdr, "\n")
  cat(strrep("-", nchar(hdr)), "\n")
  for (m in methods) {
    row <- sprintf("%-22s", m)
    for (v in vals) {
      sub <- ok[ok$method == m & ok[[col]] == v, ]
      row <- paste0(row, sprintf(paste0("  %", width, ".3f"),
                                 mean(sub$auprc, na.rm = TRUE)))
    }
    cat(row, "\n")
  }
}

print_axis_table("p_star",        p_star_grid,        "p_star (causal SNPs)", width = 8)
print_axis_table("annotation_r2", annotation_r2_grid, "annotation_r2",        width = 8)
print_axis_table("noise",         noise_grid,         "noise_fraction",       width = 8)

cat("\n  --- By x_source ---\n")
hdr <- sprintf("%-22s", "method")
for (v in names(x_list)) hdr <- paste0(hdr, sprintf("  %-16s", v))
cat(hdr, "\n")
cat(strrep("-", nchar(hdr)), "\n")
for (m in methods) {
  row <- sprintf("%-22s", m)
  for (v in names(x_list)) {
    sub <- ok[ok$method == m & ok$x_source == v, ]
    row <- paste0(row, sprintf("  %16.3f", mean(sub$auprc, na.rm = TRUE)))
  }
  cat(row, "\n")
}

cat("\n")
cat("=======================================================================\n")
cat(sprintf("SECTION C -- Per-effect distribution: c, sigma_0_2, |mu_eff|\n"))
cat(sprintf("  (%d observations per method: %d runs x %d effects)\n",
            total_runs * L, total_runs, L))
cat("=======================================================================\n")
cat(sprintf("%-22s  %8s  %8s  %8s  %8s  %9s  %9s\n",
            "method", "c_mean", "c_sd", "sig_mean", "sig_sd", "|mu|_mean", "|mu|_sd"))
cat(strrep("-", 82), "\n")
for (m in methods) {
  sub <- effect_results[effect_results$method == m, ]
  cat(sprintf("%-22s  %8.3f  %8.3f  %8.4f  %8.4f  %9.4f  %9.4f\n",
              m,
              mean(sub$c_l,          na.rm = TRUE),
              sd(sub$c_l,            na.rm = TRUE),
              mean(sub$sigma_0_2_l,  na.rm = TRUE),
              sd(sub$sigma_0_2_l,    na.rm = TRUE),
              mean(abs(sub$mu_eff_l), na.rm = TRUE),
              sd(abs(sub$mu_eff_l),  na.rm = TRUE)))
}

cat("\n  --- Per-effect sigma_0_2 mean by p_star ---\n")
hdr <- sprintf("%-22s", "method")
for (v in p_star_grid) hdr <- paste0(hdr, sprintf("  sig@p*=%2d", v))
cat(hdr, "\n", strrep("-", nchar(hdr)), "\n")
for (m in methods) {
  row <- sprintf("%-22s", m)
  for (v in p_star_grid) {
    sub <- effect_results[effect_results$method == m & effect_results$p_star == v, ]
    row <- paste0(row, sprintf("  %9.4f", mean(sub$sigma_0_2_l, na.rm = TRUE)))
  }
  cat(row, "\n")
}

cat("\n  --- Per-effect mean |c| by annotation_r2 ---\n")
hdr <- sprintf("%-22s", "method")
for (v in annotation_r2_grid) hdr <- paste0(hdr, sprintf("  |c|@r2=%.1f", v))
cat(hdr, "\n", strrep("-", nchar(hdr)), "\n")
for (m in methods) {
  row <- sprintf("%-22s", m)
  for (v in annotation_r2_grid) {
    sub <- effect_results[effect_results$method == m & effect_results$annotation_r2 == v, ]
    row <- paste0(row, sprintf("  %11.3f", mean(abs(sub$c_l), na.rm = TRUE)))
  }
  cat(row, "\n")
}

cat("\n")
cat("=======================================================================\n")
cat("SECTION D -- Errors\n")
cat("=======================================================================\n")
if (nrow(err) == 0) {
  cat("No errors.\n")
} else {
  cat(sprintf("%d error(s):\n", nrow(err)))
  for (i in seq_len(nrow(err))) {
    cat(sprintf("  run=%d  x=%s  p*=%d  r2=%.1f  noise=%.1f  method=%s\n    %s\n",
                err$run_id[i], err$x_source[i], err$p_star[i],
                err$annotation_r2[i], err$noise[i], err$method[i],
                err$error[i]))
  }
}

cat("\nDone.\n")
