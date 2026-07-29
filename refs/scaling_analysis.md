# SuSiNE Scaling Analysis

*Updated 2026-03-05. This file now includes completed prior-variance debugging from the provided vanilla output files, plus `mu_0` and functional-pi checks.*

## Status

- Prior-variance root cause: `CONFIRMED`
- Prior-variance fix: `CONFIRMED`
- `mu_0` scaling consistency: `CONFIRMED`
- Functional-pi with `c=0` behavior: `EXPECTED (not a bug)`

---

## 1. Executive Summary

1. The prior-variance bug was real and was enough to explain the vanilla `susine` vs `susie` mismatch.
2. The bug is exactly this: multiplying fixed `sigma_0_2` by `X_sd^2` in `initialize_priors` made the prior variance SNP-specific.
3. Removing that multiplication resolves the vanilla mismatch on the top-3 worst cases selected from:
   - `test_susine/output/run_table_vanilla.csv`
   - `test_susine/output/model_metrics_vanilla.csv`
4. `mu_0` scaling by `X_sd` is internally consistent and should stay.
5. `functional_pi` with `c=0` should still differ from vanilla, because `c` only affects `mu_0`, while functional pi changes inclusion weights `pi`.

---

## 2. Git Trace for the Bug

Relevant `susine` history for `R/initialize.R`:

- `efbca63`: introduced `sigma_0_2` scaling by `X_sd^2`.
- current working tree: removes that scaling and keeps only `mu_0 * X_sd`.

The old buggy line from `efbca63`:

```r
sigma_0_2 = sigma_0_2 * (X_sd^2)
```

Current intended behavior:

- `sigma_0_2` (fixed mode) is interpreted like `susieR::scaled_prior_variance`, i.e., a scalar multiplier of `var(y)`, not additionally SNP-scaled.
- `mu_0` is annotation-based and in original-X effect units, so `mu_0 * X_sd` conversion to standardized-X space is correct.

---

## 3. Top-3 Discrepant Pairs (From Provided Files)

Computed on `filtering == "unfiltered"` with matching `dataset_bundle_id` + `sigma_0_2_scalar`:

| dataset_bundle_id | sigma_0_2_scalar | run_id_susine | run_id_susie | AUPRC_susine | AUPRC_susie | abs diff |
|---:|---:|---:|---:|---:|---:|---:|
| 77  | 0.4 | 17787 | 17567 | 0.25 | 0.50 | 0.25 |
| 77  | 0.2 | 17786 | 17566 | 0.25 | 0.50 | 0.25 |
| 351 | 0.2 | 81080 | 80860 | 0.425641026 | 0.24496124 | 0.180679786 |

---

## 4. Reproduction and Validation

I added and ran:

- `test_susine/scripts/validate_prior_variance_scaling_from_outputs.R`

It does all of the following:

1. Finds top-3 discrepant pairs directly from the two provided vanilla CSVs.
2. Reconstructs the exact datasets (same `matrix_id`, `phenotype_seed`, `p_star`, `y_noise`, `sigma_0_2_scalar`).
3. Runs:
   - `susie_vanilla` baseline.
   - `susine` with fixed behavior (no `X_sd^2` prior-variance scaling).
   - `susine` with emulated old bug (injecting `sigma_0_2 * X_sd^2`).
4. Evaluates AUPRC using the same `evaluate_model` logic.

Output written to:

- `test_susine/output/prior_variance_top3_validation_results.csv`

### Key Result

For all 3 cases:

- `susine_fixed` matches `susie` (or is extremely close, numeric tolerance).
- `susine_bug` reproduces the recorded discrepant `susine` AUPRC values.

Example from output CSV:

- Bundle 77, sigma 0.2:
  - recorded: susie `0.5`, susine `0.25`
  - rerun fixed: susie `0.5`, susine_fixed `0.5`
  - rerun bug: susine_bug `0.25`
- Bundle 351, sigma 0.2:
  - recorded gap `0.180679786`
  - rerun fixed gap `0.0005167959`
  - rerun bug gap `0.1806797853`

Also confirms prior-variance distortion:

- `sigma0_ratio_fixed_max_over_min = 1`
- `sigma0_ratio_bug_max_over_min` about `27x` to `28x`

So the bug is causal, and the fix resolves it.

---

## 5. `mu_0` Scaling Consistency Check

Current logic:

- In `run_model.R`, functional mean is built as:
  - `mu_0 = c_value * annotation_vec` (original-X effect space).
- In `susine/R/initialize.R`, if `X` is provided:
  - `mu_0 <- mu_0 * X_sd` (convert to standardized-X effect space).

This is coherent with how `beta_hat`/`shat2` are computed in standardized-X space downstream. There is no analogous `mu_0` scaling bug found here.

Conclusion: do **not** remove `mu_0 * X_sd`; this part is correct.

---

## 6. Why `functional_pi` with `c=0` Differs from Vanilla

From `test_susine/R/run_model.R`:

- `c_value` only affects `mu_0` when `prior_mean_strategy == "functional_mu"`.
- `functional_pi` modifies inclusion weights via:
  - `pw <- softmax(abs(annotation_vec) / tau_value)`

So with `functional_pi`, setting `c=0` does not revert to vanilla unless inclusion priors are also set back to uniform. Different outputs vs vanilla are expected in that scenario.

---

## 7. Final Recommendation

1. Keep the prior-variance fix (no `X_sd^2` scaling for fixed scalar `sigma_0_2`).
2. Keep `mu_0 * X_sd` scaling as-is.
3. For comparisons intended to isolate `mu`, ensure `pi` is held uniform.
4. For comparisons intended to isolate `pi`, ignore `c` (it is not active in pure `functional_pi` use cases).

---

## 8. Follow-up Check: functional_mu with c=0 vs vanilla

Using:

- `test_susine/output/run_table_susine_vanilla_f_mu.csv`
- `test_susine/output/model_metrics_vanilla_f_mu.csv`

I ran:

- `test_susine/scripts/validate_fmu_c0_vs_vanilla_from_outputs.R`

Result:

- Matched pairs (`dataset_bundle_id`, `sigma_0_2_scalar`, `filtering`): `24,300`
- Maximum absolute difference across all checked metrics (`AUPRC`, `power`,
  `mean_size`, `mean_purity`, `mean_coverage`, `L_effective`,
  `cross_entropy`, `hg2`, `elbo_final`): `0`
- Offending c=0 pairs found: `none`

Diagnostic outputs:

- `test_susine/output/fmu_c0_vs_vanilla_pairs.csv`
- `test_susine/output/fmu_c0_vs_vanilla_offenders.csv`
- `test_susine/output/fmu_use_case_cvalue_summary.csv`

Interpretation:

- In these files, `functional_mu` with `c=0` is an exact replicate of vanilla.
- A likely source of apparent mismatch is collapsing `functional_mu` across all
  `c` values (including `0.3`, `0.6`, and rows with `c_value = NA`) rather than
  filtering to `c=0`.

---

## 9. Q6 Plot Root Cause and Provenance

### 9.1 Baseline mismatch in Q6

The original Q6 plotting code used `susie_vanilla` as the dashed baseline.
That creates an apparent gap at `c=0` whenever `susine_vanilla != susie_vanilla`
(which happened in this run because of the prior-variance bug).

With baseline sensitivity checks:

- `functional_mu(c=0) - susine_vanilla`: exactly `0` across all matched keys.
- `functional_mu(c=0) - susie_vanilla`: exactly equals
  `susine_vanilla - susie_vanilla` for every matched key.

So the visual discrepancy is a comparator artifact, not a `c=0` model failure.

### 9.2 Why `c_value = NA` appears

`c_value = NA` rows are introduced by run-control generation for
`exploration_methods = "single"`:

- In `test_susine/R/run_controls.R`, `axis_table_for_method("single", ...)`
  returns an empty tibble (no `c_value` column).
- Later, `make_job_config()` ensures missing columns exist and sets
  `runs$c_value <- NA_real_` if absent.

So functional-mu single-mode rows get `c_value = NA` by construction.

### 9.3 Why this is dangerous

In `test_susine/R/run_model.R`, `c_value` was previously parsed and defaulted
to `1` when non-finite. Therefore, `c_value = NA` can silently run as `c=1`
unless filtered or asserted in analysis.

### 9.4 Implemented analysis/plot safeguards

1. Q6 now uses `susine_vanilla` baseline (requested).
2. Q6 hard-filters `!is.na(c_value)` before plotting.
3. Q6 prints raw `c_value` counts and writes a numeric audit table.
4. Validator now writes explicit `c=NA` provenance and baseline-sensitivity
   CSVs:
   - `test_susine/output/fmu_c_na_provenance.csv`
   - `test_susine/output/q6_fmu_auprc_by_c_with_baselines.csv`
   - `test_susine/output/fmu_baseline_sensitivity_susine_vs_susie.csv`
   - `test_susine/output/fmu_baseline_sensitivity_c0_vs_susie.csv`
   - `test_susine/output/fmu_baseline_sensitivity_identity_check.csv`
