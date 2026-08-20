test_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
test_dir <- dirname(normalizePath(test_file))
source(file.path(test_dir, "../src/burden_index.R"))
source(file.path(test_dir, "sim_cna.R"))

pearson_spearman_bias <- function(estimate, truth) {
  c(
    pearson = cor(estimate, truth, method = "pearson"),
    spearman = cor(estimate, truth, method = "spearman"),
    median_abs_relative_bias = median(abs(estimate - truth) / pmax(truth, 1e-8))
  )
}

cat("V1 amplitude recovery\n")
design <- expand.grid(
  f = c(0.02, 0.05, 0.10, 0.20, 0.40),
  events = c(4L, 12L, 30L),
  sigma = c(0.08, 0.15)
)
results <- lapply(seq_len(nrow(design)), function(index) {
  simulation <- simulate_cna(
    design$f[index], design$events[index], design$sigma[index],
    seed = 100 + index
  )
  amplitude <- burden_amplitude(simulation$chr, simulation$logR)
  single <- sqrt(max(
    0, var(winsorise(simulation$logR)) - amplitude$sigma_noise_dlrs^2
  ))
  c(
    truth = truth_amplitude(simulation),
    variance_components = amplitude$amp,
    single_subtraction = single,
    uncorrected_sd = sd(simulation$logR)
  )
})
results <- as.data.frame(do.call(rbind, results))
for (name in c("variance_components", "single_subtraction", "uncorrected_sd")) {
  metrics <- pearson_spearman_bias(results[[name]], results$truth)
  cat(sprintf(
    "%s: Pearson %.3f; Spearman %.3f; median |relative bias| %.3f\n",
    name, metrics[1L], metrics[2L], metrics[3L]
  ))
}

cat("\nV2 noise invariance\n")
noise_sweep <- c(0.05, 0.08, 0.12, 0.15, 0.20, 0.25)
noise_results <- vapply(noise_sweep, function(sigma) {
  simulation <- simulate_cna(0.2, 12, sigma, seed = 210)
  amplitude <- burden_amplitude(simulation$chr, simulation$logR)
  c(burden_amp = amplitude$amp, sigma_recovered = amplitude$sigma_noise)
}, numeric(2))
for (index in seq_along(noise_sweep)) {
  cat(sprintf(
    "sigma %.2f -> recovered %.3f; burden_amp %.3f\n",
    noise_sweep[index], noise_results[2L, index], noise_results[1L, index]
  ))
}

cat("\nV3 offset invariance\n")
offset_simulation <- simulate_cna(0.2, 12, 0.12, seed = 310)
offsets <- c(0, log2(3 / 2), -1, 1)
offset_amps <- vapply(offsets, function(offset) {
  burden_amplitude(
    offset_simulation$chr, offset_simulation$logR + offset
  )$amp
}, numeric(1))
cat(sprintf("maximum absolute difference: %.3g\n",
            max(abs(offset_amps - offset_amps[1L]))))

cat("\nV4 type I error and p-value calibration\n")
type_i_replicates <- 300L
p_values <- vapply(seq_len(type_i_replicates), function(index) {
  simulation <- simulate_flat(0.12, seed = 4000 + index)
  permutation_test(
    simulation$chr, simulation$logR, n_perm = 99, seed = 5000 + index
  )$p
}, numeric(1))
cat(sprintf(
  "FPR 0.05: %.3f; FPR 0.10: %.3f; KS p: %.3f\n",
  mean(p_values <= 0.05), mean(p_values <= 0.10),
  ks.test(p_values, "punif")$p.value
))

cat("\nV5 power and LOD\n")
power_design <- expand.grid(
  f = c(0.02, 0.05, 0.10, 0.20),
  events = c(4L, 12L, 30L)
)
for (row in seq_len(nrow(power_design))) {
  tests <- lapply(seq_len(20L), function(replicate) {
    simulation <- simulate_cna(
      power_design$f[row], power_design$events[row], 0.12,
      seed = 6000 + row * 100 + replicate
    )
    permutation_test(
      simulation$chr, simulation$logR, n_perm = 99,
      seed = 7000 + row * 100 + replicate
    )
  })
  p <- vapply(tests, `[[`, numeric(1), "p")
  lod <- vapply(tests, `[[`, numeric(1), "lod95_amp")
  cat(sprintf(
    "f %.2f events %d: power %.2f; median LOD %.3f\n",
    power_design$f[row], power_design$events[row], mean(p < 0.05),
    median(lod)
  ))
}

cat("\nV6 focal versus broad landscapes\n")
for (landscape in c("broad", "focal")) {
  simulation <- simulate_cna(
    0.2, 12, 0.12, seed = 8100, focal = landscape == "focal"
  )
  estimate <- burden_amplitude(simulation$chr, simulation$logR)$amp
  cat(sprintf(
    "%s: truth %.3f; estimate %.3f\n",
    landscape, truth_amplitude(simulation), estimate
  ))
}

cat("\nV7 block-bootstrap coverage\n")
covered <- vapply(seq_len(100L), function(index) {
  simulation <- simulate_cna(0.2, 12, 0.12, seed = 9000 + index)
  truth <- truth_amplitude(simulation)
  interval <- bootstrap_amplitude(
    simulation$chr, simulation$logR, n_boot = 100, block = 20,
    seed = 10000 + index
  )
  interval$lo95 <= truth && interval$hi95 >= truth
}, logical(1))
cat(sprintf("coverage: %.3f\n", mean(covered)))
