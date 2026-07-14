#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(data.table); library(digest); library(jsonlite) })

parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Use --name value pairs")
    out[[substring(x[[i]], 3L)]] <- x[[i + 1L]]; i <- i + 2L
  }
  out
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("result", "signatures", "generator", "out", "program")
if (length(setdiff(required, names(args))) || length(setdiff(names(args), required))) stop("Argument mismatch")

result_dir <- normalizePath(args$result); out_dir <- normalizePath(args$out, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out_dir, "validation_report.tsv"))) stop("Refusing to overwrite validation output")
sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
checks <- list()
add <- function(name, pass, observed, expected) checks[[length(checks) + 1L]] <<-
  data.table(check = name, passed = isTRUE(pass), observed = as.character(observed), expected = as.character(expected))

required_files <- c("FROZEN_SHA256.tsv", "tier_manifest.tsv", "diagnostic_summary.json",
                    "label_blind_audit.json", "input_access_manifest.tsv", "seed_manifest.tsv",
                    "chain_diagnostics.tsv", "target_position_strata.tsv", "covariance_diagnostic.tsv",
                    "selected_null_sets_summary.tsv", "selected_null_gene_sets.tsv",
                    "retained_unique_sets_summary.tsv", "retained_unique_gene_sets.tsv",
                    "jaccard_summary.tsv", "sessionInfo.txt")
missing <- required_files[!file.exists(file.path(result_dir, required_files))]
add("required_files_present", !length(missing), paste(missing, collapse = ";"), "none")
if (length(missing)) stop("Required result files missing")

frozen <- fread(file.path(result_dir, "FROZEN_SHA256.tsv"))
actual <- vapply(file.path(result_dir, frozen$file), sha256_file, character(1))
add("frozen_result_hashes", all(actual == frozen$sha256), sum(actual == frozen$sha256), nrow(frozen))

tier <- fread(file.path(result_dir, "tier_manifest.tsv"))
add("tier_single_row", nrow(tier) == 1L, nrow(tier), 1L)
add("program_identity", tier$program[[1]] == args$program, tier$program[[1]], args$program)
add("tier_tolerance", tier$tolerance_pc1[[1]] == 0.05 && tier$tolerance_mean_r[[1]] == 0.05,
    paste(tier$tolerance_pc1[[1]], tier$tolerance_mean_r[[1]], sep = "/"), "0.05/0.05")
add("twenty_chains", tier$chains_attempted[[1]] == 20L, tier$chains_attempted[[1]], 20L)

estimable <- tier$status[[1]] == "estimable_at_0.05"
valid_status <- estimable || tier$status[[1]] == "not_estimable_at_0.05"
add("valid_tier_status", valid_status, tier$status[[1]], "estimable_at_0.05 or not_estimable_at_0.05")
if (estimable) {
  add("estimable_gate", tier$unique_sets[[1]] >= 1000L && tier$selected_sets[[1]] == 1000L &&
        tier$chains_contributing[[1]] >= 5L && isTRUE(tier$selected_for_scoring[[1]]),
      paste(tier$unique_sets[[1]], tier$selected_sets[[1]], tier$chains_contributing[[1]], sep = "/"),
      ">=1000/1000/>=5")
} else {
  legitimate_failure <- tier$chains_feasible[[1]] < 5L || tier$chains_contributing[[1]] < 5L || tier$unique_sets[[1]] < 1000L
  add("nonestimable_gate", legitimate_failure && !isTRUE(tier$selected_for_scoring[[1]]) && tier$selected_sets[[1]] == 0L,
      paste(tier$chains_feasible[[1]], tier$chains_contributing[[1]], tier$unique_sets[[1]], sep = "/"),
      "prespecified failure condition and zero selected")
}

target <- fread(file.path(result_dir, "target_position_strata.tsv"))
selected_summary <- fread(file.path(result_dir, "selected_null_sets_summary.tsv"))
selected_long <- fread(file.path(result_dir, "selected_null_gene_sets.tsv"))
if (estimable) {
  add("selected_set_count", nrow(selected_summary) == 1000L && uniqueN(selected_summary$gene_key) == 1000L,
      paste(nrow(selected_summary), uniqueN(selected_summary$gene_key), sep = "/"), "1000/1000")
  sizes <- selected_long[, .N, by = set_id]
  add("exact_gene_count", nrow(sizes) == 1000L && all(sizes$N == nrow(target)),
      paste(range(sizes$N), collapse = "-"), nrow(target))
  expected_pos <- target[, .(position, expected_target = target_gene, expected_direction = direction,
                             expected_mean = mean_decile, expected_sd = sd_decile, expected_stratum = joint_stratum)]
  mapped <- merge(selected_long, expected_pos, by = "position", all.x = TRUE)
  position_ok <- all(mapped$target_position_gene == mapped$expected_target &
                       mapped$direction == mapped$expected_direction &
                       mapped$mean_decile == mapped$expected_mean & mapped$sd_decile == mapped$expected_sd &
                       mapped$joint_stratum == mapped$expected_stratum)
  add("position_constraints", position_ok, position_ok, TRUE)
  max_pc <- max(selected_summary$abs_delta_pc1); max_r <- max(selected_summary$abs_delta_mean_r)
  add("pc1_tolerance", max_pc <= 0.05 + 1e-12, format(max_pc, digits = 15), "<=0.05")
  add("mean_r_tolerance", max_r <= 0.05 + 1e-12, format(max_r, digits = 15), "<=0.05")
} else {
  max_pc <- NA_real_; max_r <- NA_real_
  add("no_selected_sets_when_not_estimable", nrow(selected_summary) == 0L && nrow(selected_long) == 0L,
      paste(nrow(selected_summary), nrow(selected_long), sep = "/"), "0/0")
}

registry <- fread(args$signatures)
retained_long <- fread(file.path(result_dir, "retained_unique_gene_sets.tsv"))
tested_genes <- unique(c(selected_long$gene, retained_long$gene))
overlap <- intersect(tested_genes, unique(registry$gene))
add("nine_program_union_excluded", !length(overlap), paste(overlap, collapse = ";"), "none")
add("complete_registry", uniqueN(registry$signature_id) == 9L, uniqueN(registry$signature_id), 9L)

chain <- fread(file.path(result_dir, "chain_diagnostics.tsv"))
add("chain_diagnostics", nrow(chain) == 20L && uniqueN(chain$chain_id) == 20L,
    paste(nrow(chain), uniqueN(chain$chain_id), sep = "/"), "20/20")
if (estimable) add("selected_chain_balance", sum(chain$selected_unique_sets) == 1000L && sum(chain$selected_unique_sets > 0L) >= 5L,
                   paste(sum(chain$selected_unique_sets), sum(chain$selected_unique_sets > 0L), sep = "/"), "1000/>=5")

label <- read_json(file.path(result_dir, "label_blind_audit.json"), simplifyVector = TRUE)
add("label_blind", isTRUE(label$passed) && !isTRUE(label$clinical_label_files_loaded) && !isTRUE(label$endpoint_scoring_performed),
    paste(label$passed, label$clinical_label_files_loaded, label$endpoint_scoring_performed, sep = "/"), "TRUE/FALSE/FALSE")
inputs <- fread(file.path(result_dir, "input_access_manifest.tsv"))
expected_roles <- c("expression_features", "complete_frozen_program_registry")
add("declared_inputs_only", nrow(inputs) == 2L && setequal(inputs$input_role, expected_roles) && all(!inputs$contains_clinical_labels),
    paste(inputs$input_role, collapse = ";"), paste(expected_roles, collapse = ";"))
add("input_hashes", all(vapply(inputs$path, sha256_file, character(1)) == inputs$sha256), "computed", "all match")

gene_file <- if (estimable) file.path(result_dir, "selected_null_gene_sets.tsv") else file.path(result_dir, "retained_unique_gene_sets.tsv")
add("gene_manifest_hash", sha256_file(gene_file) == tier$gene_set_manifest_sha256[[1]], sha256_file(gene_file), tier$gene_set_manifest_sha256[[1]])
add("generator_hash", sha256_file(args$generator) == tier$script_sha256[[1]], sha256_file(args$generator), tier$script_sha256[[1]])

jaccard <- fread(file.path(result_dir, "jaccard_summary.tsv"))
expected_sets <- if (estimable) 1000L else tier$unique_sets[[1]]
add("jaccard_set_count", jaccard$n_sets[[1]] == expected_sets, jaccard$n_sets[[1]], expected_sets)

report <- rbindlist(checks); fwrite(report, file.path(out_dir, "validation_report.tsv"), sep = "\t")
passed <- all(report$passed)
summary <- list(passed = passed, program = args$program, tier_status = tier$status[[1]],
                checks_passed = sum(report$passed), checks_total = nrow(report),
                failed_checks = report[passed == FALSE, check], result_dir = result_dir,
                maximum_absolute_pc1_deviation = max_pc, maximum_absolute_mean_r_deviation = max_r)
write_json(summary, file.path(out_dir, "validation_summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = 15)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
paths <- c(file.path(out_dir, "validation_report.tsv"), file.path(out_dir, "validation_summary.json"), script_path)
paths <- paths[!is.na(paths) & file.exists(paths)]
fwrite(data.table(file = normalizePath(paths), sha256 = vapply(paths, sha256_file, character(1))),
       file.path(out_dir, "VALIDATION_SHA256.tsv"), sep = "\t")
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE, digits = 15), "\n")
if (!passed) quit(status = 2L)
