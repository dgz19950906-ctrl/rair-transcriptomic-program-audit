#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
})

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Use --name value pairs")
    out[[substring(x[[i]], 3L)]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("result", "signatures", "generator", "out")
if (length(setdiff(required, names(args)))) stop("Missing required arguments")
if (length(setdiff(names(args), required))) stop("Unsupported arguments")

result_dir <- normalizePath(args$result)
out_dir <- normalizePath(args$out, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(out_dir, "validation_report.tsv"))) stop("Refusing to overwrite validation output")

sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
checks <- list()
add_check <- function(name, passed, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name,
    passed = isTRUE(passed),
    observed = as.character(observed),
    expected = as.character(expected)
  )
}

required_files <- c(
  "FROZEN_SHA256.tsv", "tier_manifest.tsv", "diagnostic_summary.json",
  "label_blind_audit.json", "input_access_manifest.tsv", "seed_manifest.tsv",
  "chain_diagnostics.tsv", "target_position_strata.tsv", "covariance_diagnostic.tsv",
  "selected_null_sets_summary.tsv", "selected_null_gene_sets.tsv",
  "retained_unique_sets_summary.tsv", "retained_unique_gene_sets.tsv",
  "jaccard_summary.tsv", "sessionInfo.txt"
)
missing_files <- required_files[!file.exists(file.path(result_dir, required_files))]
add_check("required_files_present", !length(missing_files), paste(missing_files, collapse = ";"), "none missing")
if (length(missing_files)) stop("Missing required result files")

frozen <- fread(file.path(result_dir, "FROZEN_SHA256.tsv"))
hash_actual <- vapply(file.path(result_dir, frozen$file), sha256_file, character(1))
hash_ok <- identical(unname(hash_actual), frozen$sha256)
add_check("frozen_result_hashes", hash_ok, sum(hash_actual == frozen$sha256), nrow(frozen))

tier <- fread(file.path(result_dir, "tier_manifest.tsv"))
add_check("tier_manifest_single_row", nrow(tier) == 1L, nrow(tier), 1L)
add_check("tier_status", tier$status[[1]] == "estimable_at_0.05", tier$status[[1]], "estimable_at_0.05")
add_check("tier_reason", tier$reason_code[[1]] == "success_1000_unique_5plus_chains", tier$reason_code[[1]], "success_1000_unique_5plus_chains")
add_check("tier_tolerances", tier$tolerance_pc1[[1]] == 0.05 && tier$tolerance_mean_r[[1]] == 0.05,
          paste(tier$tolerance_pc1[[1]], tier$tolerance_mean_r[[1]], sep = "/"), "0.05/0.05")
add_check("tier_chain_gate", tier$chains_feasible[[1]] == 20L && tier$chains_contributing[[1]] >= 5L,
          paste(tier$chains_feasible[[1]], tier$chains_contributing[[1]], sep = "/"), "20 feasible and >=5 contributing")
add_check("tier_set_counts", tier$unique_sets[[1]] >= 1000L && tier$selected_sets[[1]] == 1000L,
          paste(tier$unique_sets[[1]], tier$selected_sets[[1]], sep = "/"), ">=1000 unique / 1000 selected")
add_check("selected_for_scoring", isTRUE(tier$selected_for_scoring[[1]]), tier$selected_for_scoring[[1]], TRUE)

selected_summary <- fread(file.path(result_dir, "selected_null_sets_summary.tsv"))
selected_long <- fread(file.path(result_dir, "selected_null_gene_sets.tsv"))
add_check("selected_summary_rows", nrow(selected_summary) == 1000L, nrow(selected_summary), 1000L)
add_check("selected_set_ids_unique", uniqueN(selected_summary$set_id) == 1000L, uniqueN(selected_summary$set_id), 1000L)
add_check("unordered_gene_keys_unique", uniqueN(selected_summary$gene_key) == 1000L, uniqueN(selected_summary$gene_key), 1000L)

set_sizes <- selected_long[, .N, by = set_id]
add_check("six_genes_per_set", nrow(set_sizes) == 1000L && all(set_sizes$N == 6L),
          paste(range(set_sizes$N), collapse = "-"), "6")
direction_counts <- selected_long[, .(positive = sum(direction == 1L), negative = sum(direction == -1L)), by = set_id]
add_check("direction_distribution", all(direction_counts$positive == 2L & direction_counts$negative == 4L),
          paste(unique(paste(direction_counts$positive, direction_counts$negative, sep = "/")), collapse = ";"), "2 positive / 4 negative")

signatures <- fread(args$signatures)
overlap <- intersect(unique(selected_long$gene), unique(signatures$gene))
add_check("real_program_union_excluded", !length(overlap), paste(overlap, collapse = ";"), "no overlap")

target_positions <- fread(file.path(result_dir, "target_position_strata.tsv"))
position_check <- merge(
  selected_long,
  target_positions[, .(position, expected_gene = target_gene, expected_direction = direction,
                        expected_mean_decile = mean_decile, expected_sd_decile = sd_decile,
                        expected_stratum = joint_stratum)],
  by = "position", all.x = TRUE
)
position_ok <- all(
  position_check$target_position_gene == position_check$expected_gene &
    position_check$direction == position_check$expected_direction &
    position_check$mean_decile == position_check$expected_mean_decile &
    position_check$sd_decile == position_check$expected_sd_decile &
    position_check$joint_stratum == position_check$expected_stratum
)
add_check("position_stratum_and_direction_preserved", position_ok, position_ok, TRUE)

max_pc <- max(selected_summary$abs_delta_pc1)
max_r <- max(selected_summary$abs_delta_mean_r)
add_check("pc1_tolerance", max_pc <= 0.05 + 1e-12, format(max_pc, digits = 15), "<=0.05")
add_check("mean_r_tolerance", max_r <= 0.05 + 1e-12, format(max_r, digits = 15), "<=0.05")

chain_diag <- fread(file.path(result_dir, "chain_diagnostics.tsv"))
add_check("twenty_chain_diagnostics", nrow(chain_diag) == 20L && uniqueN(chain_diag$chain_id) == 20L,
          paste(nrow(chain_diag), uniqueN(chain_diag$chain_id), sep = "/"), "20/20")
add_check("all_chains_feasible", all(chain_diag$entered_feasible_region), sum(chain_diag$entered_feasible_region), 20L)
add_check("selected_chain_contributions", sum(chain_diag$selected_unique_sets) == 1000L && sum(chain_diag$selected_unique_sets > 0L) >= 5L,
          paste(sum(chain_diag$selected_unique_sets), sum(chain_diag$selected_unique_sets > 0L), sep = "/"), "1000 sets / >=5 chains")

label_audit <- read_json(file.path(result_dir, "label_blind_audit.json"), simplifyVector = TRUE)
add_check("label_blind_gate", isTRUE(label_audit$passed) && !isTRUE(label_audit$clinical_label_files_loaded) &&
            !isTRUE(label_audit$endpoint_scoring_performed),
          paste(label_audit$passed, label_audit$clinical_label_files_loaded, label_audit$endpoint_scoring_performed, sep = "/"),
          "TRUE/FALSE/FALSE")
inputs <- fread(file.path(result_dir, "input_access_manifest.tsv"))
expected_roles <- c("expression_features", "frozen_program_definitions")
add_check("declared_inputs_only", nrow(inputs) == 2L && setequal(inputs$input_role, expected_roles) &&
            all(!inputs$contains_clinical_labels),
          paste(inputs$input_role, collapse = ";"), paste(expected_roles, collapse = ";"))
input_hash_ok <- all(vapply(inputs$path, sha256_file, character(1)) == inputs$sha256)
add_check("input_hashes_current", input_hash_ok, input_hash_ok, TRUE)

gene_hash <- sha256_file(file.path(result_dir, "selected_null_gene_sets.tsv"))
add_check("selected_gene_manifest_hash", gene_hash == tier$gene_set_manifest_sha256[[1]], gene_hash, tier$gene_set_manifest_sha256[[1]])
generator_hash <- sha256_file(args$generator)
add_check("generator_script_hash", generator_hash == tier$script_sha256[[1]], generator_hash, tier$script_sha256[[1]])

jaccard <- fread(file.path(result_dir, "jaccard_summary.tsv"))
add_check("jaccard_complete", jaccard$n_sets[[1]] == 1000L && jaccard$n_pairs[[1]] == 499500L,
          paste(jaccard$n_sets[[1]], jaccard$n_pairs[[1]], sep = "/"), "1000/499500")

report <- rbindlist(checks)
fwrite(report, file.path(out_dir, "validation_report.tsv"), sep = "\t")
passed <- all(report$passed)
summary <- list(
  passed = passed,
  checks_passed = sum(report$passed),
  checks_total = nrow(report),
  failed_checks = report[!passed, check],
  result_dir = result_dir,
  maximum_absolute_pc1_deviation = max_pc,
  maximum_absolute_mean_r_deviation = max_r
)
write_json(summary, file.path(out_dir, "validation_summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = 15)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
hash_paths <- c(file.path(out_dir, "validation_report.tsv"), file.path(out_dir, "validation_summary.json"), script_path)
hash_paths <- hash_paths[!is.na(hash_paths) & file.exists(hash_paths)]
hash_manifest <- data.table(file = normalizePath(hash_paths), sha256 = vapply(hash_paths, sha256_file, character(1)))
fwrite(hash_manifest, file.path(out_dir, "VALIDATION_SHA256.tsv"), sep = "\t")

cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE, digits = 15), "\n")
if (!passed) quit(status = 2L)
