#!/usr/bin/env Rscript

# General label-blind covariance-matched program-identity null generator.
# Frozen for the nine-program GSE151179 audit at tolerance 0.05.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
  library(parallel)
  library(Matrix)
})

parse_args <- function(x) {
  out <- list(); i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) stop("Use --name value pairs")
    out[[substring(x[[i]], 3L)]] <- x[[i + 1L]]; i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("expr", "signatures", "out", "program", "tolerance", "workers")
if (length(setdiff(required, names(args))) || length(setdiff(names(args), required))) stop("Argument mismatch")
if (basename(args$expr) != "GSE151179_primary_preRAI_gene_expression.tsv.gz") stop("Unregistered expression matrix")
if (basename(args$signatures) != "frozen_programs_all9.tsv") stop("The complete nine-program registry is required")
if (any(grepl("sample|metadata|phenotype|clinical|endpoint|group|response|uptake|label",
              basename(c(args$expr, args$signatures)), ignore.case = TRUE))) stop("Label-blind filename gate failed")
if (!file.exists(args$expr) || !file.exists(args$signatures)) stop("Input missing")

tolerance <- as.numeric(args$tolerance); workers <- as.integer(args$workers)
if (!isTRUE(all.equal(tolerance, 0.05, tolerance = 1e-12))) stop("Only the frozen 0.05 tier is permitted")
if (!is.finite(workers) || workers < 1L || workers > 10L) stop("workers must be 1-10")

allowed_programs <- c(
  "TDS_16", "IODIDE_HANDLING_11", "CONDELLO_2025_SIX",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS", "HALLMARK_INFLAMMATORY_RESPONSE"
)
if (!args$program %in% allowed_programs) stop("Program is outside the frozen family")

out_dir <- normalizePath(args$out, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (any(file.exists(file.path(out_dir, c("tier_manifest.tsv", "diagnostic_summary.json"))))) {
  stop("Refusing to overwrite a prior tier run")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

# Amendment 01 constants.
n_chains <- 20L; max_anneal_proposals <- 20000L
burnin_accepted <- 5000L; thin_accepted <- 50L
max_within_proposals <- 100000L; round_proposals <- 10000L
target_unique <- 1000L; minimum_contributing_chains <- 5L
seed_map <- c(
  TDS_16 = 151179051L,
  IODIDE_HANDLING_11 = 151179052L,
  CONDELLO_2025_SIX = 151179050L,
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = 151179061L,
  HALLMARK_HYPOXIA = 151179062L,
  HALLMARK_ANGIOGENESIS = 151179063L,
  HALLMARK_G2M_CHECKPOINT = 151179064L,
  HALLMARK_E2F_TARGETS = 151179065L,
  HALLMARK_INFLAMMATORY_RESPONSE = 151179066L
)
master_seed <- unname(seed_map[[args$program]])
chain_seeds <- master_seed + seq_len(n_chains) * 100003L

expr_dt <- fread(args$expr, check.names = FALSE)
if (ncol(expr_dt) != 18L || tolower(names(expr_dt)[1]) != "symbol") stop("Expression schema mismatch")
symbols <- as.character(expr_dt[[1]])
if (anyNA(symbols) || anyDuplicated(symbols)) stop("Expression symbols must be unique")
sample_ids <- names(expr_dt)[-1L]
if (!all(grepl("^GSM[0-9]+$", sample_ids))) stop("Expression columns are not GEO accessions only")
expr <- as.matrix(expr_dt[, -1L]); storage.mode(expr) <- "double"; rownames(expr) <- symbols
if (any(!is.finite(expr))) stop("Expression contains non-finite values")
rm(expr_dt)

registry <- fread(args$signatures)
required_cols <- c("signature_id", "gene", "direction", "lock_status")
if (length(setdiff(required_cols, names(registry)))) stop("Registry schema mismatch")
if (!setequal(unique(registry$signature_id), allowed_programs) || uniqueN(registry$signature_id) != 9L) {
  stop("Registry must contain exactly the frozen nine-program family")
}
registry[, direction := as.integer(direction)]
if (any(registry$lock_status != "locked") || any(!registry$direction %in% c(-1L, 1L))) stop("Registry lock failure")

gene_mean <- rowMeans(expr); gene_sd <- apply(expr, 1L, sd)
eligible_numeric <- is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0
target_requested <- registry[signature_id == args$program]
target <- target_requested[gene %in% rownames(expr) & eligible_numeric[match(gene, rownames(expr))]]
if (nrow(target) < max(3L, ceiling(0.70 * nrow(target_requested)))) stop("Insufficient mapped target genes")
target[, direction := as.integer(direction)]

real_gene_union <- unique(registry$gene)
universe_keep <- eligible_numeric & !rownames(expr) %in% real_gene_union
universe_genes <- rownames(expr)[universe_keep]
if (length(universe_genes) < 1000L) stop("Eligible universe unexpectedly small")

ecdf_decile <- function(x, reference) {
  p <- vapply(x, function(v) mean(reference <= v), numeric(1))
  pmax(1L, pmin(10L, as.integer(ceiling(10 * p))))
}
u_mean <- gene_mean[universe_genes]; u_sd <- gene_sd[universe_genes]
u_mean_decile <- ecdf_decile(u_mean, u_mean); u_sd_decile <- ecdf_decile(u_sd, u_sd)
universe_stratum <- sprintf("M%02d_S%02d", u_mean_decile, u_sd_decile)
t_mean <- gene_mean[target$gene]; t_sd <- gene_sd[target$gene]
t_mean_decile <- ecdf_decile(t_mean, u_mean); t_sd_decile <- ecdf_decile(t_sd, u_sd)
target_stratum <- sprintf("M%02d_S%02d", t_mean_decile, t_sd_decile)

universe_expr <- expr[universe_genes, , drop = FALSE]
universe_z <- (universe_expr - rowMeans(universe_expr)) / apply(universe_expr, 1L, sd)
target_expr <- expr[target$gene, , drop = FALSE]
target_z <- (target_expr - rowMeans(target_expr)) / apply(target_expr, 1L, sd)
directions <- target$direction; n_samples <- ncol(expr); p_genes <- nrow(target)

features_from_moments <- function(gram, sum_vector) {
  ev <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
  pc1 <- max(ev) / sum(ev)
  mean_r <- (sum(sum_vector^2) / (n_samples - 1) - p_genes) / (p_genes * (p_genes - 1))
  c(pc1 = pc1, mean_r = mean_r)
}

make_state_object <- function(indices) {
  oriented <- universe_z[indices, , drop = FALSE] * directions
  gram <- crossprod(oriented) / (n_samples - 1)
  sum_vector <- colSums(oriented)
  list(indices = indices, gram = gram, sum_vector = sum_vector,
       features = features_from_moments(gram, sum_vector))
}

target_oriented <- target_z * directions
target_gram <- crossprod(target_oriented) / (n_samples - 1)
target_sum <- colSums(target_oriented)
target_cov <- features_from_moments(target_gram, target_sum)
if (any(!is.finite(target_cov))) stop("Target covariance features are not finite")

position_pools <- lapply(target_stratum, function(s) which(universe_stratum == s))
needed <- table(target_stratum); available <- table(universe_stratum)
available_needed <- setNames(rep(0L, length(needed)), names(needed))
available_needed[names(available)[names(available) %in% names(needed)]] <- available[names(available) %in% names(needed)]
if (any(needed > available_needed)) stop("At least one joint stratum has an insufficient candidate pool")

make_initial_indices <- function() {
  selected <- integer(length(position_pools))
  for (s in unique(target_stratum)) {
    pos <- which(target_stratum == s); pool <- which(universe_stratum == s)
    selected[pos] <- sample(pool, length(pos), replace = FALSE)
  }
  selected
}

is_feasible <- function(features) {
  is.finite(features[["pc1"]]) && is.finite(features[["mean_r"]]) &&
    abs(features[["pc1"]] - target_cov[["pc1"]]) <= tolerance &&
    abs(features[["mean_r"]] - target_cov[["mean_r"]]) <= tolerance
}
objective <- function(features) {
  ((features[["pc1"]] - target_cov[["pc1"]]) / tolerance)^2 +
    ((features[["mean_r"]] - target_cov[["mean_r"]]) / tolerance)^2
}

propose_object <- function(obj) {
  position <- sample.int(length(obj$indices), 1L)
  choices <- setdiff(position_pools[[position]], obj$indices)
  if (!length(choices)) return(NULL)
  new_index <- sample(choices, 1L); old_index <- obj$indices[[position]]
  old_vector <- universe_z[old_index, ] * directions[[position]]
  new_vector <- universe_z[new_index, ] * directions[[position]]
  gram <- obj$gram + (tcrossprod(new_vector) - tcrossprod(old_vector)) / (n_samples - 1)
  sum_vector <- obj$sum_vector + new_vector - old_vector
  indices <- obj$indices; indices[[position]] <- new_index
  list(indices = indices, gram = gram, sum_vector = sum_vector,
       features = features_from_moments(gram, sum_vector))
}

anneal_chain <- function(chain_id) {
  set.seed(chain_seeds[[chain_id]])
  obj <- make_state_object(make_initial_indices())
  current_objective <- objective(obj$features); initial_objective <- current_objective
  accepted <- 0L; first_feasible <- if (is_feasible(obj$features)) 0L else NA_integer_
  if (is.na(first_feasible)) {
    initial_temperature <- max(1, current_objective)
    for (step in seq_len(max_anneal_proposals)) {
      proposed <- propose_object(obj); if (is.null(proposed)) next
      proposed_objective <- objective(proposed$features)
      temperature <- initial_temperature * (0.001 ^ (step / max_anneal_proposals))
      if (proposed_objective <= current_objective || runif(1) < exp((current_objective - proposed_objective) / temperature)) {
        obj <- proposed; current_objective <- proposed_objective; accepted <- accepted + 1L
      }
      if (is_feasible(obj$features)) { first_feasible <- step; break }
    }
  }
  list(chain_id = chain_id, seed = chain_seeds[[chain_id]], feasible = !is.na(first_feasible),
       first_feasible_proposal = first_feasible,
       anneal_proposals = if (is.na(first_feasible)) max_anneal_proposals else first_feasible,
       anneal_accepted = accepted, initial_objective = initial_objective, obj = obj,
       within_proposals = 0L, feasible_accepted = 0L, accepted_since_last_retain = 0L,
       retained_total = 0L, rng_state = .Random.seed)
}

run_sampling_round <- function(chain_state, n_proposals) {
  assign(".Random.seed", chain_state$rng_state, envir = .GlobalEnv)
  records <- list()
  if (!chain_state$feasible) return(list(state = chain_state, records = records))
  for (step in seq_len(n_proposals)) {
    proposed <- propose_object(chain_state$obj)
    chain_state$within_proposals <- chain_state$within_proposals + 1L
    if (is.null(proposed) || !is_feasible(proposed$features)) next
    chain_state$obj <- proposed; chain_state$feasible_accepted <- chain_state$feasible_accepted + 1L
    if (chain_state$feasible_accepted <= burnin_accepted) next
    chain_state$accepted_since_last_retain <- chain_state$accepted_since_last_retain + 1L
    if (chain_state$accepted_since_last_retain >= thin_accepted) {
      chain_state$accepted_since_last_retain <- 0L
      chain_state$retained_total <- chain_state$retained_total + 1L
      records[[length(records) + 1L]] <- list(
        chain_id = chain_state$chain_id, chain_retain_index = chain_state$retained_total,
        indices = chain_state$obj$indices,
        pc1 = unname(chain_state$obj$features[["pc1"]]),
        mean_r = unname(chain_state$obj$features[["mean_r"]])
      )
    }
  }
  chain_state$rng_state <- .Random.seed
  list(state = chain_state, records = records)
}

chains <- mclapply(seq_len(n_chains), anneal_chain, mc.cores = workers, mc.preschedule = FALSE)
if (any(vapply(chains, inherits, logical(1), "try-error"))) stop("Annealing chain failure")
feasible_chain_count <- sum(vapply(chains, function(x) x$feasible, logical(1)))
all_records <- list(); unique_records <- list()
if (feasible_chain_count >= minimum_contributing_chains) {
  for (round_id in seq_len(ceiling(max_within_proposals / round_proposals))) {
    remaining <- max_within_proposals - vapply(chains, function(x) x$within_proposals, integer(1))
    n_this <- pmin(round_proposals, pmax(0L, remaining))
    active <- which(vapply(chains, function(x) x$feasible, logical(1)) & n_this > 0L)
    if (!length(active)) break
    rr <- mclapply(active, function(i) run_sampling_round(chains[[i]], n_this[[i]]),
                   mc.cores = min(workers, length(active)), mc.preschedule = FALSE)
    if (any(vapply(rr, inherits, logical(1), "try-error"))) stop("Sampling chain failure")
    for (j in seq_along(active)) { chains[[active[[j]]]] <- rr[[j]]$state; all_records <- c(all_records, rr[[j]]$records) }
    if (length(all_records)) {
      keys <- vapply(all_records, function(rec) paste(sort(universe_genes[rec$indices]), collapse = ";"), character(1))
      unique_records <- all_records[!duplicated(keys)]
      contributors <- uniqueN(vapply(unique_records, function(rec) rec$chain_id, integer(1)))
      if (length(unique_records) >= target_unique && contributors >= minimum_contributing_chains) break
    }
  }
}
if (length(all_records)) {
  keys <- vapply(all_records, function(rec) paste(sort(universe_genes[rec$indices]), collapse = ";"), character(1))
  unique_records <- all_records[!duplicated(keys)]
}
unique_contributors <- if (length(unique_records)) uniqueN(vapply(unique_records, function(rec) rec$chain_id, integer(1))) else 0L
estimable <- length(unique_records) >= target_unique && unique_contributors >= minimum_contributing_chains
status <- if (estimable) "estimable_at_0.05" else "not_estimable_at_0.05"
reason_code <- if (estimable) "success_1000_unique_5plus_chains" else if (feasible_chain_count < 5L) {
  "fewer_than_5_feasible_chains"
} else if (unique_contributors < 5L) "fewer_than_5_contributing_chains" else "fewer_than_1000_unique_sets"

selected_records <- list()
if (estimable) {
  ord <- order(vapply(unique_records, function(x) x$chain_retain_index, integer(1)),
               vapply(unique_records, function(x) x$chain_id, integer(1)))
  selected_records <- unique_records[ord[seq_len(target_unique)]]
}

records_summary <- function(records, prefix) {
  if (!length(records)) return(data.table())
  rbindlist(lapply(seq_along(records), function(i) {
    x <- records[[i]]
    data.table(set_id = sprintf("%s_%04d", prefix, i), chain_id = x$chain_id,
               chain_retain_index = x$chain_retain_index, pc1 = x$pc1, mean_r = x$mean_r,
               abs_delta_pc1 = abs(x$pc1 - target_cov[["pc1"]]),
               abs_delta_mean_r = abs(x$mean_r - target_cov[["mean_r"]]),
               gene_key = paste(sort(universe_genes[x$indices]), collapse = ";"))
  }))
}
records_long <- function(records, prefix) {
  if (!length(records)) return(data.table())
  rbindlist(lapply(seq_along(records), function(i) {
    x <- records[[i]]
    data.table(set_id = sprintf("%s_%04d", prefix, i), chain_id = x$chain_id,
               position = seq_along(x$indices), target_position_gene = target$gene,
               gene = universe_genes[x$indices], direction = directions,
               mean_decile = t_mean_decile, sd_decile = t_sd_decile, joint_stratum = target_stratum)
  }))
}
unique_summary <- records_summary(unique_records, "candidate"); unique_long <- records_long(unique_records, "candidate")
selected_summary <- records_summary(selected_records, "null"); selected_long <- records_long(selected_records, "null")
fwrite(unique_summary, file.path(out_dir, "retained_unique_sets_summary.tsv"), sep = "\t")
fwrite(unique_long, file.path(out_dir, "retained_unique_gene_sets.tsv"), sep = "\t")
fwrite(selected_summary, file.path(out_dir, "selected_null_sets_summary.tsv"), sep = "\t")
fwrite(selected_long, file.path(out_dir, "selected_null_gene_sets.tsv"), sep = "\t")

selected_counts <- if (nrow(selected_summary)) selected_summary[, .N, by = chain_id] else data.table(chain_id = integer(), N = integer())
chain_diag <- rbindlist(lapply(chains, function(x) data.table(
  chain_id = x$chain_id, seed = x$seed, entered_feasible_region = x$feasible,
  first_feasible_proposal = x$first_feasible_proposal, anneal_proposals = x$anneal_proposals,
  anneal_accepted = x$anneal_accepted, initial_objective = x$initial_objective,
  within_region_proposals = x$within_proposals, accepted_feasible_swaps = x$feasible_accepted,
  within_region_acceptance_rate = if (x$within_proposals) x$feasible_accepted / x$within_proposals else NA_real_,
  retained_raw = x$retained_total
)))
chain_diag <- merge(chain_diag, selected_counts, by = "chain_id", all.x = TRUE)
setnames(chain_diag, "N", "selected_unique_sets"); chain_diag[is.na(selected_unique_sets), selected_unique_sets := 0L]
fwrite(chain_diag, file.path(out_dir, "chain_diagnostics.tsv"), sep = "\t")

target_positions <- data.table(
  position = seq_len(nrow(target)), target_gene = target$gene, direction = directions,
  expression_mean = unname(t_mean), expression_sd = unname(t_sd), mean_decile = t_mean_decile,
  sd_decile = t_sd_decile, joint_stratum = target_stratum, candidate_pool_size = lengths(position_pools)
)
fwrite(target_positions, file.path(out_dir, "target_position_strata.tsv"), sep = "\t")

jaccard_records <- if (length(selected_records)) selected_records else unique_records
jaccard_values <- numeric()
if (length(jaccard_records) >= 2L) {
  gene_sets <- lapply(jaccard_records, function(x) universe_genes[x$indices])
  gene_levels <- unique(unlist(gene_sets, use.names = FALSE))
  membership <- sparseMatrix(
    i = rep(seq_along(gene_sets), lengths(gene_sets)),
    j = match(unlist(gene_sets, use.names = FALSE), gene_levels),
    x = 1,
    dims = c(length(gene_sets), length(gene_levels))
  )
  intersection <- as.matrix(tcrossprod(membership))
  sizes <- lengths(gene_sets); union_size <- outer(sizes, sizes, "+") - intersection
  jaccard_values <- (intersection / union_size)[upper.tri(intersection)]
}
jaccard_summary <- data.table(
  n_sets = length(jaccard_records), n_pairs = length(jaccard_values),
  minimum = if (length(jaccard_values)) min(jaccard_values) else NA_real_,
  q25 = if (length(jaccard_values)) unname(quantile(jaccard_values, 0.25)) else NA_real_,
  median = if (length(jaccard_values)) median(jaccard_values) else NA_real_,
  q75 = if (length(jaccard_values)) unname(quantile(jaccard_values, 0.75)) else NA_real_,
  maximum = if (length(jaccard_values)) max(jaccard_values) else NA_real_
)
fwrite(jaccard_summary, file.path(out_dir, "jaccard_summary.tsv"), sep = "\t")

fwrite(data.table(chain_id = seq_len(n_chains), seed = chain_seeds), file.path(out_dir, "seed_manifest.tsv"), sep = "\t")
input_access <- data.table(
  input_role = c("expression_features", "complete_frozen_program_registry"),
  path = normalizePath(c(args$expr, args$signatures)),
  sha256 = vapply(c(args$expr, args$signatures), sha256_file, character(1)),
  contains_clinical_labels = FALSE
)
fwrite(input_access, file.path(out_dir, "input_access_manifest.tsv"), sep = "\t")

selected_file <- file.path(out_dir, "selected_null_gene_sets.tsv")
retained_file <- file.path(out_dir, "retained_unique_gene_sets.tsv")
gene_file <- if (estimable) selected_file else retained_file
script_hash <- if (!is.na(script_path) && file.exists(script_path)) sha256_file(script_path) else NA_character_
tier_manifest <- data.table(
  cohort = "GSE151179", program = args$program, tier = tolerance,
  tolerance_pc1 = tolerance, tolerance_mean_r = tolerance, status = status,
  genes_requested = nrow(target_requested), genes_present = nrow(target),
  positive_directions = sum(directions == 1L), negative_directions = sum(directions == -1L),
  universe_genes = length(universe_genes), chains_attempted = n_chains,
  chains_feasible = feasible_chain_count, chains_contributing = unique_contributors,
  annealing_proposals = sum(chain_diag$anneal_proposals),
  within_region_proposals = sum(chain_diag$within_region_proposals),
  accepted_swaps = sum(chain_diag$accepted_feasible_swaps), retained_raw = length(all_records),
  unique_sets = length(unique_records), selected_sets = length(selected_records),
  selected_for_scoring = estimable, reason_code = reason_code,
  seed_manifest = file.path(out_dir, "seed_manifest.tsv"),
  gene_set_manifest_sha256 = sha256_file(gene_file), script_sha256 = script_hash
)
fwrite(tier_manifest, file.path(out_dir, "tier_manifest.tsv"), sep = "\t")

cov_source <- if (nrow(selected_summary)) selected_summary else unique_summary
cov_diag <- data.table(
  feature = c("pc1_variance_proportion", "mean_pairwise_pearson_r"), target = unname(target_cov),
  null_median = c(if (nrow(cov_source)) median(cov_source$pc1) else NA_real_,
                  if (nrow(cov_source)) median(cov_source$mean_r) else NA_real_),
  null_q025 = c(if (nrow(cov_source)) unname(quantile(cov_source$pc1, 0.025)) else NA_real_,
                if (nrow(cov_source)) unname(quantile(cov_source$mean_r, 0.025)) else NA_real_),
  null_q975 = c(if (nrow(cov_source)) unname(quantile(cov_source$pc1, 0.975)) else NA_real_,
                if (nrow(cov_source)) unname(quantile(cov_source$mean_r, 0.975)) else NA_real_),
  maximum_absolute_deviation = c(if (nrow(cov_source)) max(cov_source$abs_delta_pc1) else NA_real_,
                                 if (nrow(cov_source)) max(cov_source$abs_delta_mean_r) else NA_real_),
  tolerance = tolerance
)
fwrite(cov_diag, file.path(out_dir, "covariance_diagnostic.tsv"), sep = "\t")

label_audit <- list(
  passed = TRUE, clinical_label_files_loaded = FALSE,
  accepted_input_roles = c("expression_features", "complete_frozen_program_registry"),
  expression_columns_are_geo_accessions_only = TRUE, n_expression_samples = ncol(expr),
  null_generation_completed_before_endpoint_scoring = TRUE, endpoint_scoring_performed = FALSE
)
write_json(label_audit, file.path(out_dir, "label_blind_audit.json"), pretty = TRUE, auto_unbox = TRUE)
diagnostic <- list(
  cohort = "GSE151179", program = args$program, tolerance = tolerance,
  status = status, reason_code = reason_code,
  target = list(genes_requested = nrow(target_requested), genes_present = nrow(target),
                positive_directions = sum(directions == 1L), negative_directions = sum(directions == -1L),
                pc1_variance_proportion = unname(target_cov[["pc1"]]),
                mean_pairwise_pearson_r = unname(target_cov[["mean_r"]])),
  universe_gene_count = length(universe_genes),
  chains = list(attempted = n_chains, feasible = feasible_chain_count, contributing = unique_contributors),
  sets = list(retained_raw = length(all_records), retained_unique = length(unique_records),
              selected = length(selected_records),
              duplicate_rate = if (length(all_records)) 1 - length(unique_records) / length(all_records) else NA_real_),
  jaccard = as.list(jaccard_summary[1]),
  computational_budget = list(max_anneal_proposals_per_chain = max_anneal_proposals,
                              burnin_accepted_swaps = burnin_accepted,
                              thinning_accepted_swaps = thin_accepted,
                              max_within_region_proposals_per_chain = max_within_proposals,
                              workers = workers),
  label_blind = label_audit
)
write_json(diagnostic, file.path(out_dir, "diagnostic_summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = 15)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
files <- sort(list.files(out_dir, full.names = TRUE)); files <- files[basename(files) != "FROZEN_SHA256.tsv"]
fwrite(data.table(file = basename(files), sha256 = vapply(files, sha256_file, character(1))),
       file.path(out_dir, "FROZEN_SHA256.tsv"), sep = "\t")

cat(toJSON(list(status = status, reason_code = reason_code, program = args$program,
                 genes_present = nrow(target), target_pc1 = unname(target_cov[["pc1"]]),
                 target_mean_r = unname(target_cov[["mean_r"]]), chains_feasible = feasible_chain_count,
                 chains_contributing = unique_contributors, unique_sets = length(unique_records),
                 selected_sets = length(selected_records), out_dir = out_dir),
           pretty = TRUE, auto_unbox = TRUE, digits = 15), "\n")
