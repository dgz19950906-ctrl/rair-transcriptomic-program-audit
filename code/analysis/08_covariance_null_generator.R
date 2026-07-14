#!/usr/bin/env Rscript

# Label-blind covariance-matched program-identity null generator.
# Implements the frozen protocol and Amendment 01 (13 July 2026).
# This script deliberately accepts no phenotype, endpoint, group, or sample-metadata input.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
  library(parallel)
})

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) {
      stop("Arguments must be supplied as --name value pairs")
    }
    out[[substring(x[[i]], 3L)]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("expr", "signatures", "out", "program", "tolerance", "workers")
missing_args <- setdiff(required, names(args))
if (length(missing_args)) stop("Missing arguments: ", paste(missing_args, collapse = ", "))
if (length(setdiff(names(args), required))) {
  stop("Unsupported arguments: ", paste(setdiff(names(args), required), collapse = ", "))
}

forbidden_input_pattern <- "sample|metadata|phenotype|clinical|endpoint|group|response|uptake|label"
input_paths <- c(args$expr, args$signatures)
if (any(grepl(forbidden_input_pattern, basename(input_paths), ignore.case = TRUE))) {
  stop("Label-blind gate failed: an input filename resembles a phenotype or endpoint file")
}
if (basename(args$expr) != "GSE151179_primary_preRAI_gene_expression.tsv.gz") {
  stop("Unexpected expression input basename; refusing an unregistered matrix")
}

tolerance <- as.numeric(args$tolerance)
workers <- as.integer(args$workers)
if (!isTRUE(all.equal(tolerance, 0.05, tolerance = 1e-12))) {
  stop("This execution gate is frozen specifically for tolerance 0.05")
}
if (!is.finite(workers) || workers < 1L || workers > 10L) stop("workers must be between 1 and 10")
if (args$program != "CONDELLO_2025_SIX") {
  stop("This first diagnostic gate is frozen specifically for CONDELLO_2025_SIX")
}
if (!file.exists(args$expr) || !file.exists(args$signatures)) stop("An input file is missing")

out_dir <- normalizePath(args$out, mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sentinel_files <- file.path(out_dir, c("tier_manifest.tsv", "diagnostic_summary.json"))
if (any(file.exists(sentinel_files))) stop("Output sentinel exists; refusing to overwrite a prior tier run")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_

# Frozen Amendment 01 constants.
n_chains <- 20L
max_anneal_proposals <- 20000L
burnin_accepted <- 5000L
thin_accepted <- 50L
max_within_proposals <- 100000L
round_proposals <- 10000L
target_unique <- 1000L
minimum_contributing_chains <- 5L
master_seed <- 151179050L
chain_seeds <- master_seed + seq_len(n_chains) * 100003L

sha256_file <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

expr_dt <- fread(args$expr, check.names = FALSE)
if (ncol(expr_dt) != 18L) stop("Expected one symbol column plus 17 expression columns")
if (tolower(names(expr_dt)[[1]]) != "symbol") stop("First expression column must be symbol")
symbols <- as.character(expr_dt[[1]])
if (anyNA(symbols) || anyDuplicated(symbols)) stop("Expression symbols must be complete and unique")
sample_ids <- names(expr_dt)[-1L]
if (!all(grepl("^GSM[0-9]+$", sample_ids))) {
  stop("Label-blind gate failed: expression columns are not GEO sample accessions only")
}
expr <- as.matrix(expr_dt[, -1L])
storage.mode(expr) <- "double"
rownames(expr) <- symbols
if (any(!is.finite(expr))) stop("Expression matrix contains non-finite values")
rm(expr_dt)

signatures <- fread(args$signatures)
signature_required <- c("signature_id", "gene", "direction", "lock_status")
if (length(setdiff(signature_required, names(signatures)))) stop("Signature table schema mismatch")
if (any(signatures$lock_status != "locked")) stop("All program rows must be locked")
signatures[, direction := as.integer(direction)]
if (any(!signatures$direction %in% c(-1L, 1L))) stop("Directions must be -1 or +1")
target <- signatures[signature_id == args$program]
if (nrow(target) != 6L) stop("Condello-6 must contain exactly six frozen genes")
if (sum(target$direction == 1L) != 2L || sum(target$direction == -1L) != 4L) {
  stop("Condello-6 direction count mismatch")
}
if (length(setdiff(target$gene, rownames(expr)))) {
  stop("Target genes absent from expression matrix: ", paste(setdiff(target$gene, rownames(expr)), collapse = ","))
}

gene_mean <- rowMeans(expr)
gene_sd <- apply(expr, 1L, sd)
eligible_numeric <- is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0
real_gene_union <- unique(signatures$gene)
universe_keep <- eligible_numeric & !rownames(expr) %in% real_gene_union
universe_genes <- rownames(expr)[universe_keep]
if (length(universe_genes) < 1000L) stop("Eligible gene universe is unexpectedly small")

ecdf_decile <- function(x, reference) {
  p <- vapply(x, function(v) mean(reference <= v), numeric(1))
  pmax(1L, pmin(10L, as.integer(ceiling(10 * p))))
}

u_mean <- gene_mean[universe_genes]
u_sd <- gene_sd[universe_genes]
u_mean_decile <- ecdf_decile(u_mean, u_mean)
u_sd_decile <- ecdf_decile(u_sd, u_sd)
universe_stratum <- sprintf("M%02d_S%02d", u_mean_decile, u_sd_decile)

t_mean <- gene_mean[target$gene]
t_sd <- gene_sd[target$gene]
t_mean_decile <- ecdf_decile(t_mean, u_mean)
t_sd_decile <- ecdf_decile(t_sd, u_sd)
target_stratum <- sprintf("M%02d_S%02d", t_mean_decile, t_sd_decile)

universe_expr <- expr[universe_genes, , drop = FALSE]
universe_z <- (universe_expr - rowMeans(universe_expr)) / apply(universe_expr, 1L, sd)
target_expr <- expr[target$gene, , drop = FALSE]
target_z <- (target_expr - rowMeans(target_expr)) / apply(target_expr, 1L, sd)
directions <- target$direction

cov_features_matrix <- function(z_rows, dirs) {
  oriented <- z_rows * dirs
  cm <- cor(t(oriented))
  if (any(!is.finite(cm))) return(c(pc1 = NA_real_, mean_r = NA_real_))
  ev <- eigen(cm, symmetric = TRUE, only.values = TRUE)$values
  c(pc1 = max(ev) / sum(ev), mean_r = mean(cm[upper.tri(cm)]))
}

target_cov <- cov_features_matrix(target_z, directions)
if (any(!is.finite(target_cov))) stop("Target covariance features are not finite")

position_pools <- lapply(target_stratum, function(s) which(universe_stratum == s))
pool_sizes <- lengths(position_pools)
needed_by_stratum <- table(target_stratum)
available_by_stratum <- table(universe_stratum)
insufficient <- names(needed_by_stratum)[needed_by_stratum > available_by_stratum[names(needed_by_stratum)]]
if (length(insufficient)) {
  stop("Insufficient candidate genes in target strata: ", paste(insufficient, collapse = ","))
}

make_initial_state <- function() {
  selected <- integer(length(position_pools))
  for (s in unique(target_stratum)) {
    pos <- which(target_stratum == s)
    pool <- which(universe_stratum == s)
    selected[pos] <- sample(pool, length(pos), replace = FALSE)
  }
  selected
}

features_for_state <- function(state) cov_features_matrix(universe_z[state, , drop = FALSE], directions)

is_feasible <- function(features) {
  is.finite(features[["pc1"]]) && is.finite(features[["mean_r"]]) &&
    abs(features[["pc1"]] - target_cov[["pc1"]]) <= tolerance &&
    abs(features[["mean_r"]] - target_cov[["mean_r"]]) <= tolerance
}

objective <- function(features) {
  ((features[["pc1"]] - target_cov[["pc1"]]) / tolerance)^2 +
    ((features[["mean_r"]] - target_cov[["mean_r"]]) / tolerance)^2
}

propose_state <- function(state) {
  position <- sample.int(length(state), 1L)
  choices <- setdiff(position_pools[[position]], state)
  if (!length(choices)) return(NULL)
  proposal <- state
  proposal[[position]] <- sample(choices, 1L)
  proposal
}

anneal_chain <- function(chain_id) {
  set.seed(chain_seeds[[chain_id]])
  state <- make_initial_state()
  features <- features_for_state(state)
  current_obj <- objective(features)
  initial_obj <- current_obj
  accepted <- 0L
  first_feasible <- NA_integer_
  if (is_feasible(features)) first_feasible <- 0L

  if (is.na(first_feasible)) {
    initial_temp <- max(1, current_obj)
    for (step in seq_len(max_anneal_proposals)) {
      proposal <- propose_state(state)
      if (is.null(proposal)) next
      proposed_features <- features_for_state(proposal)
      proposed_obj <- objective(proposed_features)
      temperature <- initial_temp * (0.001 ^ (step / max_anneal_proposals))
      accept <- proposed_obj <= current_obj || runif(1) < exp((current_obj - proposed_obj) / temperature)
      if (accept) {
        state <- proposal
        features <- proposed_features
        current_obj <- proposed_obj
        accepted <- accepted + 1L
      }
      if (is_feasible(features)) {
        first_feasible <- step
        break
      }
    }
  }

  list(
    chain_id = chain_id,
    seed = chain_seeds[[chain_id]],
    feasible = !is.na(first_feasible),
    first_feasible_proposal = first_feasible,
    anneal_proposals = if (is.na(first_feasible)) max_anneal_proposals else first_feasible,
    anneal_accepted = accepted,
    initial_objective = initial_obj,
    state = state,
    features = features,
    within_proposals = 0L,
    feasible_accepted = 0L,
    accepted_since_last_retain = 0L,
    retained_total = 0L,
    rng_state = .Random.seed
  )
}

run_sampling_round <- function(chain_state, proposals_this_round) {
  assign(".Random.seed", chain_state$rng_state, envir = .GlobalEnv)
  new_records <- list()
  if (!chain_state$feasible) return(list(state = chain_state, records = new_records))

  for (step in seq_len(proposals_this_round)) {
    proposal <- propose_state(chain_state$state)
    chain_state$within_proposals <- chain_state$within_proposals + 1L
    if (is.null(proposal)) next
    proposed_features <- features_for_state(proposal)
    if (!is_feasible(proposed_features)) next
    chain_state$state <- proposal
    chain_state$features <- proposed_features
    chain_state$feasible_accepted <- chain_state$feasible_accepted + 1L

    if (chain_state$feasible_accepted <= burnin_accepted) next
    chain_state$accepted_since_last_retain <- chain_state$accepted_since_last_retain + 1L
    if (chain_state$accepted_since_last_retain >= thin_accepted) {
      chain_state$accepted_since_last_retain <- 0L
      chain_state$retained_total <- chain_state$retained_total + 1L
      new_records[[length(new_records) + 1L]] <- list(
        chain_id = chain_state$chain_id,
        chain_retain_index = chain_state$retained_total,
        state = chain_state$state,
        pc1 = unname(chain_state$features[["pc1"]]),
        mean_r = unname(chain_state$features[["mean_r"]])
      )
    }
  }
  chain_state$rng_state <- .Random.seed
  list(state = chain_state, records = new_records)
}

chains <- mclapply(seq_len(n_chains), anneal_chain, mc.cores = workers, mc.preschedule = FALSE)
if (any(vapply(chains, inherits, logical(1), "try-error"))) stop("One or more annealing chains failed")
feasible_chain_count <- sum(vapply(chains, function(x) x$feasible, logical(1)))

all_records <- list()
unique_records <- list()
if (feasible_chain_count >= minimum_contributing_chains) {
  rounds <- ceiling(max_within_proposals / round_proposals)
  for (round_id in seq_len(rounds)) {
    remaining <- max_within_proposals - vapply(chains, function(x) x$within_proposals, integer(1))
    n_this <- pmin(round_proposals, pmax(0L, remaining))
    active <- which(vapply(chains, function(x) x$feasible, logical(1)) & n_this > 0L)
    if (!length(active)) break
    round_results <- mclapply(
      active,
      function(i) run_sampling_round(chains[[i]], n_this[[i]]),
      mc.cores = min(workers, length(active)),
      mc.preschedule = FALSE
    )
    if (any(vapply(round_results, inherits, logical(1), "try-error"))) stop("A sampling chain failed")
    for (j in seq_along(active)) {
      chains[[active[[j]]]] <- round_results[[j]]$state
      all_records <- c(all_records, round_results[[j]]$records)
    }

    if (length(all_records)) {
      keys <- vapply(all_records, function(rec) paste(sort(universe_genes[rec$state]), collapse = ";"), character(1))
      unique_records <- all_records[!duplicated(keys)]
      contributing <- length(unique(vapply(unique_records, function(rec) rec$chain_id, integer(1))))
      if (length(unique_records) >= target_unique && contributing >= minimum_contributing_chains) break
    }
  }
}

if (length(all_records)) {
  all_keys <- vapply(all_records, function(rec) paste(sort(universe_genes[rec$state]), collapse = ";"), character(1))
  unique_records <- all_records[!duplicated(all_keys)]
} else {
  all_keys <- character()
  unique_records <- list()
}

unique_contributing_chains <- if (length(unique_records)) {
  length(unique(vapply(unique_records, function(rec) rec$chain_id, integer(1))))
} else 0L

tier_estimable <- length(unique_records) >= target_unique && unique_contributing_chains >= minimum_contributing_chains
status <- if (tier_estimable) "estimable_at_0.05" else "not_estimable_at_0.05"
reason_code <- if (tier_estimable) {
  "success_1000_unique_5plus_chains"
} else if (feasible_chain_count < minimum_contributing_chains) {
  "fewer_than_5_feasible_chains"
} else if (unique_contributing_chains < minimum_contributing_chains) {
  "fewer_than_5_contributing_chains"
} else {
  "fewer_than_1000_unique_sets"
}

# Deterministic round-robin selection limits the final scoring manifest to exactly 1,000 sets.
selected_records <- list()
if (tier_estimable) {
  ord <- order(
    vapply(unique_records, function(rec) rec$chain_retain_index, integer(1)),
    vapply(unique_records, function(rec) rec$chain_id, integer(1))
  )
  selected_records <- unique_records[ord[seq_len(target_unique)]]
}

records_to_summary <- function(records, prefix) {
  if (!length(records)) return(data.table())
  rbindlist(lapply(seq_along(records), function(i) {
    rec <- records[[i]]
    data.table(
      set_id = sprintf("%s_%04d", prefix, i),
      chain_id = rec$chain_id,
      chain_retain_index = rec$chain_retain_index,
      pc1 = rec$pc1,
      mean_r = rec$mean_r,
      abs_delta_pc1 = abs(rec$pc1 - target_cov[["pc1"]]),
      abs_delta_mean_r = abs(rec$mean_r - target_cov[["mean_r"]]),
      gene_key = paste(sort(universe_genes[rec$state]), collapse = ";")
    )
  }))
}

records_to_long <- function(records, prefix) {
  if (!length(records)) return(data.table())
  rbindlist(lapply(seq_along(records), function(i) {
    rec <- records[[i]]
    data.table(
      set_id = sprintf("%s_%04d", prefix, i),
      chain_id = rec$chain_id,
      position = seq_along(rec$state),
      target_position_gene = target$gene,
      gene = universe_genes[rec$state],
      direction = directions,
      mean_decile = t_mean_decile,
      sd_decile = t_sd_decile,
      joint_stratum = target_stratum
    )
  }))
}

unique_summary <- records_to_summary(unique_records, "candidate")
unique_long <- records_to_long(unique_records, "candidate")
selected_summary <- records_to_summary(selected_records, "null")
selected_long <- records_to_long(selected_records, "null")

fwrite(unique_summary, file.path(out_dir, "retained_unique_sets_summary.tsv"), sep = "\t")
fwrite(unique_long, file.path(out_dir, "retained_unique_gene_sets.tsv"), sep = "\t")
fwrite(selected_summary, file.path(out_dir, "selected_null_sets_summary.tsv"), sep = "\t")
fwrite(selected_long, file.path(out_dir, "selected_null_gene_sets.tsv"), sep = "\t")

selected_chain_counts <- if (nrow(selected_summary)) selected_summary[, .N, by = chain_id] else data.table(chain_id = integer(), N = integer())
chain_diagnostics <- rbindlist(lapply(chains, function(x) {
  data.table(
    chain_id = x$chain_id,
    seed = x$seed,
    entered_feasible_region = x$feasible,
    first_feasible_proposal = x$first_feasible_proposal,
    anneal_proposals = x$anneal_proposals,
    anneal_accepted = x$anneal_accepted,
    initial_objective = x$initial_objective,
    within_region_proposals = x$within_proposals,
    accepted_feasible_swaps = x$feasible_accepted,
    within_region_acceptance_rate = if (x$within_proposals) x$feasible_accepted / x$within_proposals else NA_real_,
    retained_raw = x$retained_total
  )
}))
chain_diagnostics <- merge(chain_diagnostics, selected_chain_counts, by = "chain_id", all.x = TRUE)
setnames(chain_diagnostics, "N", "selected_unique_sets")
chain_diagnostics[is.na(selected_unique_sets), selected_unique_sets := 0L]
fwrite(chain_diagnostics, file.path(out_dir, "chain_diagnostics.tsv"), sep = "\t")

target_positions <- data.table(
  position = seq_len(nrow(target)),
  target_gene = target$gene,
  direction = directions,
  expression_mean = unname(t_mean),
  expression_sd = unname(t_sd),
  mean_decile = t_mean_decile,
  sd_decile = t_sd_decile,
  joint_stratum = target_stratum,
  candidate_pool_size = pool_sizes
)
fwrite(target_positions, file.path(out_dir, "target_position_strata.tsv"), sep = "\t")

jaccard_values <- numeric()
jaccard_records <- if (length(selected_records)) selected_records else unique_records
if (length(jaccard_records) >= 2L) {
  n_pair <- length(jaccard_records) * (length(jaccard_records) - 1L) / 2L
  jaccard_values <- numeric(n_pair)
  cursor <- 1L
  gene_sets <- lapply(jaccard_records, function(rec) universe_genes[rec$state])
  for (i in seq_len(length(gene_sets) - 1L)) {
    for (j in (i + 1L):length(gene_sets)) {
      a <- gene_sets[[i]]
      b <- gene_sets[[j]]
      jaccard_values[[cursor]] <- length(intersect(a, b)) / length(union(a, b))
      cursor <- cursor + 1L
    }
  }
}

jaccard_summary <- data.table(
  n_sets = length(jaccard_records),
  n_pairs = length(jaccard_values),
  minimum = if (length(jaccard_values)) min(jaccard_values) else NA_real_,
  q25 = if (length(jaccard_values)) unname(quantile(jaccard_values, 0.25)) else NA_real_,
  median = if (length(jaccard_values)) median(jaccard_values) else NA_real_,
  q75 = if (length(jaccard_values)) unname(quantile(jaccard_values, 0.75)) else NA_real_,
  maximum = if (length(jaccard_values)) max(jaccard_values) else NA_real_
)
fwrite(jaccard_summary, file.path(out_dir, "jaccard_summary.tsv"), sep = "\t")

seed_manifest <- data.table(chain_id = seq_len(n_chains), seed = chain_seeds)
fwrite(seed_manifest, file.path(out_dir, "seed_manifest.tsv"), sep = "\t")

input_access <- data.table(
  input_role = c("expression_features", "frozen_program_definitions"),
  path = normalizePath(input_paths),
  sha256 = vapply(input_paths, sha256_file, character(1)),
  contains_clinical_labels = FALSE
)
fwrite(input_access, file.path(out_dir, "input_access_manifest.tsv"), sep = "\t")

selected_gene_file <- file.path(out_dir, "selected_null_gene_sets.tsv")
retained_gene_file <- file.path(out_dir, "retained_unique_gene_sets.tsv")
script_hash <- if (!is.na(script_path) && file.exists(script_path)) sha256_file(script_path) else NA_character_
gene_manifest_for_hash <- if (tier_estimable) selected_gene_file else retained_gene_file

tier_manifest <- data.table(
  cohort = "GSE151179",
  program = args$program,
  tier = tolerance,
  tolerance_pc1 = tolerance,
  tolerance_mean_r = tolerance,
  status = status,
  chains_attempted = n_chains,
  chains_feasible = feasible_chain_count,
  chains_contributing = unique_contributing_chains,
  annealing_proposals = sum(chain_diagnostics$anneal_proposals),
  within_region_proposals = sum(chain_diagnostics$within_region_proposals),
  accepted_swaps = sum(chain_diagnostics$accepted_feasible_swaps),
  retained_raw = length(all_records),
  unique_sets = length(unique_records),
  selected_sets = length(selected_records),
  selected_for_scoring = tier_estimable,
  reason_code = reason_code,
  seed_manifest = file.path(out_dir, "seed_manifest.tsv"),
  gene_set_manifest_sha256 = sha256_file(gene_manifest_for_hash),
  script_sha256 = script_hash
)
fwrite(tier_manifest, file.path(out_dir, "tier_manifest.tsv"), sep = "\t")

cov_source <- if (nrow(selected_summary)) selected_summary else unique_summary
covariance_diagnostic <- data.table(
  feature = c("pc1_variance_proportion", "mean_pairwise_pearson_r"),
  target = unname(target_cov),
  null_median = c(
    if (nrow(cov_source)) median(cov_source$pc1) else NA_real_,
    if (nrow(cov_source)) median(cov_source$mean_r) else NA_real_
  ),
  null_q025 = c(
    if (nrow(cov_source)) unname(quantile(cov_source$pc1, 0.025)) else NA_real_,
    if (nrow(cov_source)) unname(quantile(cov_source$mean_r, 0.025)) else NA_real_
  ),
  null_q975 = c(
    if (nrow(cov_source)) unname(quantile(cov_source$pc1, 0.975)) else NA_real_,
    if (nrow(cov_source)) unname(quantile(cov_source$mean_r, 0.975)) else NA_real_
  ),
  maximum_absolute_deviation = c(
    if (nrow(cov_source)) max(cov_source$abs_delta_pc1) else NA_real_,
    if (nrow(cov_source)) max(cov_source$abs_delta_mean_r) else NA_real_
  ),
  tolerance = tolerance
)
fwrite(covariance_diagnostic, file.path(out_dir, "covariance_diagnostic.tsv"), sep = "\t")

label_blind_audit <- list(
  passed = TRUE,
  clinical_label_files_loaded = FALSE,
  accepted_input_roles = c("expression_features", "frozen_program_definitions"),
  expression_columns_are_geo_accessions_only = TRUE,
  n_expression_samples = ncol(expr),
  null_generation_completed_before_endpoint_scoring = TRUE,
  endpoint_scoring_performed = FALSE
)
write_json(label_blind_audit, file.path(out_dir, "label_blind_audit.json"), pretty = TRUE, auto_unbox = TRUE)

diagnostic_summary <- list(
  cohort = "GSE151179",
  program = args$program,
  tolerance = tolerance,
  status = status,
  reason_code = reason_code,
  target = list(
    gene_count = nrow(target),
    positive_directions = sum(directions == 1L),
    negative_directions = sum(directions == -1L),
    pc1_variance_proportion = unname(target_cov[["pc1"]]),
    mean_pairwise_pearson_r = unname(target_cov[["mean_r"]])
  ),
  universe_gene_count = length(universe_genes),
  chains = list(
    attempted = n_chains,
    feasible = feasible_chain_count,
    contributing = unique_contributing_chains
  ),
  sets = list(
    retained_raw = length(all_records),
    retained_unique = length(unique_records),
    selected = length(selected_records),
    duplicate_rate = if (length(all_records)) 1 - length(unique_records) / length(all_records) else NA_real_
  ),
  jaccard = as.list(jaccard_summary[1]),
  computational_budget = list(
    max_anneal_proposals_per_chain = max_anneal_proposals,
    burnin_accepted_swaps = burnin_accepted,
    thinning_accepted_swaps = thin_accepted,
    max_within_region_proposals_per_chain = max_within_proposals,
    workers = workers
  ),
  label_blind = label_blind_audit
)
write_json(diagnostic_summary, file.path(out_dir, "diagnostic_summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = 15)

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

files_to_hash <- sort(list.files(out_dir, full.names = TRUE))
files_to_hash <- files_to_hash[basename(files_to_hash) != "FROZEN_SHA256.tsv"]
hash_manifest <- data.table(
  file = basename(files_to_hash),
  sha256 = vapply(files_to_hash, sha256_file, character(1))
)
fwrite(hash_manifest, file.path(out_dir, "FROZEN_SHA256.tsv"), sep = "\t")

cat(toJSON(list(
  status = status,
  reason_code = reason_code,
  target_pc1 = unname(target_cov[["pc1"]]),
  target_mean_r = unname(target_cov[["mean_r"]]),
  chains_feasible = feasible_chain_count,
  chains_contributing = unique_contributing_chains,
  retained_raw = length(all_records),
  unique_sets = length(unique_records),
  selected_sets = length(selected_records),
  out_dir = out_dir
), auto_unbox = TRUE, pretty = TRUE, digits = 15), "\n")
