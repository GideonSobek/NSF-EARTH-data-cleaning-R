# ============================================================
# NSF ERC EATH Project 1.4 (Thrust 1) 
# Gideon Sobek RA, Todd Watkins PI
# August 2026
# ============================================================


# SECTION 1: Parse with FIPS
# ============================================================

lines <- readLines("master2.tsv", warn = FALSE)

rows <- lapply(lines, function(x) {
  parts <- strsplit(paste0(x, "\tEND"), "\t", fixed = TRUE)[[1]]
  parts[-length(parts)]
})
length(rows)

is_header <- function(row) {
  length(row) >= 3 && row[2] == "City" && row[3] == "State"
}

is_noise <- function(row) {
  if (length(row) < 2) return(TRUE)
  city_field <- row[2]
  is.na(city_field) || city_field == "" || startsWith(city_field, "###")
}

first_detailed_idx <- which(sapply(rows, function(r) is_header(r) && "Mechanical_Code_Model" %in% r))[1]
first_scored_idx   <- which(sapply(rows, function(r) is_header(r) && !("Mechanical_Code_Model" %in% r)))[1]

detailed_schema <- rows[[first_detailed_idx]]
scored_schema   <- rows[[first_scored_idx]]

is_header(detailed_schema)
is_header(scored_schema)

# Confirm how many "notes" rows exist — these are duplicates of "detail"
# rows and are intentionally skipped by the loop below (row[5] == "detail"
# only, never "notes")
length(which(sapply(rows, function(r) length(r) >= 5 && r[5] == "notes")))

detailed_records <- list()
scored_records <- list()

for (i in 1:length(rows)) {
  row <- rows[[i]]
  
  if (is_header(row)) {
    if ("Mechanical_Code_Model" %in% row) {
      detailed_schema <- row
    } else {
      scored_schema <- row
    }
  } else if (is_noise(row)) {
    # skip
  } else if (length(row) >= 5 && tolower(trimws(row[5])) == "score") {
    
    # Pad short rows with blank cells if necessary
    if (length(row) < length(scored_schema)) {
      row <- c(row, rep("", length(scored_schema) - length(row)))
    }
    
    row[4] <- trimws(row[4])  # Rater
    named_row <- row[1:length(scored_schema)]
    names(named_row) <- scored_schema
    scored_records[[length(scored_records) + 1]] <- named_row
    
  } else if (length(row) >= 5 && tolower(trimws(row[5])) == "detail") {
    
    # Pad short rows with blank cells if necessary
    if (length(row) < length(detailed_schema)) {
      row <- c(row, rep("", length(detailed_schema) - length(row)))
    }
    
    row[4] <- trimws(row[4])  # Rater
    named_row <- row[1:length(detailed_schema)]
    names(named_row) <- detailed_schema
    detailed_records[[length(detailed_records) + 1]] <- named_row
  }
  # "notes" rows (duplicates of "detail") and anything else fall through
  # untouched — intentional, not a bug
}

length(detailed_records)
length(scored_records)




detailed_df <- as.data.frame(do.call(rbind, detailed_records), stringsAsFactors = FALSE)
scored_df <- as.data.frame(do.call(rbind, scored_records), stringsAsFactors = FALSE)

# Rename column 4 by POSITION, not name — scored header has a "Model" typo
colnames(detailed_df)[4] <- "Rater"
colnames(scored_df)[4] <- "Rater"

# Sort by FIPS now, not City name
detailed_df <- detailed_df[order(detailed_df$FIPS), ]
scored_df <- scored_df[order(scored_df$FIPS), ]

# Duplicate check — keyed on FIPS + Rater, not City + Rater
dup_check <- function(df) {
  keys <- paste(df$Rater, df$FIPS)
  df[duplicated(keys) | duplicated(keys, fromLast = TRUE), ]
}
dup_check(detailed_df)
dup_check(scored_df)

# Coverage check — FIPS-based
setdiff(paste(detailed_df$Rater, detailed_df$FIPS), paste(scored_df$Rater, scored_df$FIPS))
setdiff(paste(scored_df$Rater, scored_df$FIPS), paste(detailed_df$Rater, detailed_df$FIPS))

write.table(detailed_df, "detailed_df_export.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(scored_df, "scored_df_export.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# SECTION 2: Reliability analysis (FIPS-based units)
# ============================================================

library(irr)
library(irrCAC)

score_cols <- c("Mech_Restrictiveness","Fire_Restrictiveness","A2L_Pathway",
                "Mechanic_Score","Operator_Score","Training_Score",
                "BEPS_Strength","Gas_Restriction","Decarb_Plan")

for (col in score_cols) {
  scored_df[[col]] <- suppressWarnings(as.numeric(scored_df[[col]]))
}
scored_df$Mech_Code_Year <- suppressWarnings(as.numeric(scored_df$Mech_Code_Year))
scored_df$Fire_Code_Year <- suppressWarnings(as.numeric(scored_df$Fire_Code_Year))

# Units are now FIPS codes, not city names — no more collision risk
build_unit_matrix <- function(df, col) {
  raters <- sort(unique(df$Rater))
  units  <- sort(unique(df$FIPS))
  mat <- matrix(NA, nrow = length(units), ncol = length(raters),
                dimnames = list(units, raters))
  for (i in 1:nrow(df)) {
    mat[df$FIPS[i], df$Rater[i]] <- df[[col]][i]
  }
  mat
}

# Small reference table so FIPS-keyed output stays human-readable
fips_lookup <- unique(scored_df[, c("FIPS", "City", "State")])

dimension_methods <- list(
  Mech_Restrictiveness   = "ordinal", Fire_Restrictiveness = "ordinal",
  A2L_Pathway             = "ordinal", Mechanic_Score        = "ordinal",
  Operator_Score          = "ordinal", Training_Score        = "ordinal",
  BEPS_Strength           = "ordinal", Gas_Restriction       = "ordinal",
  Decarb_Plan             = "ordinal", Mech_Code_Year        = "interval",
  Fire_Code_Year          = "interval", Gas_Restriction_Status = "nominal"
)

compute_alpha <- function(df, col, method) {
  mat <- build_unit_matrix(df, col)
  kripp.alpha(t(mat), method = method)
}

alpha_results <- data.frame(Dimension = character(), Alpha = numeric())
for (dim in names(dimension_methods)) {
  res <- compute_alpha(scored_df, dim, dimension_methods[[dim]])
  alpha_results <- rbind(alpha_results, data.frame(Dimension = dim, Alpha = round(res$value, 3)))
}
alpha_results

compute_gwet <- function(df, col, weighting) {
  mat <- build_unit_matrix(df, col)
  gwet.ac1.raw(mat, weights = weighting)$est$coeff.val
}

gwet_results <- data.frame(Dimension = character(), Gwet_AC2 = numeric())
for (dim in names(dimension_methods)) {
  w <- if (dimension_methods[[dim]] == "nominal") "unweighted" else "ordinal"
  val <- tryCatch(compute_gwet(scored_df, dim, w), error = function(e) NA)
  gwet_results <- rbind(gwet_results, data.frame(Dimension = dim, Gwet_AC2 = round(val, 3)))
}
gwet_results

# Jackknife — now drops one FIPS unit at a time, joins City/State back in for readability
jackknife_alpha <- function(df, col, method) {
  units <- sort(unique(df$FIPS))
  full_alpha <- compute_alpha(df, col, method)$value
  out <- data.frame(FIPS = character(), Alpha_Without = numeric(), Delta = numeric())
  for (u in units) {
    sub <- df[df$FIPS != u, ]
    a <- tryCatch(compute_alpha(sub, col, method)$value, error = function(e) NA)
    out <- rbind(out, data.frame(FIPS = u, Alpha_Without = round(a, 3), Delta = round(a - full_alpha, 3)))
  }
  out <- merge(out, fips_lookup, by = "FIPS")
  out[order(-abs(out$Delta)), c("FIPS", "City", "State", "Alpha_Without", "Delta")]
}

jackknife_all <- data.frame()
for (dim in names(dimension_methods)) {
  res <- jackknife_alpha(scored_df, dim, dimension_methods[[dim]])
  res$Dimension <- dim
  jackknife_all <- rbind(jackknife_all, res)
}
jackknife_all <- jackknife_all[, c("Dimension", "FIPS", "City", "State", "Alpha_Without", "Delta")]

unit_disagreement <- function(df, col) {
  mat <- build_unit_matrix(df, col)
  vals <- apply(mat, 1, function(row) {
    v <- row[!is.na(row)]
    if (length(v) < 2) return(NA)
    pairs <- combn(v, 2)
    mean(abs(pairs[1, ] - pairs[2, ]))
  })
  sort(vals, decreasing = TRUE)
}

# SECTION 3: Sampling stratum + stratified agreement stats
# ============================================================

# --- Stratum lookup (from stratified sampling design) ---
# Ambition tier = your High/Medium/Low ambition table, mapped onto FIPS
# already present in scored_df. St. Louis and Augusta had "?" ambition
# scores in the source table; kept in their listed tier with
# Ambition_Score = NA rather than guessing a number.

stratum_lookup <- data.frame(
  FIPS = c(
    # High Ambition
    644000, 1150000, 4159000, 2507000, 3651000, 666000, 4260000, 667000, 668000, 5363000,
    # Medium Ambition
    1714000, 3436000, 2965000, 2360545, 3755000, 1077580, 937000, 1571550, 5010675, 5167000,
    # Low Ambition
    4865000, 5613900, 2938000, 4962470, 3918000, 203000, 1772000, 1836003, 1304204, 3240000
  ),
  Ambition_Tier = c(
    rep("High",   10),
    rep("Medium", 10),
    rep("Low",    10)
  ),
  Ambition_Score = c(
    3.00, 2.67, 2.67, 2.67, 3.00, 2.33, 2.33, 2.67, 2.00, 3.00,
    1.67, 1.67, NA,   1.67, 1.67, 1.67, 1.67, 1.67, 1.67, 1.67,
    1.33, 1.00, 1.33, 1.00, 1.33, 1.00, 1.33, 1.33, NA,   1.33
  ),
  Region = c(
    "Pacific", "Mid-Atlantic", "Pacific", "New England", "NE/MA", "Pacific", "Mid-Atlantic", "Pacific", "Pacific", "Pacific",
    "Midwest", "Northeast", "Midwest", "New England", "Southeast", "Pacific NW", "NE", "Pacific", "NE", "Mid-Atlantic",
    "Texas", "Mountain", "Midwest", "Mountain", "Ohio", "Pacific/Remote", "East North Central", "Midwest", "Southeast", "Mountain West"
  ),
  Population = c(
    3878704, 702250, 635749, 673458, 8478072, 1404452, 1573916, 827526, 997368, 780995,
    2721308, 302824, NA, 70543, 499825, 73512, 121981, 344967, 44019, 233655,
    1526656, 66507, 516032, 114527, 933263, 289600, 112989, 891484, NA, 678922
  ),
  stringsAsFactors = FALSE
)
stratum_lookup$FIPS <- as.character(stratum_lookup$FIPS)

# Sanity check: every FIPS in scored_df should have a stratum, and vice versa.
# If either of these is non-empty, a city name/FIPS in the lookup table
# doesn't match what's in scored_df -- fix before trusting results below.
setdiff(unique(scored_df$FIPS), stratum_lookup$FIPS)
setdiff(stratum_lookup$FIPS, unique(scored_df$FIPS))

# Join stratum info onto scored_df
scored_df <- merge(scored_df, stratum_lookup, by = "FIPS", all.x = TRUE)
scored_df <- scored_df[order(scored_df$FIPS), ]

# --- Agreement-rate helper functions (unit = one FIPS city) ---
# All three treat a "unit" as evaluateable only if >=2 raters have a
# non-NA value for that dimension; units with 0-1 ratings are excluded
# and don't count toward the denominator (n_units reflects this).

full_agreement_rate <- function(df, col) {
  mat <- build_unit_matrix(df, col)
  agree <- apply(mat, 1, function(row) {
    v <- row[!is.na(row)]
    if (length(v) < 2) return(NA)
    length(unique(v)) == 1
  })
  agree <- agree[!is.na(agree)]
  c(rate = mean(agree), n_units = length(agree))
}

majority_agreement_rate <- function(df, col) {
  mat <- build_unit_matrix(df, col)
  agree <- apply(mat, 1, function(row) {
    v <- row[!is.na(row)]
    if (length(v) < 2) return(NA)
    tab <- table(v)
    max(tab) >= ceiling(length(v) / 2)
  })
  agree <- agree[!is.na(agree)]
  c(rate = mean(agree), n_units = length(agree))
}

mean_range <- function(df, col) {
  mat <- build_unit_matrix(df, col)
  ranges <- apply(mat, 1, function(row) {
    v <- row[!is.na(row)]
    if (length(v) < 2) return(NA)
    max(v) - min(v)
  })
  ranges <- ranges[!is.na(ranges)]
  c(mean_range = mean(ranges), n_units = length(ranges))
}

# --- Stratified alpha (reuses existing compute_alpha + dimension_methods) ---
stratified_alpha <- function(df, col, method, stratum = NULL) {
  sub <- if (is.null(stratum)) df else df[df$Ambition_Tier == stratum, ]
  if (nrow(sub) == 0) return(NA)
  tryCatch(compute_alpha(sub, col, method)$value, error = function(e) NA)
}

# --- Stratified Gwet (fixed: compute_gwet already returns a scalar) ---
stratified_gwet <- function(df, col, weighting, stratum = NULL) {
  sub <- if (is.null(stratum)) df else df[df$Ambition_Tier == stratum, ]
  if (nrow(sub) == 0) return(NA)
  tryCatch(compute_gwet(sub, col, weighting), error = function(e) NA)
}

# --- Build the full stratified summary table: dimension x stratum ---
strata_levels <- c("High", "Medium", "Low")

agreement_summary <- data.frame()

for (dim in score_cols) {
  method <- dimension_methods[[dim]]
  weighting <- if (method == "nominal") "unweighted" else "ordinal"
  
  # Overall (all strata combined)
  fa <- full_agreement_rate(scored_df, dim)
  ma <- majority_agreement_rate(scored_df, dim)
  mr <- mean_range(scored_df, dim)
  a  <- stratified_alpha(scored_df, dim, method, stratum = NULL)
  g  <- stratified_gwet(scored_df, dim, weighting, stratum = NULL)
  agreement_summary <- rbind(agreement_summary, data.frame(
    Dimension = dim, Stratum = "Overall",
    Full_Agreement_Rate = round(fa["rate"], 3),
    Majority_Agreement_Rate = round(ma["rate"], 3),
    Mean_Range = round(mr["mean_range"], 3),
    Alpha = round(a, 3),
    Gwet = round(g, 3),
    N_Units = fa["n_units"]
  ))
  
  # Per stratum
  for (s in strata_levels) {
    sub <- scored_df[scored_df$Ambition_Tier == s, ]
    fa <- full_agreement_rate(sub, dim)
    ma <- majority_agreement_rate(sub, dim)
    mr <- mean_range(sub, dim)
    a  <- stratified_alpha(scored_df, dim, method, stratum = s)
    g  <- stratified_gwet(scored_df, dim, weighting, stratum = s)
    agreement_summary <- rbind(agreement_summary, data.frame(
      Dimension = dim, Stratum = s,
      Full_Agreement_Rate = round(fa["rate"], 3),
      Majority_Agreement_Rate = round(ma["rate"], 3),
      Mean_Range = round(mr["mean_range"], 3),
      Alpha = round(a, 3),
      Gwet = round(g, 3),
      N_Units = fa["n_units"]
    ))
  }
}

rownames(agreement_summary) <- NULL
agreement_summary

# --- Add to the xlsx export (extends your existing write.xlsx call) ---
write.xlsx(
  list(
    Scored_Data = scored_df,
    Detailed_Data = detailed_df,
    Alpha_Results = alpha_results,
    Gwet_Results = gwet_results,
    Jackknife = jackknife_all,
    FIPS_Lookup = fips_lookup,
    Agreement_Summary = agreement_summary,
    Stratum_Lookup = stratum_lookup
  ),
  file = "hurl_carp_full_results.xlsx"
)
