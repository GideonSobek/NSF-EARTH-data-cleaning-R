# ----  Read the file ----
lines <- readLines("raw_data2.tsv")
rows <- strsplit(lines, "\t")
length(rows)

# identify headers and nois rows
is_header <- function(row) {
  length(row) >= 4 && row[3] == "City" && row[4] == "State"
}

is_noise <- function(row) {
  if (length(row) < 3) return(TRUE)
  city_field <- row[3]
  is.na(city_field) || city_field == "" || startsWith(city_field, "###")
}

#  Find real headers 
first_detailed_idx <- which(sapply(rows, function(r) is_header(r) && "Mechanical_Code_Model" %in% r))[1]
first_scored_idx   <- which(sapply(rows, function(r) is_header(r) && !("Mechanical_Code_Model" %in% r)))[1]

detailed_schema <- rows[[first_detailed_idx]]
scored_schema   <- rows[[first_scored_idx]]

is_header(detailed_schema)
is_header(scored_schema)

# Main loop using Type column instead, stores in lists
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
  } else if (row[2] == "score") {
    if (length(row) >= length(scored_schema)) {
      named_row <- row[1:length(scored_schema)]
      names(named_row) <- scored_schema
      scored_records[[length(scored_records) + 1]] <- named_row
    }
  } else if (row[2] == "detail") {
    if (length(row) >= length(detailed_schema)) {
      named_row <- row[1:length(detailed_schema)]
      names(named_row) <- detailed_schema
      detailed_records[[length(detailed_records) + 1]] <- named_row
    }
  }
}

length(detailed_records)
length(scored_records)

# Build data frames from lists
detailed_df <- as.data.frame(do.call(rbind, detailed_records), stringsAsFactors = FALSE)
scored_df <- as.data.frame(do.call(rbind, scored_records), stringsAsFactors = FALSE)

colnames(detailed_df)[colnames(detailed_df) == "Chat"] <- "Rater"
colnames(scored_df)[colnames(scored_df) == "Chat"] <- "Rater"

#  Sort
detailed_df <- detailed_df[order(detailed_df$City), ]
scored_df <- scored_df[order(scored_df$City), ]

# Check for duplicates
dup_check <- function(df) {
  keys <- paste(df$Rater, df$City, df$State)
  df[duplicated(keys) | duplicated(keys, fromLast = TRUE), ]
}

dup_check(detailed_df)
dup_check(scored_df)

# Coverage check
setdiff(paste(detailed_df$Rater, detailed_df$City), paste(scored_df$Rater, scored_df$City))
setdiff(paste(scored_df$Rater, scored_df$City), paste(detailed_df$Rater, detailed_df$City))

# Export 
write.table(detailed_df, "detailed_df_export.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(scored_df, "scored_df_export.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# Reliability statistics (Krippendorf, Gwet, normals Stats, )


#install.packages(c("irr", "irrCAC"))  # one-time
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

# Build a City x Rater matrix for one dimension
build_unit_matrix <- function(df, col) {
  raters <- sort(unique(df$Rater))
  cities <- sort(unique(df$City))
  mat <- matrix(NA, nrow = length(cities), ncol = length(raters),
                dimnames = list(cities, raters))
  for (i in 1:nrow(df)) {
    mat[df$City[i], df$Rater[i]] <- df[[col]][i]
  }
  mat
}

# Which method applies to which column
# ordinal:  1-3 restrictiveness/score dimensions (rank distance matters,
#          treats a 1-vs-3 gap as bigger than a 1-vs-2 gap)
# interval: years (squared numeric distance)
# nominal: free-text status categories
# this way the stats can apply to scores and to years 
dimension_methods <- list(
  Mech_Restrictiveness   = "ordinal",
  Fire_Restrictiveness   = "ordinal",
  A2L_Pathway             = "ordinal",
  Mechanic_Score          = "ordinal",
  Operator_Score          = "ordinal",
  Training_Score          = "ordinal",
  BEPS_Strength           = "ordinal",
  Gas_Restriction         = "ordinal",
  Decarb_Plan             = "ordinal",
  Mech_Code_Year          = "interval",
  Fire_Code_Year          = "interval",
  Gas_Restriction_Status  = "nominal"
)

compute_alpha <- function(df, col, method) {
  mat <- build_unit_matrix(df, col)
  kripp.alpha(t(mat), method = method)   # irr wants raters as rows, so transpose
}

# ---- 4. Run Krippendorff's alpha across every dimension ----
alpha_results <- data.frame(Dimension = character(), Alpha = numeric())
for (dim in names(dimension_methods)) {
  res <- compute_alpha(scored_df, dim, dimension_methods[[dim]])
  alpha_results <- rbind(alpha_results, data.frame(Dimension = dim, Alpha = round(res$value, 3)))
}
alpha_results

# ---- 5. Gwet's AC2 — the standard alternative to Krippendorff's alpha ----
# Less prone to the "high raw agreement but low alpha" paradox when one
# category dominates. Ordinal/interval weighting matches kripp.alpha's method.

#ensure that Gwet handles missing data like Kripp
compute_gwet <- function(df, col, weighting) {
  df_clean <- df[!is.na(df[[col]]), ]
  mat <- build_unit_matrix(df_clean, col)
  gwet.ac1.raw(mat, weights = weighting)$est$coeff.val
}

gwet_results <- data.frame(Dimension = character(), Gwet_AC2 = numeric())
for (dim in names(dimension_methods)) {
  method <- dimension_methods[[dim]]
  w <- switch(method,
              "nominal"  = "unweighted",
              "ordinal"  = "ordinal",
              "interval" = "quadratic")   # correct interval-equivalent weighting
  val <- tryCatch(compute_gwet(scored_df, dim, w), error = function(e) NA)
  gwet_results <- rbind(gwet_results, data.frame(Dimension = dim, Gwet_AC2 = round(val, 3)))
}
gwet_results

# Leave-one-out (jackknife) alpha — flags influential cities
jackknife_alpha <- function(df, col, method) {
  cities <- sort(unique(df$City))
  full_alpha <- compute_alpha(df, col, method)$value
  out <- data.frame(City = character(), Alpha_Without = numeric(), Delta = numeric())
  for (city in cities) {
    sub <- df[df$City != city, ]
    a <- tryCatch(compute_alpha(sub, col, method)$value, error = function(e) NA)
    out <- rbind(out, data.frame(City = city, Alpha_Without = round(a, 3), Delta = round(a - full_alpha, 3)))
  }
  out[order(-abs(out$Delta)), ]
}

# example usage — swap in any dimension name, this prints out
# each city with a "without" alpha, and delta from reg. alpha
# jackknife_alpha(scored_df, "Gas_Restriction", "ordinal")

# Per-unit raw disagreement (non-chance-corrected diagnostic) 
# Mean pairwise |difference| among raters within a city — cheap first pass
# before running the full jackknife on a suspect dimension.
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

unit_disagreement(scored_df, "Gas_Restriction")

# Run jackknife across every dimension, combine into one table
jackknife_all <- data.frame()
for (dim in names(dimension_methods)) {
  res <- jackknife_alpha(scored_df, dim, dimension_methods[[dim]])
  res$Dimension <- dim
  jackknife_all <- rbind(jackknife_all, res)
}

# Reorder columns so Dimension comes first
jackknife_all <- jackknife_all[, c("Dimension", "City", "Alpha_Without", "Delta")]

write.xlsx(
  list(
    Scored_Data = scored_df,
    Detailed_Data = detailed_df,
    Alpha_Results = alpha_results,
    Gwet_Results = gwet_results,
    Jackknife = jackknife_all
  ),
  file = "hurl_carp_full_results.xlsx"
)