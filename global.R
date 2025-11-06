suppressPackageStartupMessages({
  library(shiny); library(httr); library(dplyr); library(plotly); library(DT)
  library(msa); library(Biostrings); library(shinyjs); library(shinycssloaders)
  library(writexl); library(data.table)
})

combined_db_path <- "mhc_db/combined_db.rds"
IEDB_API_URL     <- "http://tools-cluster-interface.iedb.org/tools_api/mhci/"

sample_fasta <- tryCatch({
  readAAStringSet("mhc_db/sample_data.fasta")
}, error = function(e) {
  message("Warning: Could not load 'mhc_db/sample_data.fasta'. Sample data will not be available.")
  message(e)
  AAStringSet()
})

`%||%` <- function(a,b) if(!is.null(a)) a else b
is_nonempty_df <- function(x) is.data.frame(x) && nrow(x) > 0

get_consensus_vec <- function(msa_obj){
  aln <- as.matrix(msa_obj)
  apply(aln, 2, function(col){
    aa <- col[col != "-"]
    if(!length(aa)) return("-")
    names(sort(table(aa), decreasing = TRUE))[1]
  })
}

get_conserved_regions_aln <- function(msa_obj, min_length = 10){
  aln <- as.matrix(msa_obj)
  if(nrow(aln) < 2) return(data.frame(aln_start = 1L, aln_length = ncol(aln)))
  is_cons <- apply(aln, 2, function(col) all(col != "-") && length(unique(col)) == 1)
  if(!any(is_cons)) return(data.frame(aln_start = integer(0), aln_length = integer(0)))
  r <- rle(is_cons)
  starts <- cumsum(c(1, head(r$lengths, -1)))
  blocks <- data.frame(aln_start = starts[r$values], aln_length = r$lengths[r$values])
  dplyr::filter(blocks, aln_length >= min_length)
}

generate_peptides <- function(seq, len){
  n <- nchar(seq)
  if(n < len) return(character(0))
  vapply(1:(n - len + 1), function(i) substr(seq, i, i + len - 1), character(1))
}

clean_aa_string <- function(s){
  s <- toupper(as.character(s))
  s <- gsub("\\*", "", s)
  gsub("[^A-Z-]", "", s)
}

parse_iedb_tsv <- function(text){
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  lines <- lines[nzchar(lines)]
  start_idx <- grep("^allele\\t", lines)
  if(!length(start_idx)) stop("IEDB header not found")
  df <- read.table(text = paste(lines[start_idx:length(lines)], collapse = "\n"),
                   sep = "\t", header = TRUE, quote = "\"", comment.char = "", stringsAsFactors = FALSE)
  nm <- tolower(names(df)); names(df) <- nm
  ic50_col <- which(nm %in% c("ic50", "ic50(nm)"))
  out <- data.frame(Peptide = df$peptide, Allele = df$allele,
                    Rank = as.numeric(df$percentile_rank),
                    Affinity = as.numeric(df[[ic50_col[1]]]), stringsAsFactors = FALSE)
  out$AffinityClass <- cut(out$Affinity, c(-Inf, 50, 500, Inf),
                           labels = c("Strong", "Intermediate", "Weak"))
  out
}

query_iedb <- function(peps, allele, len, max_retries = 2L, timeout_sec = 180L){
  if(!length(peps)) return(NULL)
  for(i in seq_len(max_retries + 1L)){
    try({
      r <- POST(IEDB_API_URL,
                body = list(method = "netmhcpan", sequence_text = paste(peps, collapse = "\n"),
                            allele = allele, length = len, species = "human"),
                encode = "form", timeout(timeout_sec))
      if(r$status_code != 200) stop("IEDB HTTP ", r$status_code)
      return(parse_iedb_tsv(content(r, as = "text", encoding = "UTF-8")))
    }, silent = TRUE)
    if(i <= max_retries) Sys.sleep(2)
  }
  NULL
}

.load_cache <- local({env <- new.env(); env})
load_combined_db <- function(){
  if(!is.null(.load_cache$db)) return(.load_cache$db)
  if(!file.exists(combined_db_path)) return(NULL)
  db <- readRDS(combined_db_path)
  if(!inherits(db, "data.table")) db <- as.data.table(db)
  .load_cache$db <- db
  db
}

find_published_epitopes_local <- function(consensus_ungapped, uniprot_id){
  db <- load_combined_db(); if(is.null(db) || !nrow(db)) return(data.frame())
  f <- db[grepl(uniprot_id, MoleculeParentIRI, fixed = TRUE)]
  if(!nrow(f)) return(data.frame())
  neg_pat <- "negative|no\\s*binding|not\\s*detected|non.?binding|no\\s*activity"
  keep_t <- f$SourceFile == "tcell" & grepl("Positive", f$Assay, ignore.case = TRUE)
  keep_m <- f$SourceFile == "mhc" & !grepl(neg_pat, f$Assay, ignore.case = TRUE)
  f <- f[keep_t | keep_m]; if(!nrow(f)) return(data.frame())
  f[,`:=`(
    Start = vapply(Peptide, function(p){
      pos <- as.integer(regexpr(p, consensus_ungapped, fixed = TRUE))
      if(!is.na(pos) && pos > 0L) pos else NA_integer_
    }, integer(1)),
    End = vapply(Peptide, function(p){
      pos <- as.integer(regexpr(p, consensus_ungapped, fixed = TRUE))
      if(!is.na(pos) && pos > 0L) pos + nchar(p) - 1L else NA_integer_
    }, integer(1))
  )]
  m <- f[!is.na(Start)]
  if(nrow(m) > 0){
    m <- m[,.(Main_Epitope_ID = first(EpitopeID),
              EpitopeIDs = paste(unique(EpitopeID), collapse = ", "),
              Assays = paste(unique(Assay), collapse = " | "),
              `#_Assays` = .N,
              SourceFiles = paste(unique(SourceFile), collapse = "+"),
              Start = first(Start), End = first(End)), by = .(Peptide, Allele)]
  }
  m[]
}
