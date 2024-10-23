#==============================================================================
# STANDARD AUXILIARY FUNCTIONS FOR DOWNSTREAM MGS WORKFLOW ANALYSIS
#==============================================================================

#------------------------------------------------------------------------------
# Functions for processing Nextflow logs
#------------------------------------------------------------------------------
# Assumes the presence of the following columns: process, cpus, realtime

parse_duration_string <- function(dstr){
  # Parse an atomic duration string
  if(dstr == "-") return(0)
  n <- sub("([0-9,.]+).*", "\\1", dstr) %>% as.numeric
  unit <- sub("[0-9,.]+(.*)", "\\1", dstr)
  if (unit == "ms") return(n / 1000)
  if (unit == "s") return(n)
  if (unit == "m") return(n * 60)
  if (unit == "h") return(n * 3600)
}

parse_duration_vector <- function(dvec){
  # Parse a vector of duration strings
  dvec %>% sapply(parse_duration_string) %>% sum
}

parse_durations <- function(durations){
  # Parse a vector of non-atomic duration strings
  durations %>% str_split(" ") %>% sapply(parse_duration_vector)
}

middle <- function(str, sep){
  str %>% str_split(sep) %>% sapply(head, -1) %>% sapply(tail, -1) %>% sapply(paste, collapse = sep)
}

import_log <- function(path){
  col_names <- c("process", "cpus", "realtime")
  log_raw <- read_tsv(path, col_names = col_names, show_col_types = FALSE)
  log <- log_raw %>% mutate(
    realtime_s = parse_durations(realtime),
    cpu_hours = realtime_s * cpus / 3600,
    workflow = process %>% str_split(":") %>% sapply(dplyr::first),
    job = process %>% str_split(":") %>% sapply(dplyr::last),
    subworkflow = process %>% middle(":")
  )
  return(log)
}

#------------------------------------------------------------------------------
# Functions for traversing taxid trees
#------------------------------------------------------------------------------

# Define auxiliary functions
get_child_taxids <- function(taxids, reference){
  # Get all direct child taxids for a list of parent taxids in a reference DB
  reference %>% filter(parent_taxid %in% taxids) %>% pull(taxid)
}
add_child_taxids <- function(taxids, reference){
  # For a list of taxids, find all direct child taxids in a reference DB and return a list of both parents and children
  get_child_taxids(taxids, reference) %>% append(taxids) %>% sort %>% unique
}

add_descendent_taxids <- function(taxids, reference){
  # For a list of taxids, find all descendent taxids in a reference DB and return a list of all taxids
  taxids_old <- taxids
  taxids_new <- add_child_taxids(taxids_old, reference)
  while (!identical(taxids_old, taxids_new)){
    taxids_old <- taxids_new
    taxids_new <- add_child_taxids(taxids_old, reference)
  }
  return(taxids_new)
}

raise_rank <- function(read_db, taxid_db, out_rank = "species", verbose = FALSE){
  # Given a read DB and a taxid DB, promote reads to a higher taxonomic rank
  # Assumes the following columns: seq_id, taxid, parent_taxid, rank, name
  # 1. Get higher ranks than search rank
  ranks <- c("subspecies", "species", "subgenus", "genus", "subfamily", "family", "suborder", "order", "class", "subphylum", "phylum", "kingdom", "superkingdom")
  rank_match <- which.max(ranks == out_rank)
  high_ranks <- ranks[rank_match:length(ranks)]
  # 2. Merge read DB and taxid DB
  reads <- read_db %>% select(-any_of(c("parent_taxid", "rank", "name"))) %>%
    left_join(taxid_db, by="taxid")
  # 3. Extract sequences that are already at appropriate rank
  reads_rank <- filter(reads, rank == out_rank)
  # 4. Drop sequences at a higher rank and return unclassified sequences
  reads_norank <- reads %>% filter(rank != out_rank, !rank %in% high_ranks, !is.na(taxid))
  while(nrow(reads_norank) > 0){ # As long as there are unclassified sequences...
    # Promote read taxids and re-merge with taxid DB, then re-classify and filter
    reads_remaining <- reads_norank %>% mutate(taxid = parent_taxid) %>%
      select(-parent_taxid, -rank, -name) %>%
      left_join(taxid_db, by="taxid")
    reads_rank <- reads_remaining %>% filter(rank == out_rank) %>%
      bind_rows(reads_rank)
    reads_norank <- reads_remaining %>%
      filter(rank != out_rank, !rank %in% high_ranks, !is.na(taxid))
  }
  # 5. Finally, extract and append reads that were excluded during the process
  reads_dropped <- reads %>% filter(!seq_id %in% reads_rank$seq_id)
  reads_out <- reads_rank %>% bind_rows(reads_dropped) %>%
    select(-parent_taxid, -rank, -name) %>%
    left_join(taxid_db, by="taxid")
  return(reads_out)
}