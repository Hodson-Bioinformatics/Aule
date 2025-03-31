#!/usr/bin/env Rscript
#
# Usage:
#   get_bams.R <input_table> <output_table> [<library_id_column>]
#
# Notes:
#   - <library_id_column> can be a name or the number 1. If it is the
#     number 1, the script assumes that the input file has no header
#     and only has a single column. The default value is 'library_id'
#     (i.e., the script assumes a file header).
#
#   - The column names are prefixed with 'al.', 'lc.', and 'lb.' for
#     designating alignment-, libcore-, and library-specific columns.
#     It also helps with ensuring no collisions when joining the
#     obtained information with the input table.
#


# Load packages -----------------------------------------------------------

suppressWarnings(suppressPackageStartupMessages({

message("Loading packages...")

required_packages <- c(
  # Included to avoid mysterious error
  "xml2",
  # Password prompt
  "getPass",
  # Data handling
  "readr",
  "dplyr",
  "purrr",
  "tidyr",
  "lubridate",
  # Web requests
  "httr",
  "jsonlite"
)

installed_packages <- rownames(installed.packages())

for (pkg in required_packages){
  if (!pkg %in% installed_packages) {
    message(paste0("Installing ", pkg, "..."))
    install.packages(pkg, repos = "https://mirror.rcg.sfu.ca/mirror/CRAN/")
  }
  library(pkg, character.only = TRUE)
}


# Determine arguments -----------------------------------------------------

# Parse command-line arguments (when running non-interactively)
args <- commandArgs(trailingOnly = TRUE) %>% as.list()
arg_names <- c("input_table", "output_table", "library_id_column")
args <- setNames(args, arg_names[1:length(args)])

# Set default value
if (is.null(args$library_id_column)) args$library_id_column <- "library_id"

# Check if table or single-column list
is_single_column <- args$library_id_column == "1"


# Load data ---------------------------------------------------------------

message("Loading data...")

if (is_single_column) {
  libraries <- read_tsv(args$input_table, col_names = "library_id",
                        col_types = cols())
  join_vector <- c("library_id" = "library_id")
} else {
  libraries <- read_tsv(args$input_table, col_types = cols())
  join_vector <- setNames("library_id", args$library_id_column)
  join_vector_2 <- setNames("lb.name", args$library_id_column)
}

message("Prompting for GSC/GIN credentials...")

credentials <- list()
credentials$username <- getPass("Username: ", noblank = TRUE)
credentials$password <- getPass("Password: ", noblank = TRUE)


# Setup session -----------------------------------------------------------

message("Setting up session...")

base_url <- handle("https://sbs.bcgsc.ca:8100/")

session_response <- POST(handle = base_url, path = "session",
                         encode = "json", body = credentials)
session_content  <- content(session_response)

if (status_code(session_response) != 200) {
  error <- pluck(session_content, "errors", 1, "description")

  if (is.null(error)) {
    print(session_content)
    stop(paste0("Failed to authenticate with the above error."))
  } else {
    stop(paste0("Failed to authenticate with the following error: \n", error))
  }
}

session_headers <- add_headers("X-Token" = session_content$token,
                               "Content-Type" = "application/json",
                               "Accept" = "application/json")

query <- function (endpoint, ..., as_df = NULL) {

  params <- list(...)

  response <- GET(handle = base_url, path = endpoint, query = params,
                  config = session_headers)

  contents <- content(response)

  if ( !is.null(as_df) ) {

    # Check if response is empty
    is_response_empty <- length(contents) == 0

    # If the response is empty, generate a hopefully non-empty response
    # with the given set of params that are expected to work
    if (is_response_empty) {
      response <- GET(handle = base_url, path = endpoint, query = as_df,
                      config = session_headers)
    }

    # Convert the response into a data frame
    contents <-
      response %>%
      content(as = "text", encoding = "UTF-8") %>%
      fromJSON() %>%
      jsonlite::flatten()

    # If the response was empty, remove the rows (but leave the columns)
    if (is_response_empty) {
      contents <- contents[FALSE,]
    }

  }

  contents

}


# Get folders containing BAM files using /merge endpoint ------------------

message("Retrieving list of merged BAM files...")

# Extract library IDs from input file
library_ids_vec <- libraries[[names(join_vector)]]
library_ids <- paste(library_ids_vec, collapse = ",")

# Check that library IDs are valid
existing_libs <- query("library", name = library_ids, as_df = list(name = "A47818"))
if ( length(library_ids_vec) > length(existing_libs$name) ) {
  missing_libs_vec <- setdiff(library_ids_vec, existing_libs$name)
  missing_libs     <- paste(missing_libs_vec, collapse = ", ")
  warning(paste0("The following library IDs do not exist: ", missing_libs))
  if ( length(library_ids_vec) == length(missing_libs_vec) ) {
    stop("After removing non-existing libraries, none are left to query.")
  }
}


# Retrieve BAM files for all libraries using /merge endpoint
merges_raw <-
  query("merge", library = library_ids, as_df = list(library = "A47818")) %>%
  # Drop multi-library alignments
  {
    multi_lib_alignments <- map_int(.$libraries, nrow) > 1
    if (sum(multi_lib_alignments) > 0) {
      warning("Dropping some multi-library alignments...")
    }
    .[!multi_lib_alignments,]
  } %>%
  transmute(library_id = map_chr(libraries, "name"),
            gsc_external_id = map_chr(libraries, "external_identifier"),
            gsc_patient_id = map_chr(libraries, "patient_identifier"),
            construction_date = map_chr(libraries, "library_started"),
            construction_date = as_date(construction_date, tz = "PT"),
            seq_strategy = map_chr(libraries, "library_strategy"),
            seq_strategy = case_when(seq_strategy == "miRNA-Seq" ~ "mirna",
                                     seq_strategy == "RNA-Seq" ~ "mrna",
                                     seq_strategy == "WGS" ~ "genome",
                                     seq_strategy == "SPC" ~ "capture",
                                     seq_strategy == "EXC-Seq" ~ "capture",
                                     TRUE ~ paste0("ERROR:", seq_strategy)),
            merge_status = status,
            num_merged = alc_count,
            analysis_type,
            removal_date = as_date(removed, tz = "PT"),
            reference_name = lims_genome_reference.symlink_name,
            reference_name_full = lims_genome_reference.name,
            aligner_name = aligner_software.path_string,
            aligner_name_full = aligner_software.name,
            sequence_length = sequence_length,
            object_id = map(merge_xrefs, "object_id"),
            object_type = map(merge_xrefs, "object_type"),
            bam_id = row_number(),
            data_path)

# Check if any library strategies are unknown
unknown_lib_strategies <-
  merges_raw %>%
  dplyr::filter(grepl("ERROR", seq_strategy))

if(nrow(unknown_lib_strategies) > 0) {
  unknown_lib_strategies_list <-
    unique(unknown_lib_strategies$seq_strategy) %>%
    sub("ERROR:", "", .) %>%
    paste(collapse = ", ")
  warning(paste0("Unknown library strategy encountered: ",
                 unknown_lib_strategies_list))
}

# Remove merges that failed, don't have a data path, or as duplicates
merges <-
  merges_raw %>%
  dplyr::filter(merge_status != "failed",
                !is.na(data_path))

# Check that each BAM file is uniquely described by the selected attributes
merges_attr <-
  merges %>%
  select(-starts_with("bam_"), -starts_with("object_"), -data_path) %>%
  filter(merge_status != "duplicate")
merges_attr_distinct <- distinct(merges_attr)
#stopifnot(nrow(merges_attr) == nrow(merges_attr_distinct))
if(nrow(merges_attr) != nrow(merges_attr_distinct)){
  message("Warning: Not all BAM files are uniquely described by the selected attributes")
}

# Check that all libraries have merges
missing_merges <- anti_join(libraries, merges, by = join_vector)

# stopifnot(nrow(missing_merges) == 0)
if (nrow(missing_merges) > 0) {
  missing_libs <- paste(missing_merges[[names(join_vector)]], collapse = ", ")
  message(paste0("These libraries do not have merges: ", missing_libs))
  message("Will attempt to retrieve non-merge data paths for these libraries. ")
}

# Append merge info to library info
bam_files <-
  merges %>%
  rename_at(vars(construction_date, seq_strategy,
                 gsc_external_id, gsc_patient_id),
            ~ paste0("lb.", .)) %>%
  rename_at(vars(-ends_with("_id"), -bam_id, -starts_with("object_"),
                 -contains(".")),
            ~ paste0("al.", .)) %>%
  inner_join(libraries, ., by = join_vector)


# Retrieve library preparation protocol -----------------------------------

message("Retrieving library preparation protocols...")

lib_protocols <-
  existing_libs$protocol_id %>%
  paste(collapse = ",") %>%
  query("protocol", id = ., as_df = list(id = "123")) %>%
  select(id, protocol_name = name) %>%
  right_join(select(existing_libs, name, protocol_id, ffpe, strandedness),
             by = c("id" = "protocol_id")) %>%
  select(-id) %>%
  rename_all(~ paste0("lb.", .))


# Get bio QC status and comments using /aligned_libcore endpoint ----------

message("Retrieving aligned_libcore IDs...")

# Get list of object IDs used to generate the above BAM files
# (mix of aligned_libcore and reposition IDs)
object_ids <-
  bam_files %>%
  unnest(c(object_type, object_id)) %>%
  select(bam_id, object_type, object_id)

# Get aligned_libcore IDs for reposition IDs
reposition_ids <-
  object_ids %>%
  dplyr::filter(object_type == "repo.analysis")

aln_libcore_ids_for_repositions <-
  reposition_ids$object_id %>%
  paste(collapse = ",") %>%
  query("reposition", id = ., as_df = list(id = "19971")) %>%
  transmute(object_id = id,
            object_type = "repo.analysis",
            aln_libcore_id = aligned_libcore.id)

# Check that the number of aligned_libcore IDs matches the number of
# reposition IDs
stopifnot(nrow(reposition_ids) == nrow(aln_libcore_ids_for_repositions))

# Replace reposition IDs with corresponding aligned_libcore IDs
aln_libcore_ids <-
  object_ids %>%
  left_join(aln_libcore_ids_for_repositions,
            by = c("object_id", "object_type")) %>%
  mutate(aln_libcore_id = ifelse(object_type == "metadata.aligned_libcore",
                                 object_id, aln_libcore_id)) %>%
  select(bam_id, object_type, aln_libcore_id)


message("Retrieving aligned_libcore Bio QC status and comments...")

# Obtain bio QC info using aligned_libcore IDs
bio_qc_info_raw <-
  aln_libcore_ids$aln_libcore_id %>%
  na.omit() %>%
  paste(collapse = ",") %>%
  query("aligned_libcore/info", id = ., as_df = list(id = "204689")) %>%
  select(aln_libcore_id = id, bio_qc_status = libcore.bio_qc_status,
         bio_qc_comments = libcore.bio_qc_comments) %>%
  mutate(bio_qc_comments = sub("^.*[(](N/A|\\w+) -> (N/A|\\w+)[)] ?;? ?;?",
                               "", bio_qc_comments))

# Check if number of retrieved aligned_libcores matches the number provided
missing_aln_libcores <-
  aln_libcore_ids %>%
  filter(object_type != "metadata.run") %>%
  anti_join(bio_qc_info_raw, by = "aln_libcore_id")
stopifnot(nrow(missing_aln_libcores) == 0)

# Check if any libcores have BioQC status set to anything other
# than Passed/Failed
known_statuses <- bio_qc_info_raw$bio_qc_status %in% c("Passed", "Failed")
stopifnot(all(known_statuses))

# Add bam_id to bio_qc_info for later joining
bio_qc_info <-
  left_join(aln_libcore_ids, bio_qc_info_raw, by = "aln_libcore_id") %>%
  replace_na(list(bio_qc_status = "Not_Available"))

# Count the number of libcores per QC status
bio_qc_status <-
  bio_qc_info %>%
  count(bam_id, bio_qc_status) %>%
  mutate(bio_qc_status = tolower(bio_qc_status),
         bio_qc_status = paste0("lc.num_", bio_qc_status)) %>%
  pivot_wider(names_from = bio_qc_status, values_from = n,
              values_fill = list(n = 0)) %>%
  select(bam_id, lc.num_passed, everything())

# Split comments by status and collapse using libcore IDs
bio_qc_comments <-
  bio_qc_info %>%
  dplyr::filter(!is.na(bio_qc_comments), bio_qc_comments != "") %>%
  separate_rows(bio_qc_comments, sep = "; ?(?=Failed|Warning|Manual)") %>%
  mutate(bio_qc_comments = ifelse(grepl(":", bio_qc_comments),
                                  bio_qc_comments,
                                  paste0("Other:", bio_qc_comments))) %>%
  separate(bio_qc_comments, c("status", "comment"),
           sep = ":", extra = "merge") %>%
  mutate(status = sub(" ", "_", status),
         status = tolower(status),
         status = paste0("lc.comments_", status),
         comment = sub('"', "", comment),
         comment = trimws(comment),
         comment = paste0(aln_libcore_id, "={", comment, "}")) %>%
  group_by(status) %>%
  mutate(row = row_number()) %>%
  pivot_wider(names_from = status, values_from = comment) %>%
  dplyr::select(-row) %>%
  group_by(bam_id) %>%
  summarise_at(vars(contains("comments_")),
               ~ paste(sort(na.omit(unique(.))), collapse = "|"))

# Add bio QC comments to bio QC status
bio_qc_combined <-
  left_join(bio_qc_status, bio_qc_comments, by = "bam_id") %>%
  mutate_at(vars(contains("comments_")), ~ ifelse(is.na(.), "", .))

# Combine library protocol and bio QC status and comments with merge info
bam_files_with_bioqc <-
  bam_files %>%
  left_join(bio_qc_combined, by = "bam_id") %>%
  left_join(lib_protocols, by = join_vector_2) %>%
  select(-bam_id, -starts_with("object_")) %>%
  select(one_of(colnames(libraries)), starts_with("al."), starts_with("lc."),
         starts_with("lb."), starts_with("sp."), starts_with("bp."),
         starts_with("pt."), everything())

# If no mising merges, output this table ------------------------------

if (nrow(missing_merges) == 0){
  write_tsv(bam_files_with_bioqc, args$output_table)
}

# Retrieve non-merge data paths if any---------------------------------

if (nrow(missing_merges) > 0) {
message("Retrieving data paths and bio QC for libraries without merges...")

no_merge_raw <-
  missing_merges[[args$library_id_column]] %>%
  paste(collapse = ",") %>%
  query("aligned_libcore/info", library = ., as_df = list(library = "A34795")) %>%
  mutate(bam_id = row_number())

no_merge_tidy <-
  no_merge_raw %>%
  filter(successful, !is.na(data_path), libcore.billable) %>%
  transmute(library_id = libcore.library.name,
            gsc_external_id = libcore.library.external_identifier,
            gsc_patient_id = libcore.library.patient_identifier,
            construction_date = libcore.library.library_started,
            construction_date = as_date(construction_date, tz = "PT"),
            seq_strategy = libcore.library.library_strategy,
            seq_strategy = case_when(seq_strategy == "miRNA-Seq" ~ "mirna",
                                     seq_strategy == "RNA-Seq" ~ "mrna",
                                     seq_strategy == "WGS" ~ "genome",
                                     seq_strategy == "SPC" ~ "capture",
                                     TRUE ~ paste0("ERROR:", seq_strategy)),
            merge_status = "no_merge",
            num_merged = NA,
            analysis_type = NA,
            removal_date = as_date(removed, tz = "PT"),
            reference_name = lims_genome_reference.symlink_name,
            reference_name_full = lims_genome_reference.name,
            aligner_name = analysis_software.path_string,
            aligner_name_full = analysis_software.name,
            sequence_length = sequence_length,
            bam_id,
            data_path) %>%
  rename_at(vars(construction_date, seq_strategy,
                 gsc_external_id, gsc_patient_id),
            ~ paste0("lb.", .)) %>%
  rename_at(vars(-ends_with("_id"), -bam_id, -contains(".")),
            ~ paste0("al.", .))

no_merge_bio_qc_comments <-
  no_merge_raw %>%
  transmute(bam_id,
            library_id = libcore.library.name,
            aln_libcore_id = id,
            bio_qc_comments = libcore.bio_qc_comments,
            bio_qc_comments = sub("^.*[(](N/A|\\w+) -> (N/A|\\w+)[)] ?;? ?;?",
                                  "", bio_qc_comments)) %>%
  separate_rows(bio_qc_comments, sep = "; ?(?=Failed|Warning|Manual)") %>%
  mutate(bio_qc_comments = ifelse(grepl(":", bio_qc_comments),
                                  bio_qc_comments,
                                  paste0("Other:", bio_qc_comments))) %>%
  separate(bio_qc_comments, c("status", "comment"),
           sep = ":", extra = "merge") %>%
  mutate(status = sub(" ", "_", status),
         status = tolower(status),
         status = paste0("lc.comments_", status),
         comment = sub('"', "", comment),
         comment = trimws(comment),
         comment = paste0(aln_libcore_id, "={", comment, "}")) %>%
  pivot_wider(names_from = status, values_from = comment) %>%
  group_by(bam_id) %>%
  summarise_at(vars(contains("comments_")),
               ~ paste(na.omit(unique(.)), collapse = "|"))

no_merge_bio_qc_status <-
  no_merge_raw %>%
  select(bam_id,
         bio_qc_status = libcore.bio_qc_status) %>%
  count(bam_id, bio_qc_status) %>%
  mutate(bio_qc_status = tolower(bio_qc_status),
         bio_qc_status = paste0("lc.num_", bio_qc_status)) %>%
  pivot_wider(names_from = bio_qc_status, values_from = n,
              values_fill = list(n = 0)) %>%
  select(bam_id, lc.num_passed, everything())

# Add bio QC comments to bio QC status
no_merge_bio_qc_combined <-
  left_join(no_merge_bio_qc_status, no_merge_bio_qc_comments, by = "bam_id") %>%
  mutate_at(vars(contains("comments_")), ~ ifelse(is.na(.), "", .))

# Combine library protocol and bio QC status and comments with merge info
no_merge_bams_with_bioqc <-
  missing_merges %>%
  left_join(no_merge_tidy, by = join_vector) %>%
  left_join(no_merge_bio_qc_combined, by = "bam_id") %>%
  left_join(lib_protocols, by = join_vector_2) %>%
  select(-bam_id) %>%
  select(one_of(colnames(libraries)), starts_with("al."), starts_with("lc."),
         starts_with("lb."), starts_with("sp."), starts_with("bp."),
         starts_with("pt."), everything())

# Combine merge ane no-merge tables --------------------------------------

all_bams_with_bioqc <-
  bind_rows(no_merge_bams_with_bioqc, bam_files_with_bioqc)



# Output table ------------------------------------------------------------

message("Outputting final table...")

write_tsv(all_bams_with_bioqc, args$output_table)

}

}))
