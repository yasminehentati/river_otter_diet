########## clean invert diet data

# packages
library(pacman)
p_load(tidyr,devtools,here,remotes,readr,dplyr,stringr,purrr)


# download taxize 
# since it's orphaned on cran, we'll need to use the github version
# make sure you have a github PAT 
#remotes::install_github("ropensci/bold")
#remotes::install_github("ropensci/taxize")
View(allL_diet)
library(taxize)

# get NCBI key for taxize 
#taxize::use_entrez()
#usethis::edit_r_environ()

#read in OTUs 
### set up data 

invert_seq <- read_csv(here("data_from_ellie", "16S_OTU.csv"))
colnames(invert_seq)
# merged variants with some taxo data 
merged_vars <- read_csv(here("data_from_ellie", "16S_vars_summed.csv")) 

# pivot longer data and add a new column for "abundance"
invert_seq <- invert_seq %>%
  pivot_longer(
    cols = 84:1038, 
    names_to = "sample_ID",                     
    values_to = "abundance"                     
  ) %>%
  dplyr::select("sample_ID", "Alignment name", "Original Sequence", "Modified sequence", "abundance") 

ncol(merged_vars)

# pivot longer summed variants sheet 
merged_vars <- merged_vars %>%
  pivot_longer(
    cols = 91:713, 
    names_to = "sample_ID",                     
    values_to = "abundance"                     
  ) %>%
  dplyr::select("sample_ID", "Alignment name", "Possible taxon", "Original Sequence", "Modified sequence",  "abundance") # Keep relevant columns

# let's clean up the sample IDs 

merged_vars$sample_ID <- gsub("\\.\\.\\..*", "", merged_vars$sample_ID)

### merge similar sequence variants 

## we'll use the alignment name column to join data with merged variants sheet from ellie 

# first add the possible taxon column to the data 
merged_data <- invert_seq %>%
  left_join(merged_vars %>% dplyr::select("Alignment name", "Possible taxon"), 
            by = c("Alignment name" = "Alignment name"))

# change col names to be more R friendly
colnames(merged_data) <- c("sample_ID", "alignment_name", "original_seq", "modified_seq", 
                           "abundance", "possible_taxon")

# # remove mammal sequences 
merged_data <- merged_data %>% 
  filter(!grepl('Mammal', alignment_name))

# remove fish sequences (cottidae mainly)
merged_data <- merged_data %>% 
  filter(!grepl('Teleostei', alignment_name)) %>% 
  filter(!grepl('Possible_fish', alignment_name)) %>%
  filter(!grepl('Aves', alignment_name)) %>% 
  filter(!grepl('otter', alignment_name)) %>% 
  filter(!grepl('rodent', alignment_name)) %>%
  filter(!grepl('plant', alignment_name)) 

# aggregate data by "possible taxon" (assigned by ellie to be the same species or taxo group)
merged_seqs <- merged_data %>%
  group_by(sample_ID, possible_taxon) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

merged_seqs <- merged_data %>%
  group_by(sample_ID, possible_taxon) %>%
  summarise(
    abundance = sum(abundance, na.rm = TRUE), 
    across(everything(), ~ if(is.numeric(.)) sum(., na.rm = TRUE) else first(.)), 
    .groups = "drop"
  )



### get final sequences

# get a total sequence read number from each triplicate

filtered_data <- merged_seqs %>%
  group_by(sample_ID) %>%
  mutate(total_reads = sum(abundance)) %>%
  ungroup()
filtered_data$total_reads

# for the next few steps we need to remove the extraction blanks and PCR negatives
# so they don't get manipulated 
# Filter out EB and PCRNEG rows before processing
real_samples <- filtered_data %>%
  filter(!str_detect(sample_ID, "^(EB|PCRNEG)"))

# remove triplicates that have less than 100 total reads across all 3
real_samples <- real_samples %>%
  mutate(base_sample = sub("_\\d+$", "", sample_ID)) %>%
  group_by(base_sample) %>%
  mutate(total_sample_reads = sum(total_reads)) %>%
  filter(total_sample_reads >= 100) %>%
  ungroup() %>%
  select(-total_sample_reads)


# remove scats where only 1/3 triplicate have reads 
real_samples <- real_samples %>%
  mutate(base_sample = sub("_\\d+$", "", sample_ID)) %>%
  group_by(base_sample, possible_taxon) %>%
  mutate(non_zero_count = sum(abundance > 0)) %>%
  filter(non_zero_count >= 2) %>%
  ungroup() %>%
  select(-non_zero_count)

nrow(real_samples)

# Re-add EB and PCRNEG rows to the dataset (no manipulation)
filtered_data <- bind_rows(
  real_samples,
  filtered_data %>% filter(str_detect(sample_ID, "^(EB|PCRNEG)"))
)


#replace all occurrences of 118 with "Bird mite"
filtered_data <- filtered_data %>%
  mutate(possible_taxon = if_else(possible_taxon == 118, "Bird mite", possible_taxon))

# average all reads for all rows (samples, extraction blanks, PCR negatives)
avg_data <- filtered_data %>%
  mutate(base_sample = sub("_\\d+$", "", sample_ID)) %>%
  group_by(base_sample, possible_taxon) %>%
  summarise(
    abundance_avg = mean(abundance, na.rm = TRUE),
    total_reads_avg = mean(total_reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(sample_ID = base_sample)



#### clean up contamination and final cleaning 

## assign extraction blanks and PCR neg to each sample based on genetics labwork spreadsheet

#read in labwork spreadsheet  
ngs_info <- read_csv(here("data_from_ellie", "NGS_plate_data.csv"))

ngs_info <- ngs_info %>%
  filter(!is.na(ScatID)) %>% 
  mutate(`NGS-Plate-Name` = paste0("PCRNEG-", `NGS-Plate-Name`)) %>% 
  rename(sample_ID = ScatID, run_PCRNEG = "NGS-Plate-Name", run_EB = "Extraction blank")

## remove extraction blank average and PCR neg average from sample average

# first, EB 
clean_dat <- avg_data %>%
  left_join(
    ngs_info %>%
      select(sample_ID, run_EB),
    by = "sample_ID"
  )

# now PCR negs -- create column and join 
clean_dat <- clean_dat %>%
  left_join(
    ngs_info %>%
      select(sample_ID, run_PCRNEG),  
    by = "sample_ID"
  )


## now subtract the matching average EB and PCRNEG reads for each species/sample combo 
# Step 1: Create EB and PCRNEG lookup tables by taxon
run_EB_lookup <- clean_dat %>%
  filter(str_detect(sample_ID, "^EB")) %>%
  select(run_EB = sample_ID, possible_taxon, abundance_avg_EB = abundance_avg) %>%
  distinct()

run_PCRNEG_lookup <- clean_dat %>%
  filter(str_detect(sample_ID, "^PCRNEG")) %>%
  select(run_PCRNEG = sample_ID, possible_taxon, abundance_avg_PCRNEG = abundance_avg) %>%
  distinct() %>%
  mutate(run_PCRNEG = str_replace(run_PCRNEG, "PCRNEG-LC", "PCRNEG-LC_NGS_"))

# Step 2: Subtract EB and PCRNEG values from abundance_avg
clean_dat2 <- clean_dat %>%
  left_join(run_EB_lookup, by = c("run_EB", "possible_taxon")) %>%
  left_join(run_PCRNEG_lookup, by = c("run_PCRNEG", "possible_taxon")) %>%
  mutate(
    abundance_avg = if_else(
      str_detect(sample_ID, "^(EB|PCRNEG)"), 
      abundance_avg,  # Keep original for EB and PCRNEG rows
      abundance_avg - abundance_avg_EB - abundance_avg_PCRNEG  # Subtract for other samples
    )
  ) #%>%
#  select(-abundance_avg_EB, -abundance_avg_PCRNEG)  # Drop intermediate columns


# let's remove any rows that now have negative numbers -- this means they had too much contamination 
clean_dat2 <- clean_dat2 %>%
  filter(abundance_avg >= 0)

# recalculate total_reads_avg
clean_dat2 <- clean_dat2 %>%
  group_by(sample_ID) %>%
  mutate(total_reads_avg = sum(abundance_avg, na.rm = TRUE)) %>%
  ungroup()

# let's save what we have so far
write_csv(clean_dat2, "cleaned_invert_data.csv")
clean_dat2 <- read_csv("cleaned_invert_data.csv")

# remove scientific notation and round decimals
options(scipen = 999)
#clean_dat$abundance_avg <- round(clean_dat$abundance_avg, 4)

# Preprocess possible_taxon to handle " sp", " spp", "Possible ", and numbers
tax_dat <- clean_dat2 %>%
  mutate(
    processed_taxon = gsub("^Possible\\s+", "", possible_taxon),  # Remove "Possible " at the beginning
    processed_taxon = gsub("\\s+sp\\b.*", "", processed_taxon),   # Remove " sp" and everything after it
    processed_taxon = gsub("\\s+spp\\b.*", "", processed_taxon),  # Remove " spp" and everything after it
    processed_taxon = gsub("\\d.*", "", processed_taxon), # Remove numbers and everything after them
    processed_taxon = gsub("\\s+or\\s+.*", "", processed_taxon),  # Keep only the first word before " or "
    processed_taxon = trimws(processed_taxon) 
  ) %>% 
  filter(!is.na(processed_taxon))

tax_dat$processed_taxon


# now let's clean up any common names to the corresponding taxonomic name 
# since we have some vague common names we'll have to do it manually

taxon_mapping <- c(
  "copepod" = "Copepoda",
  "fly" = "Diptera",
  "insect" = "Insecta",
  "midge" = "Chironomidae",
  "shrimp" = "Decapoda",
  "springtail" = "Collembola",
  "water flea" = "Daphnia",
  "true bug" = "Hemiptera",
  "amphipod" = "Amphipoda",
  "Amphipod" = "Amphipoda",
  "caddisfly" = "Trichoptera",
  "beetle" = "Coleoptera",
  "booklouse" = "Psocoptera",
  "bird mite" = "Dermanyssus",
  "Bird mite" = "Dermanyssus",
  "bacteria" = "Bacteria",
  "freshwater snail" = "Gastropoda",
  "Unknown arthropod" = "Arthropoda",
  "decapod" = "Decapoda",
  "cricket" = "Gryllidae",
  "Rove beetle" = "Staphylinidae",
  "greenhouse millipede" = "Oxidus gracilis",
  "mite" = "Acari",
  "hydrophilid" = "Hydrophilidae",
  "gastropod" = "Gastropoda",
  "aphid" = "Aphidoidea",
  "isopod" = "Isopoda",
  "Brachypanorpa oregensis" = "Brachypanorpa oregonensis",
  "Culiseta longiareolata OQ" = "Culiseta longiareolata",
  "lapathi" = "Cryptorhynchus lapathi",
  "Lapathi" = "Cryptorhynchus lapathi",
  "Licerus" = "Lirceus"
  
)

# Replace common names with corresponding Latin names 
tax_dat <- tax_dat %>%
  mutate(processed_taxon = recode(processed_taxon, !!!taxon_mapping))


# now we need to add and fill in our taxonomic columns 
# based on the data we already ahve 
# this code will fetchup upstream taxonomy using NCBI based on the most specific 
# data you have for each row. you'll need an API key 
# the code will retry a few times if it runs into API/server issues,
# pause for 5 seconds if it needs to, then time out after 3 tries


# function to fetch taxonomy using taxize with our retry mechanism
fill_taxonomy_with_retry <- function(taxon, max_retries = 3, pause_time = 5) {
  retries <- 0
  tax_data <- NULL
  
  while (retries < max_retries && is.null(tax_data)) {
    tax_data <- tryCatch({
      taxon_info <- get_uid(taxon)  # get NCBI UID first
      if (length(taxon_info) > 0) {
        # retrieve full taxonomy based on the NCBI UID
        taxonomy <- classification(taxon_info)
        return(taxonomy)
      } else {
        return(NULL)
      }
    }, error = function(e) {
      message("Error fetching taxonomy for ", taxon, ": ", e$message)
      return(NULL)
    })
    
    # if tax_data is NULL, increment retry counter and pause
    if (is.null(tax_data)) {
      retries <- retries + 1
      message("Retry ", retries, " for taxon: ", taxon)
      Sys.sleep(pause_time) 
    }
  }
  
  return(tax_data)
}


# get unique processed_taxon values
unique_taxa <- unique(tax_dat$processed_taxon)

# fetch taxonomy data for each unique taxon
taxonomy_results <- sapply(unique_taxa, fill_taxonomy_with_retry)

## we had 1 organism that didn't populate, so we'll fix that now
# create a data frame with the taxonomic hierarchy for missing values
euphilomedes_taxonomy <- data.frame(
  name = c("Animalia", "Arthropoda", "Crustacea", "Oligostraca", "Ostracoda", "Myodocopa", "Myodocopida"),
  rank = c("kingdom", "phylum", "subphylum", "superclass", "class", "subclass", "order"),
  stringsAsFactors = FALSE
)

taxonomy_results[["Euphilomedes producta"]] <- euphilomedes_taxonomy

# remove the original one that didn't work
taxonomy_results[["Euphilomedes producta.NA"]] <- NULL

# clean the names of the data frames in taxonomy_results
names(taxonomy_results) <- gsub("\\.[0-9]+$", "", names(taxonomy_results))

# identify all unique taxonomic ranks

unique_ranks <- unique(unlist(lapply(taxonomy_results, function(df) {
  if (!is.data.frame(df)) {
    df <- as.data.frame(df)
  }
  if (!all(c("rank", "name") %in% colnames(df))) {
    stop("Each data frame must contain 'rank' and 'name' columns.")
  }
  return(df$rank)
})))



# make an empty data frame to fill our ranks 
rank_df <- data.frame(organism = names(taxonomy_results), stringsAsFactors = FALSE)
for (rank in unique_ranks) {
  rank_df[[rank]] <- NA
}

# fill data frame with our data from the taxonomy results function
for (i in seq_along(taxonomy_results)) {
  organism <- names(taxonomy_results)[i]
  taxonomy <- taxonomy_results[[i]]
  for (j in seq_len(nrow(taxonomy))) {
    rank <- taxonomy$rank[j]
    name <- taxonomy$name[j]
    rank_df[rank_df$organism == organism, rank] <- name
  }
}

# combine the original data with the new taxonomy columns
final_dat <- left_join(tax_dat, rank_df, by = c("processed_taxon" = "organism")) %>% 
  select(-c(10))

colnames(final_dat)

# save our final data set!
write_csv(final_dat, here("clean_data", "16S_cleaned.csv"))

