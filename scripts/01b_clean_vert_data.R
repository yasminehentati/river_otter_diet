########## clean vert diet data

# packages 

library(pacman)
p_load(tidyr,devtools,here,remotes,readr,dplyr,stringr,vroom,taxize)

##### set up data 

# read in OTUs 

# for the verts, we've already removed undetermined reads and merged similar sequences 

vert_seq <- read_csv(here("data_from_ellie", "12S_OTU_clean.csv"))
ncol(vert_seq)
colnames(vert_seq)

# pivot longer data and add a new column for "abundance"
merged_vars <- vert_seq %>%
  pivot_longer(
    cols = 10:936, 
    names_to = "sample_ID",                     
    values_to = "abundance"                     
  ) %>%
  dplyr::select("sample_ID", "Sequence", "ID",  "abundance") 

merged_vars$sample_ID
colnames(merged_vars)

# let's clean up the sample IDs 

# remove the additional strings from sample ID 
merged_vars$sample_ID <- gsub("_S.*", "", merged_vars$sample_ID) 

# remove the additional strings from ID 
merged_vars$ID <- sub("^[^_]+_([^_]+_[^_]+).*", "\\1", merged_vars$ID)

# rename the columns to be less confusing

colnames(merged_vars) <- c("sample_ID", "sequence", "processed_taxon", "abundance")

#### remove contamination 

# let's see what unique values we have 
unique_values <- unique(merged_vars$processed_taxon)

# just from a quick look at the data, LC070 seems to be a dog or coyote so we'll remove it 
merged_vars <- merged_vars[!grepl("LC070", merged_vars$sample_ID), ]

colnames(merged_vars)


# get a total sequence read number from each triplicate

filtered_data <- merged_vars %>%
  group_by(sample_ID) %>%
  mutate(total_reads = sum(abundance)) %>%
  ungroup()

filtered_data$total_reads


# for the next few steps we need to remove the extraction blanks and PCR negatives
# so they don't get manipulated 
# Filter out EB and PCRNEG rows before processing
real_samples <- filtered_data %>%
  filter(!str_detect(sample_ID, "^(EB|PCRNEG)"))

# get a total sequence read number from each scat 

# remove scats that have less than 100 total reads 
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
  group_by(base_sample, processed_taxon) %>%
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


# average all reads for all rows (samples, extraction blanks, PCR negatives)
avg_data <- filtered_data %>%
  mutate(base_sample = sub("_\\d+$", "", sample_ID)) %>%
  group_by(base_sample, processed_taxon) %>%
  summarise(
    abundance_avg = mean(abundance, na.rm = TRUE),
    total_reads_avg = mean(total_reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(sample_ID = base_sample)


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
  select(run_EB = sample_ID, processed_taxon, abundance_avg_EB = abundance_avg) %>%
  distinct()

run_PCRNEG_lookup <- clean_dat %>%
  filter(str_detect(sample_ID, "^PCRNEG")) %>%
  select(run_PCRNEG = sample_ID, processed_taxon, abundance_avg_PCRNEG = abundance_avg) %>%
  distinct() %>%
  mutate(run_PCRNEG = str_replace(run_PCRNEG, "PCRNEG-LC", "PCRNEG-LC_NGS_"))


# Step 2: Subtract EB and PCRNEG values from abundance_avg
clean_dat2 <- clean_dat %>%
  left_join(run_EB_lookup, by = c("run_EB", "processed_taxon")) %>%
  left_join(run_PCRNEG_lookup, by = c("run_PCRNEG", "processed_taxon")) %>%
  mutate(
    abundance_avg = if_else(
      str_detect(sample_ID, "^(EB|PCRNEG)"), 
      abundance_avg,  # Keep original for EB and PCRNEG rows
      abundance_avg - coalesce(abundance_avg_EB, 0) - coalesce(abundance_avg_PCRNEG, 0)  # Subtract for other samples, treat NA as 0
    )
  ) #%>%
#  select(-abundance_avg_EB, -abundance_avg_PCRNEG)  # Drop intermediate columns


# let's remove any rows that now have negative numbers -- this means they had too much contamination 
clean_dat2 <- clean_dat2 %>%
  filter(abundance_avg >= 0)

# final contamination clean -- we'll need to take out remaining reads of:
# coyote, bobcat, human, moutnain lion, otter, mink
# the last plate didn't use an otter blocking primer so we'll need to compare that to the others
# for now we'll keep minks

clean_dat2 <- clean_dat2[!grepl("Canis|Homo|Lontra|Lynx|Puma|Procyon|Vulpes", clean_dat2$processed_taxon), ]


# recalculate total_reads_avg
clean_dat2 <- clean_dat2 %>%
  group_by(sample_ID) %>%
  mutate(total_reads_avg = sum(abundance_avg, na.rm = TRUE)) %>%
  ungroup()

# let's save what we have so far
write_csv(clean_dat2, "cleaned_vert_data.csv")
clean_dat2 <- read_csv("cleaned_vert_data.csv")


# le'ts do some quick chekcs 

# Calculate means of abundance_avg by processed_taxon
summary_stats <- clean_dat2 %>%
  filter(str_detect(sample_ID, "LC")) %>%
  group_by(processed_taxon) %>%
  summarise(
    mean_abundance_avg = mean(abundance_avg, na.rm = TRUE)  # Mean of abundance_avg by processed_taxon
  )

# Calculate mean and standard deviation for total_reads_avg for the entire dataset
total_reads_stats <- clean_dat2 %>%
  filter(str_detect(sample_ID, "LC")) %>%
  summarise(
    mean_total_reads_avg = mean(total_reads_avg, na.rm = TRUE),  # Mean of total_reads_avg for the entire dataset
    sd_total_reads_avg = sd(total_reads_avg, na.rm = TRUE)       # Standard deviation of total_reads_avg for the entire dataset
  )



######## add taxonomic units 

unique(clean_dat2$processed_taxon)

# preprocess to handle " sp", " spp" and make clasification consistent
tax_dat <- clean_dat2 %>%
  mutate(
    processed_taxon = gsub("_sp\\b.*", "", processed_taxon),   # Remove " sp" and everything after it
    processed_taxon = gsub("_spp\\b.*", "", processed_taxon),  # Remove " spp" and everything after it
    processed_taxon = gsub("\\bMustela\\b", "Neogale", processed_taxon), #make mink consistent 
    processed_taxon = trimws(processed_taxon) 
  ) %>% 
  filter(!is.na(processed_taxon))



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


# create a data frame with the taxonomic hierarchy for missing values
nfrenata_taxonomy <- data.frame(
  name = c("cellular organisms",   "Eukaryota",            "Opisthokonta",         "Metazoa",              "Eumetazoa",           
           "Bilateria",            "Deuterostomia",        "Chordata" ,            "Craniata",             "Vertebrata",          
           "Gnathostomata",        "Teleostomi",           "Euteleostomi",         "Sarcopterygii",        "Dipnotetrapodomorpha",
           "Tetrapoda"  ,          "Amniota" ,             "Mammalia",             "Theria",               "Eutheria",           
           "Boreoeutheria",        "Laurasiatheria" ,      "Carnivora",            "Caniformia",           "Musteloidea",         
           "Mustelidae" ,          "Mustelinae" ,          "Neogale",              "Neogale frenata"),
  rank =  
    c("no rank",      "superkingdom", "clade"  ,      "kingdom" ,     "clade" ,       "clade" ,       "clade"   ,    
  "phylum",       "subphylum",    "clade"    ,    "clade"    ,    "clade"     ,   "clade"     ,   "superclass"  ,
  "clade" ,       "clade"   ,     "clade"   ,     "class"    ,    "clade"     ,   "clade"     ,   "clade"    ,   
  "superorder" ,  "order"  ,      "suborder" ,    "superfamily" , "family"    ,   "subfamily"  ,  "genus"  ,     
  "species"  ),
  stringsAsFactors = FALSE
)

taxonomy_results[["Neogale_frenata"]] <- nfrenata_taxonomy

# remove the original one that didn't work
taxonomy_results[["Neogale_frenata.NA"]] <- NULL

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
  select(-c(9))

# save our final data set!
write_csv(final_dat, here("clean_data", "12S_cleaned.csv"))
