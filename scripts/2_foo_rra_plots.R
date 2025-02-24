library(pacman)
p_load(here, dplyr, ggplot2, treemap, tidyr, vegan, MicroNiche, RInSp, iNEXT, ggforce, goeveg, strex, glmmTMB,
       RColorBrewer)

#################


# read in data
# verts first  
# let's keep only our actual otter samples -- delete controls 

vert_dat <- read_csv(here("clean_data", "12S_cleaned.csv")) %>%
  filter(!grepl("PCR|EB|Otter", sample_ID))


# let's also remove 6 scats we suspect to be raccoon scats based on manual analysis
remove_samples <- c("LC288A", "LC225", "LC208", "LC19", "LC197", "LC005")

vert_dat <- vert_dat %>%
  filter(!sample_ID %in% remove_samples)

unique(vert_dat$sample_ID)

# let's cut all rows that are less than 0.05% of total average reads 
vert_dat <- vert_dat %>%
  mutate(threshold = 0.0005 * total_reads_avg) %>%  # Calculate 0.05% of total_reads_avg
  filter(abundance_avg >= threshold) %>%  # Keep rows where abundance_avg is greater than or equal to the threshold
  select(-threshold)  # Remove the threshold column after filtering


# fix missing taxonomy 

vert_dat <- vert_dat %>%
  mutate(order = ifelse(species == "Cymatogaster aggregata", "Blenniiformes", order))

vert_dat <- vert_dat %>%
  mutate(order = ifelse(family == "Cottidae", "Perciformes", order))

vert_dat <- vert_dat %>%
  mutate(order = ifelse(family == "Psychrolutidae", "Perciformes", order))



##### 

# now let's read in  inverts 

# let's keep only our actual otter samples -- delete controls 

invert_dat <- read_csv(here("clean_data", "16S_cleaned.csv"))%>%
  filter(!grepl("PCR|EB|Otter", sample_ID))

# let's also remove 6 scats we suspect to be raccoon scats based on manual analysis
remove_samples <- c("LC288A", "LC225", "LC208", "LC19", "LC197", "LC005")

invert_dat <- invert_dat %>%
  filter(!sample_ID %in% remove_samples)

# also delete bacteria
invert_dat<- invert_dat[!grepl("Bacteria", invert_dat$superkingdom), ] 


# let's cut all rows that are less than 1% of total average reads 
invert_dat <- invert_dat %>%
  mutate(threshold = 0.01 * total_reads_avg) %>%  # Calculate 1% of total_reads_avg
  filter(abundance_avg >= threshold) %>%  # Keep rows where abundance_avg is greater than or equal to the threshold
  select(-threshold)  # Remove the threshold column after filtering

# fix one family that is incorrect 
invert_dat  <-invert_dat  %>%
  mutate(family = ifelse(family == "Atelecyclidae", "Cheiragonidae", family))



# rewrite the files 

# vert_dat <- write_csv(vert_dat, here("clean_data", "12S_final.csv"))

# invert_dat <- write_csv(invert_dat, here("clean_data", "16S_final.csv"))

#######################
# start from here for analyses 
#####  RRA plots --- vertebrates 

# Calculate the total summed abundance_avg for all families (total abundance)
total_abundance_sum <- sum(vert_dat$abundance_avg, na.rm = TRUE)

# Calculate the total abundance_avg for each family
family_abundance_sum <- vert_dat %>%
  group_by(family,class) %>%
  summarise(total_abundance_avg = sum(abundance_avg, na.rm = TRUE)) %>%
  ungroup()

# Calculate RRA by dividing the total abundance_avg for each family by the total summed abundance_avg
family_abundance_sum <- family_abundance_sum %>%
  mutate(RRA = total_abundance_avg / total_abundance_sum) #%>% 
 #filter(RRA >= 0.0005)  # Filter families with RRA >= 0.05%

# Check that the sum of RRA equals almost 1 (minus the ones we took out)
sum(family_abundance_sum$RRA, na.rm = TRUE)

family_abundance_sum <- family_abundance_sum[order(family_abundance_sum$RRA, decreasing=T),]
family_abundance_sum$family <- factor(family_abundance_sum$family, 
                                      levels=rev(as.character(family_abundance_sum$family)))


# Define custom labels for the legend
custom_labels_vert <- c(
  "Actinopteri" = "Ray-finned fishes",
  "Mammalia" = "Mammals",
  "Aves" = "Birds",
  "Amphibia" = "Amphibians",
  "Hyperoartia" = "Lampreys"
  # Add more custom labels as needed
)



# define custom colors so they remain the same regardless of which classes show up
classes_vert <- c("Aves", "Mammalia", "Actinopteri", "Hyperoartia", "Amphibia")

# Get the first 5 colors from the Set2 palette
set2_colors_vert <- brewer.pal(5, "Set2")

# Create a named vector
class_colors_vert <- setNames(set2_colors_vert, classes_vert)


# Create the RRA bar plot
ggplot(family_abundance_sum, aes(x = family, y = RRA, fill = class)) + 
  geom_bar(stat = "identity", color = 'black') +
  geom_text(aes(label = sprintf("%.3f%%", RRA * 100)), hjust = -0.1, size = 3.5) + 
  coord_flip() +
  labs(x = "Family", y = "Relative Read Abundance (RRA)", fill = "Class") +
  scale_fill_manual(values = class_colors_vert,labels=custom_labels_vert) +  # Apply the custom colors
  theme_classic() +
  theme(axis.text = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 13))



#####  RRA plots --- invertebrates 



# Calculate the total summed abundance_avg for all families (total abundance)
total_abundance_sum2 <- sum(invert_dat$abundance_avg, na.rm = TRUE)

# Calculate the total abundance_avg for each family
family_abundance_sum2 <- invert_dat %>%
  group_by(kingdom,phylum,class,order,family) %>%
  summarise(total_abundance_avg = sum(abundance_avg, na.rm = TRUE)) %>%
  ungroup() 


# Calculate RRA by dividing the total abundance_avg for each family by the total summed abundance_avg
invert_RRA <- family_abundance_sum2 %>%
  mutate(RRA = total_abundance_avg / total_abundance_sum2) %>% 
  filter(RRA >= 0.001) %>% # Filter for RRA >= 1%
  drop_na(family) %>% # remove unidentified fam
  filter(class != "Tremellomycetes")


# let's fill in NAs 

invert_RRA <- invert_RRA %>%   
  mutate(class = case_when(
    phylum == "Arthropoda" & is.na(class) ~ "Arthropoda (undet. class)",
    TRUE ~ class
  )) %>% 
  mutate(order = case_when(
    class == "Arthropoda (undet. class)" & is.na(order) ~ "Arthropods (undet.)",
    class == "Insecta" & is.na(order) ~ "Insects (undet.)",
    class == "Arachnida" & is.na(order) ~ "Arachnids (undet.)",
    TRUE ~ order
  )) %>%   mutate(family = case_when(
    order == "Amphipoda" & is.na(family) ~ "Amphipoda (undet.)",
    order == "Arthropods (undet.)" & is.na(family) ~ "Arthropoda (undet.)",
    order == "Insecta" & is.na(family) ~ "Insecta (undet.)",
    order == "Arachnids (undet.)" & is.na(family) ~ "Arachnida (undet.)",
    class == "Copepoda" & is.na(family) ~ "Copepods (unk. family)"
    TRUE ~ family
  ))
 

# Check that the sum of RRA equals almost 1 (minus the ones we took out)
sum(invert_RRA$RRA, na.rm = TRUE)

invert_RRA <- invert_RRA[order(invert_RRA$RRA, decreasing=T),]
invert_RRA$family <- factor(invert_RRA$family, 
                                       levels=rev(as.character(invert_RRA$family)))

unique(invert_RRA$order)
unique(invert_RRA$family)
unique(invert_RRA$class)
unique(invert_RRA$phylum)


# we'll need to reset our colors because we have more representation now 
# Define custom labels for the legend -- class

custom_labels <- c(
  "Malacostraca" = "Malacostracans",
  "Insecta" = "Insects",
  "Branchiopoda" = "Branchiopods",
  #  "Arthropoda" = "Arthropods",
  # "Gastropoda" = "Gastropods",
  "Cestoda" = "Cestodes",
  "Arachnida" = "Arachnids",
  "Collembola" = "Springtails",
  "Bivalvia" = "Bivalves",
  #  "Hexanauplia" = "Copepods",
  "Diplopoda" = "Millipedes"
)

# define custom colors so they remain the same regardless of which classes show up
classes <-c(
  "Malacostraca",
  "Insecta",
  "Branchiopoda",
  #  "Arthropoda",
  # "Gastropoda",
  "Cestoda",
  "Arachnida",
  "Collembola",
  "Bivalvia",
  # "Copepoda",
  "Diplopoda"
)

# Define custom labels for the legend -- class
custom_labels <- c(
  "Malacostraca" = "Malacostracans",
  "Insecta" = "Insects",
  "Branchiopoda" = "Branchiopods",
  "Arthropoda (undet. class)" = "Arthropods (unk. class)",
  "Gastropoda" = "Gastropods (undet. family)",
  "Cestoda" = "Cestodes",
  "Arachnida" = "Arachnids",
#  "Collembola" = "Springtails",
  "Bivalvia" = "Bivalves"
)


# define custom colors so they remain the same regardless of which classes show up
classes <-c(
  "Malacostraca",
  "Insecta",
  "Branchiopoda",
  "Arthropoda (undet. class)",
  "Gastropoda",
  "Cestoda",
  "Arachnida",
#  "Collembola",
  "Bivalvia"
)


# Get the first colors from the Set1 palette
set1_colors <- brewer.pal(9, "Set1")

# Create a named vector
class_colors <- setNames(set1_colors, classes)

# Create the RRA bar plot
ggplot(invert_RRA, aes(x = family, y = RRA, fill = class)) + 
  geom_bar(stat = "identity", color = 'black') +
  geom_text(aes(label = sprintf("%.2f%%", RRA * 100)), hjust = -0.1, size = 3.5) + 
  coord_flip() +
  labs(x = "Family", y = "Relative Read Abundance (RRA)", fill = "Class") +
  scale_fill_manual(values = class_colors,labels=custom_labels) +  # Apply the custom colors
  theme_classic() +
  theme(axis.text = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        legend.text = element_text(size = 13))

########################
# let's calculate frequency of occurrence (FOO) 


### Vert FOO -- by family
# we need to remove families that occur more than once in each scat 

vert_dat <- vert_dat %>%
  distinct(sample_ID, family, .keep_all = TRUE)


# Calculate Frequency of Occurrence (FOO)
FOO_All_Samples <- vert_dat %>%
  dplyr::count(family) %>%
  #  filter(family != "") %>%
  #  drop_na() %>%
  rename(Type = family, N = n) %>%
  arrange(desc(N))


# Join with vert_dat to add 'class' and 'order' information
FOO_All_Samples <- FOO_All_Samples %>%
  left_join(vert_dat %>% select(family, class, order) %>% 
              distinct(family, .keep_all = TRUE), by = c("Type" = "family"))


# Order families by frequency
FOO_All_Samples <- FOO_All_Samples %>%
  arrange(desc(N))

# Convert 'Type' to a factor with levels ordered by frequency
FOO_All_Samples$Type <- factor(FOO_All_Samples$Type, levels = rev(FOO_All_Samples$Type))



# Create the bar plot
ggplot(FOO_All_Samples, aes(x = Type, y = N, fill = class)) + 
  geom_bar(stat = "identity", color = 'black') +
  #geom_text(aes(label = sprintf("%.1f%%", (N / sum(N)) * 100)), hjust = -0.1, size = 3.5) + 
  coord_flip() +
  labs(x = "Family", y = "Frequency of Occurrence (Number of Samples)", fill = "Class") + 
  scale_fill_manual(values = class_colors_vert,labels=custom_labels_vert) +  # Apply the custom colors
  theme_classic() +
  theme(axis.text = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        legend.text = element_text(size = 13))


######## 
# invert FOO
# Calculate Frequency of Occurrence (FOO)

# we need to remove families that occur more than once in each scat 
invert_dat <- invert_dat %>%
  distinct(sample_ID, family, .keep_all = TRUE) %>%
  filter(class != "Tremellomycetes")



# get rid of NAs 
#invert_dat2 <- invert_dat %>%   
 # mutate(class = case_when(
 #   phylum == "Arthropoda" & is.na(class) ~ "Arthropoda (undet. class)",
 #   TRUE ~ class
 # )) %>% 
#  mutate(order = case_when(
#    class == "Arthropoda (undet. class)" & is.na(order) ~ "Arthropods (undet.)",
#    class == "Insecta" & is.na(order) ~ "Insects (undet.)",
  #  class == "Arachnida" & is.na(order) ~ "Arachnids (undet.)",
#    TRUE ~ order
 #)) %>% 
#  mutate(family = case_when(
#    order == "Amphipoda" & is.na(family) ~ "Amphipoda (unk. fam)",
#    class == "Arthropoda" & is.na(family) ~ "Arthropoda (unk. fam)",
#    class == "Insecta" & is.na(family) ~ "Insecta (unk. fam)",
#    class == "Arachnida" & is.na(family) ~ "Arachnida (unk. fam)",
#    class == "Hexanauplia" & is.na(family) ~ "Copepoda (unk. fam)",
#    class == "Gastropoda" & is.na(family) ~ "Gastropoda (unk. fam)",
#    TRUE ~ family
#  ))


# Calculate Frequency of Occurrence (FOO)
FOO_All_Samples <- invert_dat %>%
  dplyr::count(family) %>%
  drop_na(family) %>%
  rename(Type = family, N = n) %>%
  arrange(desc(N))
# View(FOO_All_Samples)

FOO_All_Samples <- FOO_All_Samples %>%
  mutate(FOO = N / total_samples)

#FOO_All_Samples <- FOO_All_Samples # %>%
  #filter(FOO >= 0.009)

# Join with invert_dat to add 'class' and 'order' information
FOO_All_Samples <- FOO_All_Samples %>%
  left_join(invert_dat %>% select(family, class, order) %>% 
              distinct(family, .keep_all = TRUE), by = c("Type" = "family"))

# Fill in NAs
#FOO_All_Samples <- FOO_All_Samples %>%   
#  mutate(Type = case_when(
#    class == "Gastropoda" & is.na(Type) ~ "Gastropoda (undet. family)",
#    TRUE ~ Type
#  )) 

# Convert 'Type' to a factor with levels ordered by frequency
FOO_All_Samples$Type <- factor(FOO_All_Samples$Type, levels = rev(FOO_All_Samples$Type))
unique(FOO_All_Samples$class)

# we'll need to reset our colors because we have more representation now 
# Define custom labels for the legend -- class

custom_labels <- c(
  "Malacostraca" = "Malacostracans",
  "Insecta" = "Insects",
  "Branchiopoda" = "Branchiopods",
 #  "Arthropoda" = "Arthropods",
 # "Gastropoda" = "Gastropods",
  "Cestoda" = "Cestodes",
  "Arachnida" = "Arachnids",
  "Collembola" = "Springtails",
  "Bivalvia" = "Bivalves",
#  "Hexanauplia" = "Copepods",
  "Diplopoda" = "Millipedes"
)

# define custom colors so they remain the same regardless of which classes show up
classes <-c(
  "Malacostraca",
  "Insecta",
  "Branchiopoda",
#  "Arthropoda",
 # "Gastropoda",
  "Cestoda",
  "Arachnida",
  "Collembola",
  "Bivalvia",
 # "Copepoda",
  "Diplopoda"
)


# Get the first colors from the Set1 palette
set1_colors <- brewer.pal(8, "Set1")

# Create a named vector
class_colors <- setNames(set1_colors, classes)


# Create the bar plot
ggplot(FOO_All_Samples, aes(x = Type, y = N, fill = class)) + 
  geom_bar(stat = "identity", color = 'black') +
#  geom_text(aes(label = sprintf("%.1f%%", (N / total_samples) * 100)), hjust = -0.1, size = 3.5) + 
  coord_flip() +
  labs(x = "Family", y = "Frequency of Occurrence (Number of Samples)", fill = "Class") + 
  scale_fill_manual(values = class_colors, labels = custom_labels) +  # Apply the custom colors
  theme_classic() +
  theme(axis.text = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),,
        legend.text = element_text(size = 13))

