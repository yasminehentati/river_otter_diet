# rarefaction curves using inext 
# code adapted from Sam Kreling 


library(pacman)
p_load(here, readr, dplyr, ggplot2, treemap, tidyr, vegan, MicroNiche, RInSp, ggforce, goeveg, strex, glmmTMB)

p_load("dplyr", "tidyr", "here", "ggplot2", "randomcoloR", "stringr", "remotes", 'grid', "iNEXT")

# we'll want to make a rarefaction curve for each latrine 

# Matrix will be for each family with rows as each dietary species and columns 
# for each scat from any individual within that family. The values should be binary (presence/absence).
# This will allow us to check to see if we can make comparisons between different latrines.


# first we need to set up a data frame formatted by latrine / diet species 

Unique_Combo <- read_csv(here("clean_data", "all_diet_cov_data.csv")) %>%
  select(1,3,22,34) # only need sample ID, abundance, family, site 

# fill values greater than 0 with 1 
Unique_Combo$abundance_avg[Unique_Combo$abundance_avg>0] <-1

# now we need to filter out all sites with less than 5 scats
# otherwise we'll have a hard time calculating rarefaction 


# count unique scats per site
site_counts <- Unique_Combo %>%
  group_by(site) %>%
  summarise(unique_samples = n_distinct(sample_ID))

site_counts

# identify sites with fewer than 5 unique sample IDs
removed_sites <- site_counts %>%
  filter(unique_samples < 5) %>%
  pull(site)

# filter dataset to keep only sites with at least 5 unique sample IDs 
Unique_Combo <- Unique_Combo %>%
  filter(site %in% site_counts$site[site_counts$unique_samples >= 5])

# unique(Unique_Combo$site)
#Complete cases
# Unique_Combo <- Unique_Combo %>% drop_na(family)
Unique_Combo$family <- ifelse(is.na(Unique_Combo$family), "unknown", Unique_Combo$family)

# Pivot so that we have each species by each scat sample
## values_fn summarizes the values for any duplicate listings

species.by.scat <- tidyr::pivot_wider(data = Unique_Combo, 
                                      values_from="abundance_avg", 
                                      values_fn=sum, 
                                      names_from = c("family"), 
                                      values_fill = 0)
colnames(species.by.scat)

# now pivot longer all values from site to end of family list **Note** may have to change this range
species.by.scat <- tidyr::pivot_longer(data=species.by.scat, cols=c(3:ncol(species.by.scat))) 

# Remove all spaces from species names -- only need this if doing by species level
# species.by.scat$processed_taxon <- str_replace_all(species.by.scat$processed_taxon, " ", "")

# Convert 2's to 1's this is from the value_fn=sum
species.by.scat$value <- ifelse(species.by.scat$value>1, 1, species.by.scat$value)

# Pivot wider again
species.by.scat <- pivot_wider(data=species.by.scat, values_from = 'value', names_from = c("sample_ID"), values_fill = 0)

# Make sure the numbers are numeric and not being read as characters
# species.by.scat <- species.by.scat %>% mutate_at(c(34:ncol(species.by.scat)), as.numeric)

# funciton to remove 0 sum columns from elements of list that is matrix

remove_zero_sum_columns <- function(mat) {
  cols_to_remove <- colSums(mat) == 0
  return(mat[, !cols_to_remove])
}

# create a list of matrices for each latrine 
mat.list <- list()

for(i in 1:length(unique(species.by.scat$site))){
  
  mat.list[[i]] <- as.matrix(species.by.scat[species.by.scat$site==unique(species.by.scat$site)[[i]],34:ncol(species.by.scat)])
  
  rownames(mat.list[[i]]) <- species.by.scat$name[species.by.scat$site==unique(species.by.scat$site)[[i]]]
  
  mat.list[[i]] <- remove_zero_sum_columns(mat.list[[i]])
}

names(mat.list) <- unique(species.by.scat$site)

# check how many observations per site 
mat.list2 <- mat.list[sapply(mat.list, is.matrix)] 

######### 
# visualize accumulation curves 

#install_version("iNEXT", version = "2.0.20", repos = "http://cran.us.r-project.org")
#packageVersion("iNEXT")
#library(iNEXT)
## Dead end error given in the newer versions of this package

palette <- (distinctColorPalette(length(mat.list2)))

# Hill-Shannon Diversity
inext.out <- iNEXT(mat.list, q=1,datatype = "incidence_raw", endpoint=100)

# head(inext.out)

# Plot with type = 1, with species diversity and number of sampling units as X
g <- ggiNEXT(inext.out, type=1, facet.var="None", se=F)+theme_classic() +
  scale_fill_manual(values=palette)+ scale_shape_manual(values=rep(1,length(palette)))+
  ylim(0,30)+scale_shape_manual(values=c(rep(20,length(palette))))+
  scale_colour_manual(values=c(palette))+
  labs(x="Number of Scat Samples", y="Hill-Shannon Diversity")+
  scale_x_continuous(breaks=seq(0,50,by=5))

g2 <- ggplot_build(g + theme(legend.position ="none"))
g2 <- ggplot_build(g)
g2$data[[1]]$size <- 3
g2$data[[2]]$linewidth <- .6

grid.draw(ggplot_gtable(g2))

# Hill-Simpson Diversity plot 
inext.out <- iNEXT(mat.list, q=2,datatype = "incidence_raw", endpoint=50)

# Plot with type = 1, with species diversity and number of sampling units as X
g <- ggiNEXT(inext.out, type=1, facet.var="None", se=F)+theme_classic() +
  scale_fill_manual(values = palette) + 
  scale_shape_manual(values=rep(1,length(palette)))+ylim(0,40)+
  scale_shape_manual(values=c(rep(20,length(palette))))+
  scale_colour_manual(values=c(palette))+
  labs(x="Number of Scat Samples", 
                                              y="Hill-Simpson Diversity")+
  scale_x_continuous(breaks=seq(0,50,by=5))

g2 <- ggplot_build(g + theme(legend.position = "none"))

g2$data[[1]]$size <- 3
g2$data[[2]]$linewidth <- .6

grid.draw(ggplot_gtable(g2))



#########
# calculate hill-shannon and hill-simpson diversity 
# again here we only have sites with at least 5 samples 

D <- estimateD(mat.list2, datatype="incidence_raw", base='size')

#write.csv(D, here::here("Sequence_Data", "DiversityIndices_SeattleOnly.csv"))
HillSh <- D[D$order==1,]
mean(HillSh$qD) #mean HShD is 14.298
sd(HillSh$qD) #sd of HShD is 2.856 

#write.csv(HillSh,("WiS_HillSh_10Samp.csv"))

HillSi <- D[D$order==2,]
mean(HillSi$qD) #mean HSiD is 12.378
sd(HillSi$qD) # sd of HSiD is 2.549

SpRich <- D[D$order==0,]
mean(SpRich$qD) #mean SR is 16.917
sd(SpRich$qD) #sd of SR is 3.191 

########

# let's check for differences among sites 
# to do this we'll use a permanova
# but first we need to resample our sites to equal sizes
colnames(species.by.scat)


resample_data <- function(data, sample_size = 5) {
  data %>%
    group_by(site) %>%
    filter(n() >= sample_size) %>% 
    sample_n(sample_size, replace = TRUE) %>%  # Sample with replacement
    ungroup() %>%
    filter(rowSums(select(., 3:ncol(.)), na.rm = TRUE) > 0)  # Ensure sampled rows have at least one nonzero value
}


# Number of PERMANOVA runs
n_runs <- 1000
permanova_results <- vector("list", n_runs)

for (i in 1:n_runs) {
  resampled_data <- resample_data(species.by.scat, sample_size = min)
  dist <- vegdist(as.matrix(resampled_data[,4:105]), method="jaccard", binary=T)
  # Adjust the formula according to your actual dataset
  permanova_results[[i]] <- adonis2(dist ~ FoodDesert, data = resampled_data, method = "jaccard", binary=T)
}

View(all_diet)
# Run PERMANOVA multiple times with resampled data
for (i in 1:n_runs) {
  resampled_data <- resample_data(species.by.scat, sample_size = 5)
  
  # Check for and remove completely empty rows (all-zero species data)
  valid_rows <- rowSums(resampled_data[, 3:ncol(resampled_data)], na.rm = TRUE) > 0
  resampled_data <- resampled_data[valid_rows, ]
  
  # Only proceed if there are at least two valid rows
  if (nrow(resampled_data) > 1) {
    dist <- vegdist(as.matrix(resampled_data[, 3:ncol(resampled_data)]), method = "jaccard", binary = TRUE)
    
    permanova_results[[i]] <- adonis2(dist ~ site, 
                                      data = resampled_data,
                                      method = "jaccard", 
                                      binary = TRUE)
  } else {
    permanova_results[[i]] <- NA  # Mark as NA if resampled data is insufficient
  }
}

head(permanova_results)


# Extract and average the R2 values or p-values across runs
r2.df <- matrix(NA, nrow=n_runs, ncol=3)
p.df <- matrix(NA, nrow=n_runs, ncol=1)


for(i in 1:n_runs){
  r2.df[i,] <- permanova_results[[i]]$R2
  p.df[i,] <- permanova_results[[i]]$`Pr(>F)`[1]
}

avg.r2 <- colMeans(r2.df)
# 0.617  first value: R² for the predictor (site) -- proportion of variation explained by site
# 0.383 second value: R² for residuals → unexplained variation (error)
# 1.00 third value: Total R² → Should sum to 1 (100% of variation explained)

avg.p <- mean(p.df)
# p = 0.001536 
