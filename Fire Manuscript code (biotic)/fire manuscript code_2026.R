#######------------------Fire Freq ITS, 16S, and Plant Data Analysis ---------------------#######
#Michelle S. Henson
#6/2023-2026

#Set up
#Set working directory, clean environment, etc.

library(phyloseq)
library(ggplot2)
library(ape)
library(vegan)
library(reshape2)
library(dplyr)
library(tidyr)
library(qiime2R)
library(ggpubr)
library(ggrepel)
library(GGally)
library(factoextra)
library(tidyverse)
library(performance)
library(broom)
library(emmeans)
library(lmerTest)
library(car)
library(rstatix)
library(microbiome)
library(knitr)
library(DESeq2)

#***FUNGAL COMMUNITY ANALYSES (ITS)***----------------------------------------

#feature table
ITS.SVs<-read_qza("feature_table4_ch1fungi.qza")

#metadata
metadata<-read.table("ch1metadata_info.csv")

metadata <- readr::read_csv("ch1metadata_info.csv")
head(metadata)

#taxonomy
ITS.taxonomy<-read_qza("taxonomy4ch1fungi_unite.qza")
head(ITS.taxonomy$data)
ITS.taxonomy<-parse_taxonomy(ITS.taxonomy$data)
head(ITS.taxonomy)

#creating a phyloseq object
ITS.physeq<-qza_to_phyloseq(
  features="feature_table4_ch1fungi.qza",
  tree="rooted-tree.qza",
  "taxonomy4ch1fungi_unite.qza",
  metadata = "ch1metadata.txt")


##Alpha Diversity -------------------------------------------------------

###analysis of different alpha diversity indices from https://microbiome.github.io/tutorials/Tutorial.html
ITS.pseq <- ITS.physeq

ITS.alpha <- microbiome::alpha
ITS.tab <- alpha(ITS.pseq, index = "all")
kable(head(ITS.tab))


#Now create a new dataframe with metadata, couldn't figure out how to merge on R so did it manually
##Extract ITS.tab into .csv table to manually add metadata
ITS_alpha_diversity<-read.csv("ITS.tab.csv")

str(ITS_alpha_diversity)
ITS_alpha_diversity$TSF_days <- as.numeric(ITS_alpha_diversity$TSF_days)

###Observed Richness---------------------------------------------

#Try quick linear model for observed richness
m4<-lm(observed~fire_freq*field_site, data = ITS_alpha_diversity)
hist(m4$residuals)#left skewed
qqPlot(m4)
plot(m4)
Anova(m4)
#no sign. diff's

#TSF (1 year vs two months)
m4.1<-lm(observed~fire_freq*TSF, data = ITS_alpha_diversity)
hist(m4.1$residuals)
qqPlot(m4.1)
summary(m4.1)
Anova(m4.1)
#no sign.
#this code has misaligned coefficients, so it doesnt work

m4.2<-lm(observed~fire_freq*TSF_days, data = ITS_alpha_diversity)
hist(m4.2$residuals)
qqPlot(m4.2)
summary(m4.2)
Anova(m4.2)
#changed TSF_days to a number
#no sign.

m4.3<-lm(observed~fire_freq + TSF_days, data = ITS_alpha_diversity)
hist(m4.3$residuals)
qqPlot(m4.3)
summary(m4.3)
Anova(m4.3)
#use two main effects (+); not sign.

m4.4<-lm(observed~fire_freq + TSF, data = ITS_alpha_diversity)
hist(m4.4$residuals)
qqPlot(m4.4)
summary(m4.4)
Anova(m4.4)
#no main effects

m4.4.1<-lm(observed~P1, data = ITS_alpha_diversity)
Anova(m4.4.1)
#no effects
m4.4.1<-lm(observed~P2, data = ITS_alpha_diversity)
Anova(m4.4.1)
#no effects

###Shannon Diversity---------------------------------------------
#Shannon diversity
m11<-lm(diversity_shannon~fire_freq*field_site, data = ITS_alpha_diversity)
hist(m11$residuals) #normalish, kinda left skewed tho
qqPlot(m11)
plot(m11)
Anova(m11)
#no sign. diff's

#TSF (1 year vs two months)
m5.1<-lm(diversity_shannon~fire_freq*TSF, data = ITS_alpha_diversity)
hist(m5.1$residuals)
qqPlot(m5.1)
summary(m5.1)
Anova(m5.1)
#no sign.
#this code has misaligned coefficients, so it doesnt work

m5.2<-lm(diversity_shannon~fire_freq*TSF_days, data = ITS_alpha_diversity)
hist(m5.2$residuals)
qqPlot(m5.2)
summary(m5.2)
Anova(m5.2)
#changed TSF_days to a number
#no sign.

m5.3<-lm(diversity_shannon~fire_freq + TSF_days, data = ITS_alpha_diversity)
hist(m5.3$residuals)
qqPlot(m5.3)
summary(m5.3)
Anova(m5.3)
#use two main effects (+); not sign.

m4.5.1<-lm(diversity_shannon~P1, data = ITS_alpha_diversity)
Anova(m4.5.1)
#no sign.

###Pielou's Evenness---------------------------------------------

#Try quick linear model for pielou's evenness
m12<-lm(evenness_pielou~fire_freq*field_site, data = ITS_alpha_diversity)
hist(m12$residuals)#left skewedish
qqPlot(m12)
plot(m12)
Anova(m12)
#no sign. diffs

#TSF (1 year vs two months)
m6.1<-lm(evenness_pielou~fire_freq*TSF, data = ITS_alpha_diversity)
hist(m6.1$residuals)
qqPlot(m6.1)
summary(m6.1)
Anova(m6.1)
#no sign.


m6.2<-lm(evenness_pielou~fire_freq*TSF_days, data = ITS_alpha_diversity)
hist(m6.2$residuals)
qqPlot(m6.2)
summary(m6.2)
Anova(m6.2)
#changed TSF_days to a number
#no sign.

m6.3<-lm(evenness_pielou~fire_freq + TSF_days, data = ITS_alpha_diversity)
hist(m6.3$residuals)
qqPlot(m6.3)
summary(m6.3)
Anova(m6.3)
#use two main effects (+); not sign.

m4.6.1<-lm(evenness_pielou~P2, data = ITS_alpha_diversity)
Anova(m4.6.1)
#nothing

##Beta Diversity --------------------------------------------------------

set.seed(431)
#Calculate beta diversity distance matrix
ITS_beta_diversity <- phyloseq::distance(ITS.physeq, method = "bray")

#Perform PCoA (Principal Coordinate Analysis)
ITS.pcoa <- ordinate(ITS.physeq, method = "PCoA", distance = "bray")

#Plot PCoA
plot_ordination(ITS.physeq, ITS.pcoa)

###NMDS -----------------------------------------------------------------

set.seed(112)

ITS_beta_diversity <- phyloseq::distance(ITS.physeq, method = "bray")

ITS.pcoa <- ordinate(ITS.physeq, method = "PCoA", distance = "bray")
plot_ordination(ITS.physeq, ITS.pcoa)

#Run NMDS (3 dimensions, 500 tries)
ITS.nmds <- metaMDS(
  ITS_beta_diversity,
  distance = "bray",
  k = 3,       #3 dimensions
  try = 500)

ITS.nmds 
stressplot(ITS.nmds)  

#Extract NMDS scores and combine with metadata
ITS.nmds_scores <- as.data.frame(scores(ITS.nmds))
ITS.info3d <- cbind(metadata, ITS.nmds_scores)
ITS.info3d$site <- row.names(metadata)


#Build NMDS figure in ggplot
ggplot(ITS.info3d, aes(NMDS1, NMDS2)) +
  geom_point() +
  geom_text_repel(aes(NMDS1, NMDS2, label = fire_freq), size = 3)

#Color NMDS plot by different fire freq
#Revise plotting code from above
ggplot(ITS.info3d, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=fire_freq, shape=fire_freq)) +
  geom_text_repel(aes(NMDS1, NMDS2, label = site), size = 2)

ggplot(ITS.info3d, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = fire_freq, shape = field_site), size = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(color = "Fire Regime", shape = "Sites")



###PERMANOVA ------------------------------------------------------------

###PERMANOVA to test whether composition differs by fire frequency or site
adonis(ITS_beta_diversity~fire_freq*field_site, data = ITS.info, permutations=999)
#fire_freq, field_site, and interaction were sign.

adonis(ITS_beta_diversity~fire_freq*P1, data = ITS.info, permutations=999)
#not sign.

adonis(ITS_beta_diversity~fire_freq*TSF, data = ITS.info, permutations=999)
#TSF sign.

adonis(ITS_beta_diversity~fire_freq*field_site, data = ITS.info, permutations=9999)
#more statistical accuracy by a power of 10; fire freq, field site, and their interaction were sign.
#they were all significant!!

ITS.centroid.fire_freq <- ITS.info %>%
  group_by(fire_freq) %>%
  summarize(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

ggplot(ITS.info, aes(NMDS1, NMDS2, color = fire_freq)) +
  geom_point() +
  stat_ellipse(show.legend=FALSE) +
  geom_point(data=ITS.centroid.fire_freq, size=3, shape=21, color="black",
             aes(fill=fire_freq), show.legend=FALSE) + 
  labs(color = "Fire Frequency")
#centroid for fire_freq

ITS.centroid.site <- ITS.info %>%
  group_by(field_site) %>%
  summarize(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

ggplot(ITS.info, aes(NMDS1, NMDS2, color = field_site)) +
  geom_point() +
  stat_ellipse(show.legend=FALSE) +
  geom_point(data=ITS.centroid.site, size=3, shape=21, color="black",
             aes(fill=field_site), show.legend=FALSE) +
  labs(color = "Field Site")
#centroid for site

###Beta Dispersion ------------------------------------------------------

#Beta-dispersion
ITS.bd <- betadisper(ITS_beta_diversity, ITS.info$fire_freq)
#Beta_disper
anova(ITS.bd)
permutest(ITS.bd)
#sign.

ITS.bd.site <- betadisper(ITS_beta_diversity, ITS.info$field_site)
#Beta_disper
anova(ITS.bd.site)
permutest(ITS.bd.site)
#not sign.
#look at vegan vignette and see what they use;
#might want to use permutest for dispersion

##Indicator Species Analysis (ISA) -------------------------------------

library(indicspecies)

ITS.taxonomy$ASVs <- rownames(ITS.taxonomy)

#Remove row names from the dataframe
rownames(ITS.taxonomy) <- NULL

head(ITS.taxonomy)

#Write out your phyloseq OTU table and export it
write.csv(ITS.physeq@otu_table,'ITS.M3.SP.otus.csv')

#Import phyloseq OTU table as an OTU table/dataframe
ITS.SpOTU<-read.csv('ITS.M3.SP.otus.csv')

#do some shuffling of the OTU table
ITS.SpOTUFlip <- as.data.frame(t(ITS.SpOTU)) #makes it a dataframe and puts x into y and y into x (flips it)
names(ITS.SpOTUFlip) <- as.matrix(ITS.SpOTUFlip[1, ]) #renames columns
ITS.SpOTUFlip<- ITS.SpOTUFlip[-1, ] #removes first row
ITS.SpOTUFlip_num<-as.data.frame(lapply(ITS.SpOTUFlip, as.numeric)) #convert from character to number
ITS.SpOTUFlip_num$sample_id<-row.names(ITS.SpOTUFlip) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like


#read in metadata
metadata<-read.csv("ch1metadata_info.csv") #read in metadata
head(metadata) # check

#join based on sample IDs, assuming they're the same for both OTU table and metadata
ITS.SpOTU_Final<-left_join(ITS.SpOTUFlip_num, metadata, by = c("sample_id" = "sample_id"))
head(ITS.SpOTU_Final)

ITS.SPotus = ITS.SpOTU_Final[,1:6922]#select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
ITS.SPfire = ITS.SpOTU_Final$fire_freq #the metadata column group you care about
ITS.SPsite = ITS.SpOTU_Final$field_site

ITS.SPind=multipatt(x=ITS.SPotus, cluster=ITS.SPfire, func = "r.g", control = how(nperm=9999))
ITS.SPind.site=multipatt(x=ITS.SPotus, cluster=ITS.SPsite, func = "r.g", control = how(nperm=9999))


summary(ITS.SPind)
summary(ITS.SPind.site)
#gives ASV name, not taxonomic name

#Now join ISA dataframe with taxonomy table
ITS.PBind.df <- data.frame(
  ASV = rownames(ITS.SPind$sign), 
  stat = ITS.SPind$sign$stat,     
  p.value = ITS.SPind$sign$p.value)

str(ITS.PBind.df)
head(ITS.PBind.df)


##Functional Guilds ----------------------------------------------------
library(fungaltraits)
fungal_traits()

fun2fun<-read.csv("funtothefun.csv")
fun2fun <- tidyr::spread(fun2fun,key = trait_name,value = value,convert = TRUE)

write.csv(ITS.taxonomy, file = 'ITS_taxonomy.csv', row.names = TRUE)

##Merging my functional guilds data from FungalTraits with my metadata and OTU table
funct_guilds <- read.csv("funct_guilds.csv")
ITS.SVs<-read_qza("feature_table4_ch1fungi.qza")
ITS.SVs_df <- as.data.frame(ITS.SVs$data)

write.csv(ITS.SVs_df, file = 'ITS.SVs_df.csv', row.names = TRUE)


ITS.SVs_df$ASVs <- rownames(ITS.SVs_df)
rownames(ITS.SVs_df) <- NULL

merged_data <- left_join(funct_guilds, ITS.SVs_df, by = "ASVs")

svs_tax <- left_join(ITS.SVs_df, ITS.taxonomy, by = "ASVs")

#Trying to merge together my ITS.SVs_df (OTU table) and taxonomy table based on ASVs
merged_data <- left_join(ITS.SVs_df, ITS.taxonomy, by = "ASVs")

#Now I'm trying to merge my data together by ASVs for functional guild
merged_data <- left_join(merged_data, funct_guilds, by = "ASVs")

#Getting rid of the unneccessarty taxonomic information because all we need is Genus
merged_data <- merged_data %>%
  select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus.x, -Species)

#Renaming the Genus.y to just Genus
names(merged_data)[names(merged_data) == "Genus.y"] <- "Genus"

#Tax table + Guild table
#Then add OTU table
#then joined OTU and meta data

#Trying to go from wide data structure to long (Column for sample_Id etc)
write.csv(merged_data, file = 'merged_data.csv', row.names = TRUE)
#transposed it in excel, we'll see if this works. its look wonky

#now upload it
agg_data <- read.csv("agg_merged_data.csv")

agg_merged_data <- left_join(agg_data, metadata, by = "sample_id")
#metadata_columns <- c("sample_id", "field_site", "fire_freq", "plot", "TSF_days", "TSF", "P1", "P2")
#merged_data1 <- agg_merged_data[, c(metadata_columns, setdiff(names(agg_merged_data), metadata_columns))]

#Tax table: ITS.taxonomy1 = ASVs
#Guild table: funct_guild = ASVs
#OTU table: ITS.SVs_df = ASVs, sample_id/DL01
#metadata: sample_id, fire_freq, etc

otu_trt_agg_data <- read.csv("otu_trt_merged_data.csv")
#took the OTU table (ITS.SVs), transposed in excel, then added fire treatment in excel

#otu_trt_agg_dataframe
#Group by fire_freq and sum for each ASV column (need to add guild info) then transpose (or add guild info after)
#new dataframe <- tidyverse group by column #3-2,### then drop then pipe it into sum equation
str(otu_trt_agg_data)
#grouping fire_freq for sum of each ASV into this dataset
summarized_data <- otu_trt_agg_data %>%
  group_by(otu_trt_agg_data[[1]]) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

#Rename the column "otu_trt_agg_data[[1]]" to "fire_freq"
colnames(summarized_data)[1] <- "fire_freq"

#Transpose from long to wide format
long_data <- summarized_data %>%
  pivot_longer(cols = -fire_freq, 
               names_to = "ASV_column", 
               values_to = "sum_value")

#Remove x prefix from ASV names in ASV column
long_data <- long_data %>%
  mutate(ASV_column = gsub("^X", "", ASV_column))

colnames(long_data)[colnames(long_data) == "ASV_column"] <- "ASVs"

colnames(long_data)[colnames(long_data) == "sum_value"] <- "abundance"

#next add funct_guild info to my long dataset (treatment, ASV, abundance(sum))
merged_data_final <- merge(long_data, funct_guilds, by = "ASVs", all.x = TRUE)

#Replace NA values in Genus column with 'unidentified'
merged_data_final <- merged_data_final %>%
  mutate(Genus = coalesce(Genus, 'Unassigned'))

#calculate relative abundance of each guild within fire treatment
merged_data_final <- merged_data_final %>%
  group_by(fire_freq) %>%
  mutate(
    relative_abundance = abundance / sum(abundance) * 100
  ) %>%
  ungroup()

#summarize relative abundance by guild
merged_data_final <- merged_data_final %>%
  group_by(fire_freq, guild) %>%
  summarise(
    relative_abundance = sum(relative_abundance),
    .groups = "drop"
  )


##### GROUP FUNCTIONAL GUILDS #####

group_guild <- function(guild) {
  case_when(
    guild %in% c(
      "soil_saprotroph",
      "unspecified_saprotroph",
      "litter_saprotroph",
      "wood_saprotroph",
      "dung_saprotroph",
      "pollen_saprotroph",
      "nectar/tap_saprotroph"
    ) ~ "Saprotroph",
    
    guild %in% c("ectomycorrhizal") ~ "Ectomycorrhizal",
    
    guild %in% c(
      "mycoparasite",
      "animal_parasite",
      "algal_parasite",
      "lichen_parasite"
    ) ~ "Parasitic",
    
    guild %in% c("plant_pathogen") ~ "Pathogen",
    
    guild %in% c(
      "foliar_endophyte",
      "root_endophyte"
    ) ~ "Endophyte",
    
    guild %in% c("epiphyte") ~ "Epiphyte",
    
    guild %in% c("sooty_mold") ~ "Mold",
    
    guild %in% c("arbuscular_mycorrhizal") ~ "Arbuscular Mycorrhizal",
    
    guild %in% c("lichenized") ~ "Lichenized",
    
    TRUE ~ "Unassigned")}

#Assign grouped guild names
grouped_data <- merged_data_final %>%
  mutate(
    guild = ifelse(is.na(guild), "Unassigned", guild),
    guild = group_guild(guild))



###Chi-Square Tests -----------------------------------------------------

#clean up my dataframe to get it into a useable format for chi sq test
#Subset the dataframe to include only relevant columns
merged_data_subset <- merged_data_final[, c("ASVs", "fire_freq", "abundance", "guild")]

#Summarize the data into a contingency table
#Create contingency table by summing abundance across guilds and fire frequency treatments
contingency_table <- xtabs(abundance ~ guild + fire_freq, data = merged_data_subset)

chisq <- chisq.test(contingency_table)
chisq
#x-squared= 59452, df=20, p = < 0.0001



#***BACTERIAL COMMUNITY ANALYSES (16S rRNA)***-------------------------------

#feature table
SVs<-read_qza("feature_table5_ch1Whole.qza")

#metadata
metadata<-read.table("ch1metadata_info.csv")

metadata <- readr::read_csv("ch1metadata_info.csv")
head(metadata)
str(metadata)

#taxonomy
bact.taxonomy<-read_qza("taxonomy5ch1Whole_silva.qza")
head(bact.taxonomy$data)
bact.taxonomy<-parse_taxonomy(bact.taxonomy$data)
head(bact.taxonomy)

#creating a phyloseq object
bact.physeq<-qza_to_phyloseq(
  features="feature_table5_ch1Whole.qza",
  tree="rooted-tree5.qza",
  "taxonomy5ch1Whole_silva.qza",
  metadata = "ch1metadata.txt")

##Alpha Diversity -------------------------------------------------------
###analysis of different alpha diversity indices from https://microbiome.github.io/tutorials/Tutorial.html

bact.pseq <- bact.physeq

alpha <- microbiome::alpha
tab <- alpha(bact.pseq, index = "all")
kable(head(tab))

#Now I need to create a new dataframe with metadata, couldn't figure out how to merge on R so did it manually
#write.csv(tab, file = "tab.csv", row.names = FALSE)

bact_alpha_diversity<-read.csv("alphatab.csv")

str(bact_alpha_diversity)
bact_alpha_diversity$TSF_days <- as.numeric(bact_alpha_diversity$TSF_days)

###Observed Richness---------------------------------------------

#Try quick linear model for observed richness
m1<-lm(observed~fire_freq*field_site, data = bact_alpha_diversity)
hist(m1$residuals)#normal
qqPlot(m1)
summary(m1)
Anova(m1)
#no sign. diff's

#TSF (1 year vs two months)
m3<-lm(observed~fire_freq*TSF, data = bact_alpha_diversity)
hist(m3$residuals)
qqPlot(m3)
summary(m3)
Anova(m3)
#no sign.

m4<-lm(observed~fire_freq*TSF_days, data = bact_alpha_diversity)
hist(m4$residuals)
qqPlot(m4)
summary(m4)
Anova(m4)
#changed TSF_days to a number

#Estimated marginal means (slopes) for fire frequency across TSF days
emm <- emtrends(m4, ~ fire_freq, var = "TSF_days")

m5<-lm(observed~fire_freq + TSF_days, data = bact_alpha_diversity)
hist(m5$residuals)
qqPlot(m5)
summary(m5)
Anova(m5)
#use two main effects (+); not sign.


m6<-lm(observed~fire_freq + TSF, data = bact_alpha_diversity)
hist(m6$residuals)
qqPlot(m6)
summary(m6)
Anova(m6)
#no main effects

###Shannon Diversity---------------------------------------------

m2<-lm(diversity_shannon~fire_freq*field_site, data = bact_alpha_diversity)
hist(m2$residuals)#normalish?
qqPlot(m2)
summary(m2)
Anova(m2)
#no sign. diff's

#TSF (1 year vs two months)
m2.1<-lm(diversity_shannon~fire_freq*TSF, data = bact_alpha_diversity)
hist(m2.1$residuals)
qqPlot(m2.1)
summary(m2.1)
Anova(m2.1)
#no sign.

m2.2<-lm(diversity_shannon~fire_freq*TSF_days, data = bact_alpha_diversity)
hist(m2.2$residuals)
qqPlot(m2.2)
summary(m2.2)
Anova(m2.2)
#changed TSF_days to a number
#no sign.

m2.3<-lm(diversity_shannon~fire_freq + TSF_days, data = bact_alpha_diversity)
hist(m2.3$residuals)
qqPlot(m2.3)
summary(m2.3)
Anova(m2.3)
#use two main effects (+); not sign.

m2.4<-lm(diversity_shannon~fire_freq + TSF, data = bact_alpha_diversity)
hist(m2.4$residuals)
qqPlot(m2.4)
summary(m2.4)
Anova(m2.4)
#no main effects

###Pielou's Evenness---------------------------------------------
#Try linear model for pielou's evenness
m3<-lm(evenness_pielou~fire_freq*field_site, data = bact_alpha_diversity)
hist(m3$residuals)#normal
qqPlot(m3)
summary(m3)
Anova(m3)
#field_site p=0.0009

#Pairwise comparisons of field sites for pielou's evenness
emmeans(m3,pairwise~field_site, adjust="Tukey")
#Here we see that the difference (estimate) is largest for 21-Acre and DL for evenness
#DL-FOJ and DL-TOA


#TSF (1 year vs two months)
m3.1<-lm(evenness_pielou~fire_freq*TSF, data = bact_alpha_diversity)
hist(m3.1$residuals)
qqPlot(m3.1)
summary(m3.1)
Anova(m3.1)
#no sign.

m3.2<-lm(evenness_pielou~fire_freq*TSF_days, data = bact_alpha_diversity)
hist(m3.2$residuals)
qqPlot(m3.2)
summary(m3.2)
Anova(m3.2)
#changed TSF_days to a number
#no sign.

m3.3<-lm(evenness_pielou~fire_freq + TSF_days, data = bact_alpha_diversity)
hist(m3.3$residuals)
qqPlot(m3.3)
summary(m3.3)
Anova(m3.3)
#use two main effects (+); not sign.

m3.4<-lm(evenness_pielou~fire_freq + TSF, data = bact_alpha_diversity)
hist(m3.4$residuals)
qqPlot(m3.4)
summary(m3.4)
Anova(m3.4)
#no main effects

##Beta Diversity --------------------------------------------------------

set.seed(23)
#Calculate beta diversity distance matrix
beta_diversity <- phyloseq::distance(bact.physeq, method = "bray")

#Perform PCoA 
pcoa <- ordinate(bact.physeq, method = "PCoA", distance = "bray")

#Plot PCoA
plot_ordination(bact.physeq, pcoa)


#Run NMDS on distance matrix
bact.nmds<-metaMDS(beta_diversity,distance="bray", #use bray-curtis distance
                   k=2, #2 dimensions
                   try=500 #force it to try 100 different times (default is 20, for publication I recommend 500)
)
bact.nmds #Aim for <0.2
#stress=0.133 so we're good

#Check the fit
stressplot(bact.nmds) 
#We want this plot to show a monotonic relationship; 0.9 is a good fit; 
#R2=0.92 it's a good fit!

###NMDS -----------------------------------------------------------------
#Using basic vegan graphics
#Site labels only
ordiplot(bact.nmds, type="text", display="sites")
#Points only
ordiplot(bact.nmds, type="points", display="sites")

#Can output ordination coordinates and use ggplot for nicer graphics
bact.nmds.scores <- as.data.frame(scores(bact.nmds))  
#Add to bacteria metadata info dataframe
bact.info<-cbind.data.frame(metadata,bact.nmds.scores)
bact.info$site <- row.names(metadata)

#Build NMDS figure in ggplot
ggplot(bact.info, aes(NMDS1, NMDS2)) +
  geom_point() +
  geom_text_repel(aes(NMDS1, NMDS2, label = fire_freq), size = 3)

#Color NMDS plot by different fire freq
#Revise plotting code from above
ggplot(bact.info, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=fire_freq, shape=fire_freq)) +
  geom_text_repel(aes(NMDS1, NMDS2, label = site), size = 3)

###PERMANOVA ------------------------------------------------------------

#PERMANOVA to test whether composition differs by fire frequency or site
adonis(beta_diversity~fire_freq*field_site, data = bact.info, permutations=999)
#fire freq, field site, and interaction were sign.

adonis(beta_diversity~fire_freq*field_site, data = bact.info, permutations=9999)
#more statistical accuracy by a power of 10; fire freq, field site, and their interaction were sign.
#they were all significant 

###Beta Dispersion ------------------------------------------------------

#Beta-disper
bact.bd.firefreq <- betadisper(beta_diversity, bact.info$fire_freq)
#Beta_disper
anova(bact.bd.firefreq)
permutest(bact.bd.firefreq)
#not significant

bact.bd.site <- betadisper(beta_diversity, bact.info$field_site)
#Beta_disper
anova(bact.bd.site)
permutest(bact.bd.site)
#not significant

boxplot(bact.bd.firefreq) 
tapply(bact.bd.firefreq$distances, bact.info$fire_freq, mean)

##Indicator Species Analysis (ISA) -------------------------------------

library(indicspecies)

#Import phyloseq OTU table as an OTU table/dataframe
SpOTU<-read.csv('M3.SP.otus.csv')

#do some shuffling of the OTU table
SpOTUFlip <- as.data.frame(t(SpOTU)) #makes it a dataframe and puts x into y and y into x (flips it)
names(SpOTUFlip) <- as.matrix(SpOTUFlip[1, ]) # renames columns
SpOTUFlip<- SpOTUFlip[-1, ] #removes first row
SpOTUFlip_num<-as.data.frame(lapply(SpOTUFlip, as.numeric)) #convert from character to number
SpOTUFlip_num$SampleID<-row.names(SpOTUFlip) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like
SpOTUFlip_num$sample_id<-row.names(SpOTUFlip)

#read in metadata
metadata<-read.csv("ch1metadata_info.csv") #read in metadata
head(metadata) #check

##Join based on SampleID
SpOTU_Final<-left_join(SpOTUFlip_num, metadata, by = c("SampleID" = "sample_id"))
#join based on sample IDs, assuming they're the same for both OTU table and metadata
SpOTU_Final<-left_join(SpOTUFlip_num, metadata, by = c("sample_id" = "sample_id"))
head(SpOTU_Final)

SPotus = SpOTU_Final[,1:19036]#select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
SPfire = SpOTU_Final$fire_freq #the metadata column group you care about
SPsite = SpOTU_Final$field_site

SPind=multipatt(x=SPotus, cluster=SPfire, func = "r.g", control = how(nperm=9999))
SPind.site=multipatt(x=SPotus, cluster=SPsite, func = "r.g", control = how(nperm=9999))

summary(SPind)
summary(SPind.site)
#gives ASV name, not taxonomic name


#***PLANT COMMUNITY ANALYSES***----------------------------------------------

#Load required libraries
library(vegan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(qiime2R)
library(phyloseq)
library(lme4)
library(lmerTest)
library(performance)
library(emmeans)
library(car)
library(ggpubr)
library(patchwork)
library(see)

agsb2021 <- read.csv("asgb2021.csv")
plant_spp <- read.csv("plant_spp.csv")

str(agsb2021)

#Data wrangling
wide2021 <- pivot_wider(agsb2021, names_from=species, values_from=cover_code, values_fill = 0)
#info2020<-wide2021[,c(2,3,6,7)]
#site, plot, type

info2021 <-wide2021[,c(2,3,4)]
spp2021<-wide2021[,c(6:149)]

info2021$plot.name<-paste(info2021$site,info2021$plot,sep=".")
row.names(spp2021)<-info2021$plot.name

agsbplot <- aggregate(agsb2021$species, by = list(agsb2021$site, agsb2021$plot, agsb2021$fire_freq), length)
names(agsbplot) <- c("Site", "Plot", "Type", "Richness")

mean(agsbplot$Richness)

agsbsite <- aggregate(agsbplot$Richness, by = list(agsbplot$Site, agsbplot$Type), mean)
names(agsbsite) <- c("Site", "Type", "Richness")


##Alpha Diversity -------------------------------------------------------
#Calculating Shannon diversity, richness, and evenness
plant.div <- diversity(spp2021, index = "shannon")
plant.rich <- specnumber(spp2021)
plant.even <- diversity(spp2021, index = "shannon") / log(specnumber(spp2021)) 

#combined data set with metadata and calculated values
spp.combined <- cbind(info2021,plant.div, plant.rich, plant.even)

###Observed Richness---------------------------------------------
##Richness
mrich<-lmer(plant.rich~fire_freq*site+(1|plot), data=spp.combined) 
#mixed model with Plot as random factor to control for repeated measures
#year as random factor to control for year effects
#No need for year as continuous predictor because we already know that there are no temporal trends from above analyses
check_model(mrich)#looks great
summary(mrich)
#significant effects of all predictors and interactions
#See results in ANOVA table
Anova(mrich)
#Sign. diff between fire frequency, site for plant species richness

emmeans(mrich,pairwise~fire_freq, adjust="Tukey")
#maintenance higher spp richness
emmeans(mrich,pairwise~fire_freq|site, adjust="Tukey")

mrich1<-lmer(plant.rich~fire_freq*site+(1|plot), data=spp.combined) 
#mixed model with Plot as random factor to control for repeated measures
#year as random factor to control for year effects
#No need for year as continuous predictor because we already know that there are no temporal trends from above analyses
check_model(mrich)#looks great
summary(mrich)
#significant effects of all predictors and interactions
#See results in ANOVA table
Anova(mrich)
#Sign. diff between fire frequency, site for plant species richness

###Shannon Diversity---------------------------------------------

mdiv<-lmer(plant.div~fire_freq*site+(1|plot), data=spp.combined) 
#mixed model with Plot as random factor to control for repeated measures
#year as random factor to control for year effects
#No need for year as continuous predictor because we already know that there are no temporal trends from above analyses
check_model(mdiv)#looks great
summary(mdiv)
#significant effects of all predictors and interactions
#See results in ANOVA table
Anova(mdiv)
#Sign. diff between fire frequency, site

emmeans(mdiv,pairwise~fire_freq, adjust="Tukey")
#maintenance higher diversity
emmeans(mdiv,pairwise~fire_freq|site, adjust="Tukey")

###Pielou's Evenness---------------------------------------------

meven<-lmer(plant.even~fire_freq*site+(1|plot), data=spp.combined) 
#mixed model with Plot as random factor to control for repeated measures
#year as random factor to control for year effects
#No need for year as continuous predictor because we already know that there are no temporal trends from above analyses
check_model(meven)#looks great
summary(meven)
#significant effects of all predictors and interactions
#See results in ANOVA table
Anova(meven)
#Sign. diff between fire frequency and interaction

emmeans(meven,pairwise~fire_freq, adjust="Tukey")
#maintenance higher diversity
emmeans(meven,pairwise~fire_freq|site, adjust="Tukey")

ggline(data=spp.combined, x = "fire_freq", y = "plant.div", add = "mean_se",
       size=0.8,
       color = "site", palette = "jco")

##Beta Diversity --------------------------------------------------------

#NMDS all plots
set.seed(20)
dist1<-vegdist(spp2021, method="bray")

str(spp2021)

###NMDS -----------------------------------------------------------------

###3 dimensions for plants for a better stress value
plant_nmds3 <- metaMDS(dist1, k=3, try=500)
plant_nmds3  #check stress
#stress = 0.177
stressplot(plant_nmds3)
ordiplot(plant_nmds3, type="text", display="sites")

#Extract NMDS scores from the 3D solution
nmds.scores3 <- as.data.frame(scores(plant_nmds3))  # NMDS1, NMDS2, NMDS3
info2021 <- cbind(info2021, nmds.scores3)

#ggplot
#nmds.scores1 <- as.data.frame(scores(plant_nmds1))
#info2021<-cbind.data.frame(info2021,nmds.scores1) 
#Build NMDS figure in ggplot

###PERMANOVA ------------------------------------------------------------

#PERMANOVA
plant_dist <- adonis(dist1~fire_freq*site, data = info2021, permutations=999)

###Force it into 3 dimensions for a better stress value
plant_nmds3 <- metaMDS(dist1, k=3, try=100)
plant_nmds3  #check stress
#stress = 0.177

###Beta Dispersion ------------------------------------------------------

#Run beta dispersion analysis
#Group by fire frequency
plant_bd_firefreq <- betadisper(dist1, group = info2021$fire_freq)

#Test for differences in dispersion (spread)
anova(plant_bd_firefreq)
#not sign.
permutest(plant_bd_firefreq)
#not sign.

boxplot(plant_bd_firefreq)  #or check group averages manually
tapply(plant_bd_firefreq$distances, info2021$fire_freq, mean)

##Pub figure
ggplot(info2021, aes(NMDS1, NMDS2, color = fire_freq, shape = site)) +
  geom_point(size = 3, alpha = 0.8) +  
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  scale_shape_manual(values = c(16, 17, 18, 15)) +  
  labs(color = "Fire Frequency", shape = "Site") +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
  theme_minimal(base_size = 14) +  
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
    axis.title = element_text(face = "bold", size = 14),  
    axis.text = element_text(size = 12),  
    legend.position = "right",  
    legend.text = element_text(size = 12), 
    legend.title = element_text(face = "bold", size = 14))

#just site
ggplot(info2021, aes(NMDS1, NMDS2, color=site)) +
  geom_point(size=2)

#just frequency
ggplot(info2021, aes(NMDS1, NMDS2, color=fire_freq)) +
  geom_point(size=2)


##Indicator Species Analysis (ISA) -------------------------------------

library(indicspecies)
library(dplyr)

#Make sure rows align
all(rownames(spp2021) == info2021$plot.name)  
#should be TRUE

#Metadata with grouping variable
plant.meta <- info2021 %>%
  mutate(
    fire_freq = factor(fire_freq),
    site = factor(site))

set.seed(123)

plant.ind <- multipatt(
  spp2021,
  plant.meta$fire_freq,
  func = "IndVal.g",   
  duleg = TRUE,        
  control = how(nperm = 999))

summary(plant.ind)

summary(plant.ind, alpha = 0.05, indvalcomp = TRUE)

plant.ind.df <- plant.ind$sign %>%
  as.data.frame() %>%
  tibble::rownames_to_column("species") %>%
  filter(p.value <= 0.05) %>%
  arrange(p.value)

plant.ind.df

plant.cover.summary <- agsb2021 %>%
  group_by(species, fire_freq) %>%
  summarise(mean_cover = mean(cover_code), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = fire_freq,
    values_from = mean_cover,
    values_fill = 0)

plant.ind.final <- left_join(
  plant.ind.df,
  plant.cover.summary,
  by = "species")

plant.ind.final

##Control for site effects
library(permute)
library(indicspecies)

set.seed(123)

plant.ind.site <- multipatt(
  spp2021,
  plant.meta$fire_freq,
  func = "IndVal.g",
  duleg = TRUE,
  control = how(
    nperm = 999,
    blocks = plant.meta$site))

summary(plant.ind.site, alpha = 0.05, indvalcomp = TRUE)

write.csv(plant.ind.final, "Supplementary_Table_Indicator_Species.csv", row.names = FALSE)


#PUBLICATION FIGURES----------------------------------------
#Figures for manuscript (minus soil figs; soil data are separate)


###Fig. 2--------------------------------------------
#alpha div richness multi plot

library(ggplot2)
library(patchwork)
library(dplyr)

custom_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "none",
    plot.title = element_blank(),
    plot.tag = element_text(face = "bold", size = 13))

#Create a function to make the 9 plot fig
generate_plot <- function(data, metric, y_label = NULL, x_label = NULL, 
                          sig_data = NULL, remove_x_text = FALSE) {
  
  fire_colors <- c(
    "Annual" = "#EC4E3B",
    "Maintenance" = "#4A91D1")
  
  p <- ggplot(data, aes(x = fire_freq, y = !!sym(metric))) +
    geom_boxplot(outlier.shape = NA, fill = NA, color = "black", alpha = 0.6) +
    geom_jitter(aes(fill = fire_freq), 
                shape = 21, color = "black",
                width = 0.2, size = 2, alpha = 0.9, stroke = 0.5) +
    scale_fill_manual(values = fire_colors) +
    labs(x = x_label, y = y_label) +
    custom_theme
  
  if (remove_x_text) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  if (!is.null(sig_data)) {
    p <- p +
      geom_line(data = sig_data$line, aes(x = x, y = y), linewidth = 0.7, inherit.aes = FALSE) +
      geom_text(data = sig_data$text, aes(x = x, y = y, label = label), 
                fontface = "bold", size = 3.5, inherit.aes = FALSE)
  }
  
  return(p)}

#Adding visual sig data for plants

plant_rich_sig <- list(
  line = data.frame(x = c(1, 2), y = c(38.5, 38.5)),
  text = data.frame(x = 1.5, y = 39, label = "***"))

plant_div_sig <- list(
  line = data.frame(x = c(1, 2), y = c(3.65, 3.65)),
  text = data.frame(x = 1.5, y = 3.67, label = "***"))

plant_even_sig <- list(
  line = data.frame(x = c(1, 2), y = c(0.986, 0.986)),
  text = data.frame(x = 1.5, y = 0.9865, label = "*"))

#Generating the plots

#Fungi
fungi_rich_plot <- generate_plot(
  ITS_alpha_diversity, 
  "observed",
  y_label = "Fungal Richness", 
  remove_x_text = TRUE)

fungi_even_plot <- generate_plot(
  ITS_alpha_diversity, 
  "evenness_pielou",
  y_label = "Fungal Evenness", 
  remove_x_text = TRUE)

fungi_div_plot <- generate_plot(
  ITS_alpha_diversity, 
  "diversity_shannon",
  y_label = "Fungal Diversity", 
  remove_x_text = TRUE)

#Bacteria
bact_rich_plot <- generate_plot(
  bact_alpha_diversity, 
  "observed",
  y_label = "Bacterial Richness", 
  remove_x_text = TRUE)

bact_even_plot <- generate_plot(
  bact_alpha_diversity, 
  "evenness_pielou",
  y_label = "Bacterial Evenness", 
  remove_x_text = TRUE)

bact_div_plot <- generate_plot(
  bact_alpha_diversity, 
  "diversity_shannon",
  y_label = "Bacterial Diversity", 
  remove_x_text = TRUE)

#Plant
plant_rich_plot <- generate_plot(
  spp.combined, 
  "plant.rich", 
  y_label = "Plant Richness", 
  sig_data = plant_rich_sig)

plant_even_plot <- generate_plot(
  spp.combined, 
  "plant.even", 
  y_label = "Plant Evenness", 
  sig_data = plant_even_sig,
  x_label = "Fire Regime")

plant_div_plot <- generate_plot(
  spp.combined, 
  "plant.div",  
  y_label = "Plant Diversity", 
  sig_data = plant_div_sig)

#Combine plots
final_plot <- (
  (fungi_rich_plot | fungi_even_plot | fungi_div_plot) /
    (bact_rich_plot  | bact_even_plot  | bact_div_plot) /
    (plant_rich_plot | plant_even_plot | plant_div_plot)
) +
  plot_annotation(tag_levels = "A")

###Fig. 3--------------------------------------------
#NMDS multiplot


#Define colors for fire regimes
ellipse_colors <- c(
  "Annual" = "#E41A1C",
  "Maintenance" = "#377EB8")

#Fungi
fungi_plot <- ggplot(ITS.info3d, aes(NMDS1, NMDS2, color = fire_freq, shape = field_site)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = ellipse_colors) +
  scale_shape_manual(values = c(16, 17, 18, 15)) +
  stat_ellipse(
    aes(group = fire_freq, color = fire_freq),
    level = 0.95,
    linetype = "solid",
    size = 0.6,
    show.legend = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
    plot.tag = element_text(face = "bold", size = 12)
  ) +
  ggtitle("Fungi") +
  labs(tag = "A")

#bacteria
bacteria_plot <- ggplot(bact.info, aes(NMDS1, NMDS2, color = fire_freq, shape = field_site)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = ellipse_colors) +
  scale_shape_manual(values = c(16, 17, 18, 15)) +
  stat_ellipse(
    aes(group = fire_freq, color = fire_freq),
    level = 0.95,
    linetype = "solid",
    size = 0.6,
    show.legend = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
    plot.tag = element_text(face = "bold", size = 12)) +
  ggtitle("Bacteria") +
  labs(tag = "B")

#plant
plants_plot <- ggplot(info2021, aes(NMDS1, NMDS2, color = fire_freq, shape = site)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = ellipse_colors, name = "Fire Regime") +
  scale_shape_manual(values = c(16, 17, 18, 15), name = "Site") +
  stat_ellipse(
    aes(group = fire_freq, color = fire_freq),
    level = 0.95,
    linetype = "solid",
    size = 0.6) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.box = "vertical",
    legend.key.height = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
    plot.tag = element_text(face = "bold", size = 12)) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)) +
  ggtitle("Plant") +
  labs(tag = "C")

#Combine all three plots in one row
final_plot <- fungi_plot + bacteria_plot + plants_plot +
  plot_layout(ncol = 3) +
  plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5)))

###Fig. 4--------------------------------------------
#total abundance indicator taxa plots

library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(indicspecies)
library(permute)
library(tibble)
library(patchwork)

#Bacteria

#Significant ASVs
significant_ASVs <- rownames(
  SPind$sign[SPind$sign$p.value < 0.05, ])

#Subset abundance table
significant_OTUs <- SPotus[, significant_ASVs]

#Add in our metadata
significant_OTUs$fire_freq <- SPfire

#Convert it to long format
bact.long <- significant_OTUs %>%
  pivot_longer(
    cols = -fire_freq,
    names_to = "ASV",
    values_to = "abundance")

#summarize abundance
bact.summary <- bact.long %>%
  group_by(fire_freq, ASV) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop")

bact.summary$ASV <- gsub("^X", "", bact.summary$ASV)
bact.taxonomy$ASV <- rownames(bact.taxonomy)

bact.plot.data <- bact.summary %>%
  left_join(bact.taxonomy, by = "ASV") %>%
  mutate(
    Phylum = ifelse(is.na(Phylum), "Unclassified", Phylum))

#Fungi

#Significant ASVs
ITS.significant_ASVs <- rownames(
  ITS.SPind$sign[ITS.SPind$sign$p.value < 0.05, ])

#Subset the abundance table
ITS.significant_OTUs <- ITS.SPotus[, ITS.significant_ASVs]

#Add our metadata
ITS.significant_OTUs$fire_freq <- ITS.SPfire

#Convert it to long format
fungi.long <- ITS.significant_OTUs %>%
  pivot_longer(
    cols = -fire_freq,
    names_to = "ASVs",
    values_to = "abundance")

#summarize abundance
fungi.summary <- fungi.long %>%
  group_by(fire_freq, ASVs) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop")

#Clean ASV names
fungi.summary$ASVs <- gsub("^X", "", fungi.summary$ASVs)

#Add taxonomy (don't need this here)
#ITS.taxonomy$ASVs <- rownames(ITS.taxonomy)

fungi.plot.data <- fungi.summary %>%
  left_join(ITS.taxonomy, by = "ASVs") %>%
  mutate(
    Phylum = ifelse(is.na(Phylum), "Unclassified", Phylum))

#Plant
#check row alignment
all(rownames(spp2021) == info2021$plot.name)

#Metadata
plant.meta <- info2021 %>%
  mutate(
    fire_freq = factor(fire_freq),
    site = factor(site))

#Run ISA controlling for site
set.seed(123)

plant.ind.site <- multipatt(
  spp2021,
  plant.meta$fire_freq,
  func = "IndVal.g",
  duleg = TRUE,
  control = how(
    nperm = 999,
    blocks = plant.meta$site))

#Extract sign. indicator taxa
plant.ind.sig <- plant.ind.site$sign %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  filter(p.value <= 0.05)

#Add sci names
plant.taxonomy <- plant_spp %>%
  mutate(
    sci_name = paste(Genus, Species)) %>%
  select(Species.Code, sci_name)

plant.ind.sig.named <- plant.ind.sig %>%
  left_join(
    plant.taxonomy,
    by = c("species" = "Species.Code"))

#Summarize total abundance by fire regime

plant.ind.total <- agsb2021 %>%
  filter(species %in% plant.ind.sig.named$species) %>%
  group_by(fire_freq, species) %>%
  summarise(
    value = sum(cover_code),
    .groups = "drop") %>%
  left_join(
    plant.ind.sig.named %>%
      select(species, sci_name),
    by = "species")

#cllapse taxa contributing <2% total abundance

threshold <- 2

species.totals <- plant.ind.total %>%
  group_by(sci_name) %>%
  summarise(
    total_abundance = sum(value),
    .groups = "drop") %>%
  mutate(
    percent_total = 100 * total_abundance / sum(total_abundance))

keep.species <- species.totals %>%
  filter(percent_total >= threshold) %>%
  pull(sci_name)

other.label <- paste0("Other (< ", threshold, "%)")

plant.ind.total.plot <- plant.ind.total %>%
  mutate(
    sci_name = ifelse(
      sci_name %in% keep.species,
      sci_name,
      other.label)) %>%
  group_by(fire_freq, sci_name) %>%
  summarise(
    value = sum(value),
    .groups = "drop")

#Order legend by total abundance

legend.order <- plant.ind.total.plot %>%
  group_by(sci_name) %>%
  summarise(
    total = sum(value),
    .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(sci_name)

#Move Other to bottom for visualization organization
legend.order <- c(
  setdiff(legend.order, other.label),
  other.label)

plant.ind.total.plot$sci_name <- factor(
  plant.ind.total.plot$sci_name,
  levels = legend.order)

#Color palettes

wes_16 <- c(
  wes_palette("Rushmore1", 5),
  wes_palette("Moonrise3", 4),
  wes_palette("Royal1", 3),
  wes_palette("BottleRocket2", 4))

custom_18_colors <- c(
  wes_palette("Royal1", 4),
  wes_palette("Cavalcanti1", 5),
  wes_palette("GrandBudapest2", 4),
  wes_palette("Rushmore", 4),
  wes_palette("Chevalier1", 4))

darj <- wes_palette("Darjeeling1", 4)

fungi_palette <- c(
  wes_palette("Moonrise2", 4),
  darj[2:4],
  darj[1])

#Bacterial plot

p1 <- ggplot(bact.plot.data, aes(
    x = fire_freq,
    y = total_abundance,
    fill = Phylum)) + geom_bar(
    stat = "identity",
    position = "stack") +
  scale_fill_manual(values = custom_18_colors) +
  scale_y_continuous(
    labels = function(x)
      format(x, scientific = FALSE, big.mark = ",")) +
  labs(x = "Fire Regime", y = "Total Abundance (Reads)",
    fill = "Bacterial Phylum") + theme_minimal()

#Fungal plot

p2 <- ggplot(fungi.plot.data, aes(x = fire_freq, y = total_abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = fungi_palette) + scale_y_continuous(
    labels = function(x)
      format(
        x,
        scientific = FALSE,
        big.mark = ",")) +
  labs( x = "Fire Regime", y = "Total Abundance (Reads)", fill = "Fungal Phylum") + theme_minimal()

#plant plot

p3 <- ggplot(plant.ind.total.plot,
  aes(x = fire_freq,
    y = value,
    fill = sci_name)) + geom_bar(
    stat = "identity",
    position = "stack") +
  scale_fill_manual(values = wes_16) +
  scale_y_continuous(
    labels = function(x)
      format(
        x,
        scientific = FALSE,
        big.mark = ",")) +
  labs(x = "Fire Regime", y = "Total Abundance",
    fill = "Plant Species") + theme_minimal()

#Combined

combined_plot <- p1 + p2 + p3 +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")

combined_plot

###Fig. 5--------------------------------------------
#Mean relative abundance of fungal functional guild plot

guild_relative_abundance <- grouped_data %>%
  group_by(fire_freq, guild) %>%
  summarise(
    relative_abundance = sum(relative_abundance),
    .groups = "drop")

#Order guilds for plotting
guild_relative_abundance$guild <- factor(
  guild_relative_abundance$guild,
  levels = c(
    "Arbuscular Mycorrhizal",
    "Ectomycorrhizal",
    "Endophyte",
    "Epiphyte",
    "Lichenized",
    "Mold",
    "Parasitic",
    "Pathogen",
    "Saprotroph",
    "Unassigned"))

royal1_like_colors <- c("#F21A00", "#0969A2", "#A8C0A1", "lightblue","darkred",
                        "#E1BD6D", "#3E5A5B", "#9986A5", "#B7A399", "#6B8E6F")

#plot
ggplot(guild_relative_abundance,
  aes(
    x = fire_freq,
    y = relative_abundance, fill = guild)) +
  geom_bar(stat = "identity") +
  labs(x = "Fire Regime", y = "Relative Abundance (%)",
    fill = "Functional Guild") +
  scale_fill_manual(values = royal1_like_colors) +
  theme_minimal() +
  theme(legend.position = "right")

###Fig. S1-------------------------------------------
#Mean relative abundance of fungal and bacterial indicator genera

#Bacteria first
#Genus level data for supplementary
#Extract significant bacterial ASVs
significant_ASVs <- rownames(SPind$sign[SPind$sign$p.value < 0.05, ])

#Subset OTU table
significant_OTUs <- SPotus[, significant_ASVs]

#Calculate relative abundance
relative_abundance <- sweep(
  significant_OTUs, 1,
  rowSums(significant_OTUs),
  FUN = "/")

relative_abundance <- as.data.frame(relative_abundance)
relative_abundance$fire_freq <- SPfire

#Mean relative abundance by fire freq
relative_abundance_summary <- relative_abundance %>%
  group_by(fire_freq) %>%
  summarize(across(everything(), mean, na.rm = TRUE))

#Reshape for plotting
relative_abundance_melted <- melt(
  relative_abundance_summary,
  id.vars = "fire_freq")

colnames(relative_abundance_melted)[colnames(relative_abundance_melted) == "variable"] <- "ASVs"

relative_abundance_melted$ASVs <- gsub(
  "^X",
  "", relative_abundance_melted$ASVs)

#Merge taxonomy
bact.taxonomy$ASVs <- rownames(bact.taxonomy)

relative_abundance_with_taxonomy <- merge(relative_abundance_melted, bact.taxonomy, by = "ASVs")

#convert to percentage
relative_abundance_with_taxonomy$value <- relative_abundance_with_taxonomy$value * 100

#replace unidentified taxa
relative_abundance_with_taxonomy$Genus[
  is.na(relative_abundance_with_taxonomy$Genus) |
    grepl("Uncultured", relative_abundance_with_taxonomy$Genus,
      ignore.case = TRUE)] <- "Unassigned"

#Remove underscores
relative_abundance_with_taxonomy$Genus <- gsub(
  "_",
  " ",
  relative_abundance_with_taxonomy$Genus)

#aggregate by genus
relative_abundance_genus <- relative_abundance_with_taxonomy %>%
  group_by(fire_freq, Genus) %>%
  summarize(relative_abundance = sum(value, na.rm = TRUE), .groups = "drop")

#group low abundance taxa
threshold <- 3

relative_abundance_genus <- relative_abundance_genus %>%
  group_by(fire_freq) %>%
  mutate(total_abundance = sum(relative_abundance),
    relative_abundance_perc =
      (relative_abundance / total_abundance) * 100) %>%
  ungroup() %>%
  mutate(Genus = ifelse(relative_abundance_perc < threshold, "Other (< 3%)", Genus)) %>%
  group_by(fire_freq, Genus) %>%
  summarize(relative_abundance = sum(relative_abundance),
    .groups = "drop")

royal1_colors <- c(
  "#4B6A35",
  "#D04F4F",
  "#7CA1B4",
  "#F4C4A1",
  "#6D4687",
  "#F1D500",
  "#1D5B61",
  "#E25F89",
  "#A2A2A1")

#Bacteria plot
bacteria_plot <- ggplot(relative_abundance_genus,
  aes(x = fire_freq, y = relative_abundance,
    fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(
    title = "A", x = "Fire Regime",
    y = "Relative Abundance (%)") +
  scale_fill_manual(values = royal1_colors) +
  theme_minimal()


#Fungal genus level
#Replace missing genera
fungi.plot.data$Genus <- as.character(fungi.plot.data$Genus)

fungi.plot.data$Genus[
  is.na(fungi.plot.data$Genus) |
    grepl(
      "uncultured|unclassified|unidentified",
      fungi.plot.data$Genus,
      ignore.case = TRUE)] <- "Unassigned"

#Remove underscores
fungi.plot.data$Genus <- gsub(
  "_",
  " ",
  fungi.plot.data$Genus)

#Summarize by genus
fungi.genus <- fungi.plot.data %>%
  group_by(fire_freq, Genus) %>%
  summarise(total_abundance = sum(total_abundance),
    .groups = "drop")

#Calculate percent contribution
fungi.genus <- fungi.genus %>%
  group_by(fire_freq) %>%
  mutate(percent_abundance = 100 * total_abundance / sum(total_abundance)) %>%
  ungroup()

genus.threshold <- 3

#Collapse low abundance genera
fungi.genus <- fungi.genus %>%
  mutate(Genus = ifelse(percent_abundance < genus.threshold,
      paste0("Other (< ", genus.threshold, "%)"),
      Genus)) %>%
  group_by(fire_freq, Genus) %>%
  summarise(
    total_abundance = sum(total_abundance),
    .groups = "drop")

#Genus level color palettes
bacteria_genus_palette <- c(
  wes_palette("Royal1", 4),
  wes_palette("Cavalcanti1", 5),
  wes_palette("GrandBudapest2", 4),
  wes_palette("Rushmore", 4),
  wes_palette("Chevalier1", 4))

fungi_genus_palette <- c(
  wes_palette("Moonrise2", 4),
  darj[2:4],
  darj[1],
  wes_palette("BottleRocket2", 3))

#Bacterial genera plot
supp_p1 <- ggplot(
  relative_abundance_genus,
  aes(x = fire_freq, y = relative_abundance,
    fill = Genus)) + geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = bacteria_genus_palette) + scale_y_continuous(labels = function(x)format(
        x, scientific = FALSE,
        big.mark = ",")) +
  labs(
    x = "Fire Regime",
    y = "Total Abundance (Reads)",
    fill = "Bacterial Genus") +
  theme_minimal()

#fungal genera plot
supp_p2 <- ggplot(fungi.genus, aes(x = fire_freq, y = total_abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = fungi_genus_palette) +
  scale_y_continuous(labels = function(x)
      format(x, scientific = FALSE,
        big.mark = "," )) +
  labs(x = "Fire Regime", y = "Total Abundance (Reads)",
    fill = "Fungal Genus") + theme_minimal()

supp_combined_plot <- supp_p2 + supp_p1 + 
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

supp_combined_plot
