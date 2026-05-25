#######------------------Ch 1: Fire Freq ITS Data analysis ---------------------#######
#Michelle S. Henson
#6/2023

#Set up---------------------------
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

#######~~~~ITS~~~~~#######

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

###~~~Alpha-diversity~~~~####
###analysis of different alpha diversity indices from https://microbiome.github.io/tutorials/Tutorial.html
ITS.pseq <- ITS.physeq

ITS.alpha <- microbiome::alpha
ITS.tab <- alpha(ITS.pseq, index = "all")
kable(head(ITS.tab))


#Now I need to create a new dataframe with metadata, couldn't figure out how to merge on R so did it manually
##Extract ITS.tab into .csv table to manually add metadata
ITS_alpha_diversity<-read.csv("ITS.tab.csv")

str(ITS_alpha_diversity)
ITS_alpha_diversity$TSF_days <- as.numeric(ITS_alpha_diversity$TSF_days)

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
#nothing


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
#this code has misaligned coefficients, so it doesnt work

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

###~~~Beta-diversity~~~~####
set.seed(43)
# Calculate beta diversity distance matrix
ITS_beta_diversity <- phyloseq::distance(ITS.physeq, method = "bray")

# Perform PCoA (Principal Coordinate Analysis)
ITS.pcoa <- ordinate(ITS.physeq, method = "PCoA", distance = "bray")

# Plot PCoA
plot_ordination(ITS.physeq, ITS.pcoa)


#Run NMDS on distance matrix
ITS.nmds<-metaMDS(ITS_beta_diversity,distance="bray", #use bray-curtis distance
                  k=2, #2 dimensions
                  try=100 #force it to try 100 different times (default is 20, for publication I recommend 500)
)
ITS.nmds #Gives stress value of best solution, aim for <0.2
#stress=0.2366; mehhh

#Check the fit
stressplot(ITS.nmds) #Also called a "Shepard diagram"
#We want this plot to show a monotonic relationship; 0.9 is a good fit; 
#Linear R2=0.724

#Make ordination plot--------------------------
#Using basic vegan graphics
#Site labels only
ordiplot(ITS.nmds, type="text", display="sites")
#Points only
ordiplot(ITS.nmds, type="points", display="sites")

#Can output ordination coordinates and use ggplot for nicer graphics
ITS.nmds.scores <- as.data.frame(scores(ITS.nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
#Add to bacteria metadata info dataframe
ITS.info<-cbind.data.frame(metadata,ITS.nmds.scores)
ITS.info$site <- row.names(metadata)


#Build NMDS figure in ggplot
ggplot(ITS.info, aes(NMDS1, NMDS2)) +
  geom_point() +
  geom_text_repel(aes(NMDS1, NMDS2, label = fire_freq), size = 3)#avoids text overlapping

#Example: Color NMDS plot by different fire freq----------------------
#Revise plotting code from above
ggplot(ITS.info, aes(NMDS1, NMDS2)) +
  geom_point(aes(color=fire_freq, shape=fire_freq)) +
  geom_text_repel(aes(NMDS1, NMDS2, label = site), size = 2)#avoids text overlapping

ggplot(ITS.info, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = fire_freq, shape = field_site), size = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(color = "Fire Regime", shape = "Sites")

#PERMANOVA to test whether composition differs by fire frequency or site---------
adonis(ITS_beta_diversity~fire_freq*field_site, data = ITS.info, permutations=999)
#fire_freq, field_site, and interaction were sign.

adonis(ITS_beta_diversity~fire_freq*P1, data = ITS.info, permutations=999)
#nope

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

###Indicator Species Analysis for ITS###
#######Indicator Species Analysis#######
library(indicspecies)

ITS.taxonomy$ASVs <- rownames(ITS.taxonomy)
# Step 2: Remove row names from the dataframe
rownames(ITS.taxonomy) <- NULL

head(ITS.taxonomy)

# Write out your phyloseq OTU table and export it
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

#Now trying to join ISA dataframe with taxonomy table
ITS.PBind.df <- data.frame(
  ASV = rownames(ITS.SPind$sign),  # Extract ASV names
  stat = ITS.SPind$sign$stat,      # Extract statistic values
  p.value = ITS.SPind$sign$p.value # Extract p-values
)

str(ITS.PBind.df)
head(ITS.PBind.df)



#######Functional Guilds#######
library(fungaltraits)
fungal_traits()

fun2fun<-read.csv("funtothefun.csv")
fun2fun <- tidyr::spread(fun2fun,key = trait_name,value = value,convert = TRUE)

write.csv(ITS.taxonomy1, file = 'ITS_taxonomy.csv', row.names = TRUE)

#Used this youtube to help assign: https://www.youtube.com/watch?v=87AWRkxq8y4
write.csv(funguild, file = 'funguild.csv', row.names = TRUE)


##Merging my functional guilds data from FungalTraits with my metadata and OTU table
funct_guilds <- read.csv("funct_guilds.csv")
ITS.SVs<-read_qza("feature_table4_ch1fungi.qza")
ITS.SVs_df <- as.data.frame(ITS.SVs$data)

write.csv(ITS.SVs_df, file = 'ITS.SVs_df.csv', row.names = TRUE)


ITS.SVs_df$ASVs <- rownames(ITS.SVs_df)
rownames(ITS.SVs_df) <- NULL

merged_data <- left_join(funct_guilds, ITS.SVs_df, by = "ASVs")

svs_tax <- left_join(ITS.SVs_df, ITS.taxonomy1, by = "ASVs")

#Trying to merge together my ITS.SVs_df (OTU table) and taxonomy table based on ASVs
merged_data <- left_join(ITS.SVs_df, ITS.taxonomy1, by = "ASVs")

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

# Rename the column "otu_trt_agg_data[[1]]" to "fire_freq"
colnames(summarized_data)[1] <- "fire_freq"

# Transpose from long to wide format
long_data <- summarized_data %>%
  pivot_longer(cols = -fire_freq, 
               names_to = "ASV_column", 
               values_to = "sum_value")

# Remove 'X' prefix from ASV names in ASV_column
long_data <- long_data %>%
  mutate(ASV_column = gsub("^X", "", ASV_column))

colnames(long_data)[colnames(long_data) == "ASV_column"] <- "ASVs"

colnames(long_data)[colnames(long_data) == "sum_value"] <- "abundance"

#next add funct_guild info to my long dataset (treatment, ASV, abundance(sum))
merged_data_final <- merge(long_data, funct_guilds, by = "ASVs", all.x = TRUE)

# Replace NA values in Genus column with 'unidentified'
merged_data_final <- merged_data_final %>%
  mutate(Genus = coalesce(Genus, 'Unassigned'))


######Chi-square test#######

#clean up my dataframe to get it into a useable format for chi sq test
#Subset the dataframe to include only relevant columns
merged_data_subset <- merged_data_final[, c("ASVs", "fire_freq", "abundance", "guild")]

#Summarize the data into a contingency table
#Create contingency table by summing abundance across guilds and fire frequency treatments
contingency_table <- xtabs(abundance ~ guild + fire_freq, data = merged_data_subset)

chisq <- chisq.test(contingency_table)
chisq
#x-squared= 59452, df=20, p = < 0.0001
