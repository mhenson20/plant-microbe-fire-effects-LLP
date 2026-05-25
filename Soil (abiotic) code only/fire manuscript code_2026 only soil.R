#######------------------Soil Abiotic Data Analysis ---------------------#######
#Michelle S. Henson
#6/2021

#Set up
#Set working directory, clean environment, etc.

#Load packages
library(GGally)
library(factoextra)
library(tidyverse)
library(ggplot2)

#***SOIL NUTRIENT ANALYSES***----------------------------------------

#Import data
soil_gwc <- read.csv("soil_gwc.csv")

str(soil_gwc)

soil_gwc$plot <- as.character(soil_gwc$plot)

#1) What is the variation in soil characteristics measured among the 2 different fire frequencies?
#Use ggpairs to examine correlations among all the response variables in the dataset

#Examine correlations among the body shape response variables
#ggpairs plots for fire freq
ggpairs(soil_gwc, mapping = aes(color = fire_freq), 
        columns = c("gwc_percent", "OM", "P1", 
                    "P2", "K", "Mg", "Ca", "pH", "NO3", "CEC"))+
  scale_colour_manual(values = c("blue","darkorange")) +
  scale_fill_manual(values = c("blue","darkorange"))

#ggpairs plots for site
ggpairs(soil_gwc, mapping = aes(color = site), 
        columns = c("gwc_percent", "OM", "P1", 
                    "P2", "K", "Mg", "Ca", "pH", "NO3", "CEC"))+
  scale_colour_manual(values = c("purple","darkorange", "blue")) +
  scale_fill_manual(values = c("purple","darkorange", "blue"))

#2) Analyze all response variables using PCA with fire frequency
soil_gwc.pca <- prcomp (~ gwc_percent + OM + P1 + P2 + K + Mg + Ca + pH + NO3 + CEC,
                        data=soil_gwc,
                        scale. = TRUE) 

#Get factor loadings on principle components
soil_gwc.pca
#Visualize how much variation is explained by each principle component (Scree plot)
#Basic plot
plot(soil_gwc.pca)#however, this gives us the absolute variances, not % of total variance
#Use function from a different package (factoextra) for nicer plot
fviz_eig(soil_gwc.pca,addlabels = TRUE)
#PC1 and PC2 combined explain 74% of total variation
#PCA 1 explained ~52.7% of residual variances; explains the most variance
#PCA 1 and 2 explain 74% of the variances

#3) Visualize results in a biplot by fire frequency
fviz_pca_biplot(soil_gwc.pca, label = "var", 
                col.ind = soil_gwc$fire_freq, 
                palette = c("darkorange","blue"), 
                col.var = "black", repel = TRUE,
                legend.title = "Fire Frequency")

#Visualize results in a biplot by site
fviz_pca_biplot(soil_gwc.pca, label = "var", 
                col.ind = soil_gwc$site, palette = c("darkorange","blue", "purple"), 
                col.var = "black", repel = TRUE,
                legend.title = "Site")
#try to add shape into this code

fviz_pca_ind(soil_gwc.pca, label = "var", 
             habillage = soil_gwc$fire_freq, 
             col.var = "black", repel = TRUE,
             legend.title = "Fire Frequency",
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark1")
#adds ellipses to fire frequency; but no variables

fviz_pca_ind(soil_gwc.pca, label = "var", 
             habillage = soil_gwc$site, 
             col.var = "black", repel = TRUE,
             legend.title = "Site",
             addEllipses=TRUE, ellipse.level=0.95, palette = "Dark1")
#adds ellipses to site; but no variables


fviz_pca_biplot(soil_gwc.pca, label = "var", habillage=soil_gwc$fire_freq,
                col.var = "black", repel = TRUE,
                legend.title = "Fire Frequency",
                addEllipses=TRUE, ellipse.level=0.95,palette = "Dark3")
#by fire frequency with variables


fviz_pca_biplot(soil_gwc.pca, label = "var", habillage=soil_gwc$site,
                col.var = "black", repel = TRUE,
                legend.title = "Site",
                addEllipses=TRUE, ellipse.level=0.95,palette = "Dark3")
#works for site with variables

#Linear model section

#1: gravimetric water content (gwc)-----------------
#Prelim visualization
hist(soil_gwc$gwc_percent)
#vaguely normal
boxplot(soil_gwc$gwc_percent~soil_gwc$fire_freq*soil_gwc$site)

#lm model

#interaction
GWC <-lm(gwc_percent~fire_freq*site,data=soil_gwc)
check_model(GWC)
#normal-ish

#Check output
summary(GWC)

#Run Anova on it;
Anova(GWC) 
#site is significant

#Pairwise comparisions
emmeans(GWC,pairwise~site, adjust="Tukey")
#DL-FOJ sign. diff (slightly); FOJ-TOA sign. diff

#Results figure
ggerrorplot(data=soil_gwc, x="site", y="gwc_percent")

#1b: by site alone#
GWC1 <-lm(gwc_percent~site,data=soil_gwc)
check_model(GWC1)
#normal-ish

#Check output
summary(GWC1)

#Run Anova on it;
Anova(GWC1) 
#site is significant; like above

#Pairwise comparisions
emmeans(GWC1,pairwise~site, adjust="Tukey")
#no diffs in the interaction

#1c. by fire freq. only
GWC2 <-lm(gwc_percent~fire_freq,data=soil_gwc)
Anova(GWC2) 
#not sign.

#2: Organic matter (OM)-----------------
#Prelim visualization
hist(soil_gwc$OM)
#vaguely normal
boxplot(soil_gwc$OM~soil_gwc$fire_freq*soil_gwc$site)

#lm model
OM1 <-lm(OM~fire_freq*site,data=soil_gwc)
check_model(OM1)
#normal-ish

#Check output
summary(OM1)

#Run Anova on it;
Anova(OM1) 
#site is sign. diff.

#Pairwise comparisions
emmeans(OM1,pairwise~site, adjust="Tukey")
#DL-FOJ sign diff; FOJ - TOA sign. diff

#Results figure
ggerrorplot(data=soil_gwc, x="site", y="OM")

#3: P1 (_______)-----------------
#Prelim visualization
hist(soil_gwc$P1)
#not normal; right skewed
boxplot(soil_gwc$P1~soil_gwc$fire_freq*soil_gwc$site)

#lm model
P1 <-lm(P1~fire_freq*site,data=soil_gwc)
check_model(P1)
#normal-ish

#Check output
summary(P1)

#Run Anova on it;
Anova(P1) 
#Fire freq is sign; site is sign.

#Pairwise comparisions
####site###
emmeans(P1,pairwise~site, adjust="Tukey")
#DL-FOJ only sign diff
####fire freq###
emmeans(P1,pairwise~fire_freq, adjust="Tukey")
#annual-control sig.

emmeans(P1,pairwise~site*fire_freq, adjust="Tukey")

#Results figure
##fire freq##
ggerrorplot(data=soil_gwc, x="fire_freq",y="P1")
#sign. diff

ggerrorplot(data = soil_gwc, x = "site", y = "P1", 
            color = "fire_freq", palette = c("#40B0A6", "#E1BE6A"),
            legend.title = "Fire Frequency",
            xlab = "Fire Freq", ylab = "P")

#4: P2 (______) -----------------
#Prelim visualization
hist(soil_gwc$P2)
#Left skewed
boxplot(soil_gwc$P2~soil_gwc$fire_freq*soil_gwc$site)

#lm model
#interaction
P2 <-lm(P2~fire_freq*site,data=soil_gwc)
check_model(P2)
#normal-ish

#Check output
summary(P2)

#Run Anova on it;
Anova(P2) 
#fire freq and site is significant

#Pairwise comparisions
emmeans(P2,pairwise~site, adjust="Tukey")
#no diffs in site
emmeans(P2,pairwise~fire_freq, adjust="Tukey")
#sign. different between fire freq
emmeans(P2,pairwise~site*fire_freq, adjust="Tukey")

#Results figure
ggerrorplot(data=soil_gwc, x="fire_freq", y="P2")

#5: Potassium (K) -----------------
#Prelim visualization
hist(soil_gwc$K)
#Right skewed
boxplot(soil_gwc$P2~soil_gwc$fire_freq*soil_gwc$site)

#lm model
#interaction
K <-lm(K~fire_freq*site,data=soil_gwc)
check_model(K)
#normal-ish

#Check output
summary(K)

#Run Anova on it;
Anova(K) 
#no sign. diff anywhere

#6: Magnesium (Mg)-----------------
#Prelim visualization
hist(soil_gwc$Mg)
#right skewed
boxplot(soil_gwc$Mg~soil_gwc$fire_freq*soil_gwc$site)

#lm model
MG <-lm(Mg~fire_freq*site,data=soil_gwc)
check_model(MG)
#normal-ish

#Check output
summary(MG)

#Run Anova on it;
Anova(MG) 
#site is sign. diff.

#Pairwise comparisions
emmeans(MG,pairwise~site, adjust="Tukey")

#Results figure
ggerrorplot(data=soil_gwc, x="site", y="Mg")

#7: Calcium (Ca)-----------------
#Prelim visualization
hist(soil_gwc$Ca)
#pretty normal
boxplot(soil_gwc$Ca~soil_gwc$fire_freq*soil_gwc$site)

#lm model
CA <-lm(Ca~fire_freq*site,data=soil_gwc)
check_model(CA)
#normal-ish

#Check output
summary(CA)

#Run Anova on it;
Anova(CA) 
#site is sign. diff.

#Pairwise comparisions
emmeans(CA,pairwise~site, adjust="Tukey")
#DL-FJ is sign. diff; FOJ-TOA is very sign. diff!

#Results figure
ggerrorplot(data=soil_gwc, x="site", y="Ca")

#8: pH -----------------
#Prelim visualization
hist(soil_gwc$pH)
#left skewed
boxplot(soil_gwc$pH~soil_gwc$fire_freq*soil_gwc$site)

#log transform
ph2<-lm(log(pH)~fire_freq*site, data=soil_gwc)
hist(ph2$residuals)
qqPlot(ph2)
Anova(ph2)

#recip transform
ph3<-lm(1/(pH)~fire_freq*site, data=soil_gwc)
hist(ph3$residuals)
qqPlot(ph3)
Anova(ph3)

#lm model
PH <-lm(pH~fire_freq*site,data=soil_gwc)
check_model(PH)
#normal-ish

#Check output
summary(PH)

#Run Anova on it;
Anova(PH) 
#site is sign. diff.

#Pairwise comparisions
emmeans(PH,pairwise~site, adjust="Tukey")
#DL-TOA is slightly sign. diff; FOJ-TOA is sign. diff

#Results figure
ggerrorplot(data=soil_gwc, x="site", y="pH")

#9: Cation Exchange Capacity (CEC) -----------------
#Prelim visualization
hist(soil_gwc$CEC)
#slightly normal
boxplot(soil_gwc$CEC~soil_gwc$fire_freq*soil_gwc$site)

#lm model
CEC <-lm(CEC~fire_freq*site,data=soil_gwc)
check_model(CEC)

#Check output
summary(CEC)

#Run Anova on it;
Anova(CEC) 
#site is sign. diff.

#Pairwise comparisions
emmeans(CEC,pairwise~site, adjust="Tukey")
#DL-TOA is  sign. diff; FOJ-TOA is very diff

#Results figure
ggerrorplot(data=soil_gwc, x="site", y="CEC")



#10: Nitrate (NO3) -----------------
#Prelim visualization
hist(soil_gwc$NO3)

#need to transform data
boxplot(soil_gwc$NO3~soil_gwc$fire_freq*soil_gwc$site)

boxplot(soil_gwc$NO3~soil_gwc$fire_freq)

#lm model
NO3 <-lm(NO3~fire_freq*site,data=soil_gwc)
check_model(NO3)
hist(NO3$residuals)

#Check output
summary(NO3)

#Run Anova on it;
Anova(NO3) 
#nothing is sign. diff. here

#log transformation
m2<-lm(log(NO3)~fire_freq*site, data=soil_gwc)
hist(m2$residuals)
qqPlot(m2)
plot(m2)
Anova(m2)

#reciprocal transformation
m4<-lm(1/(NO3)~fire_freq*site, data=soil_gwc)
check_model(m4)
hist(m4$residuals)
qqPlot(m4)
plot(m4)
Anova(m4)

#PUBLICATION FIGURE----------------------------------------

#Fig. 6
library(ggplot2)
library(dplyr)
library(patchwork)

#Define colors
custom_colors <- c("Annual" = "#EC4E3B", "Maintenance" = "#4A91D1")

#plot for P1; Available Phosphorus
plot_p1 <- soil_gwc %>%
  ggplot(aes(x = fire_freq, y = P1)) +
  geom_boxplot(outlier.shape = NA, color = "black", linewidth = 0.5) +
  geom_jitter(aes(fill = fire_freq), width = 0.2, size = 2.5, shape = 21,
              color = "black", stroke = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Fire Regime", y = "Available Phosphorus (ppm)", tag = "A") +
  geom_line(data = tibble(x = c(1, 2), y = c(13, 13)),
            aes(x = x, y = y), linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = tibble(x = 1.5, y = 13.1),
            aes(x = x, y = y, label = "***"), size = 5, inherit.aes = FALSE) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(size = 14, family = "Arial"),
    axis.line = element_line(linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10))

#plot for P2; Total Phosphorus
plot_p2 <- soil_gwc %>%
  ggplot(aes(x = fire_freq, y = P2)) +
  geom_boxplot(outlier.shape = NA, color = "black", linewidth = 0.5) +
  geom_jitter(aes(fill = fire_freq), width = 0.2, size = 2.5, shape = 21,
              color = "black", stroke = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Fire Regime", y = "Total Phosphorus (ppm)", tag = "B") +
  geom_line(data = tibble(x = c(1, 2), y = c(14, 14)),
            aes(x = x, y = y), linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = tibble(x = 1.5, y = 14.1),
            aes(x = x, y = y, label = "***"), size = 5, inherit.aes = FALSE) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(size = 14, family = "Arial"),
    axis.line = element_line(linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10))

combined_plot <- plot_p1 + plot_p2 + 
  plot_layout(ncol = 2) &
  theme(plot.tag = element_text(size = 13, face = "bold", family = "Arial"))

