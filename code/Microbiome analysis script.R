library(phyloseq)
library(qiime2R)#Convert to phyloseq object
library(microbiome)
library(tidyverse)
library(BiocManager)
library(ggh4x) #Linear mixed models
library(ggtext)
library(plotly) #Make ggplot2 graphs interactive
library(ggplot2) #Graphs
library(ggpubr) #Additional ggplot2 themes
library(ANCOMBC) #Differential abundance
library(cowplot) #Combine plots into one image
library(readr)  #Write csv files
library(knitr)  #View tibble in console
library(vegan) #PERMANOVA
library(lme4) #Linear models
library(writexl) #create excel file
library(microshades) #Relative abundance graphs
library(patchwork) #Combine multiple plots
library(gridExtra) #Arrange grid based plots


### Import data into phyloseq object
ps <- qza_to_phyloseq(
  features = "data/SK_RE_merged_ASV_table.qza",
  tree = "data/SK_sepp_tree.qza",
  taxonomy = "data/SK_SILVA_138.qza",
  metadata = "data/SK_metadata.tsv"
)

#Find number of reads per sample
sample_sums(ps)
sort(sample_sums(ps))
mean(sample_sums(ps))

#Add 'sampleid' column to metadata
ps@sam_data$sample_id = rownames(ps@sam_data)

#Filter out Chloroplast and Mitochondria
ps <- subset_taxa(ps, Genus != "Chloroplast" & Order != "Chloroplast")
ps <- subset_taxa(ps, Genus != "Mitochondria" & Order != "Mitochondria")
#Check
table(tax_table(ps)[, "Genus"], exclude = NULL)

# Here we rarefy the table, rngseed = the random seed we use -- useful for reproducbility
ps_rar <- rarefy_even_depth(ps, sample.size = 9510, rngseed = 1337,verbose=TRUE)

## `set.seed(1337)` was used to initialize repeatable random subsampling, record this
## Please record this for your records so others can reproduce.

ps_rar

#Filter samples 
ps_filt <- ps_rar

# Set the relative abundance threshold
threshold <- 0.0005 # = 0.05%

# Calculate the total counts for each sample
total_counts <- colSums(ps_rar@otu_table)

# Calculate the threshold counts for each sample
threshold_counts <- total_counts * threshold

# Multiply the OTU table by a logical matrix indicating which values are above the threshold
filtered <- ps_rar@otu_table * (ps_rar@otu_table >= threshold_counts)

# Load back into out phyloseq object
ps_filt@otu_table <- otu_table(filtered, taxa_are_rows = TRUE)

# Check out how much data was removed:
scales::percent(sample_sums(ps_filt) / sample_sums(ps_rar), accuracy = 0.1)

ps_filt

#####################################Samples from Individuals across Time

#Sub-sampling: isolate the recaptures
ps_resamples <- subset_samples(ps_filt, Recapture == "Y")
alpha_diversityresamp <- alpha(ps_resamples, index = "all")

#Bring in metadata
metadataresamp <- meta(ps_resamples)
metadataresamp$name <- rownames(metadataresamp)
alpha_diversityresamp$name <- rownames(alpha_diversityresamp)
alpha_diversity_metadata_resamp <- merge(alpha_diversityresamp, metadataresamp, by = "name")


                 #Calculate Core ASVS, relative abundance of core ASVs
options(digits=3)

##How many ASVs are always present in a given wombat through time?
#Create a list of each subset (Wombat) 
Microchip_list <- c(
  ps_783872B <- subset_samples(ps_filt, Microchip=="783872B"),
  ps_783DCF3 <- subset_samples(ps_filt, Microchip=="783DCF3"),
  ps_7ABD26D <- subset_samples(ps_filt, Microchip=="7ABD26D"),
  ps_7ABE063 <- subset_samples(ps_filt, Microchip=="7ABE063"),
  ps_7AC5251 <- subset_samples(ps_filt, Microchip=="7AC5251"),
  ps_7ACA357 <- subset_samples(ps_filt, Microchip=="7ACA357"),
  ps_94320324122 <- subset_samples(ps_filt, Microchip=="94320324122")
)

#Create functions to use in loop
gen_x <- function(ps){
  print(length(taxa_names(filter_taxa(ps, function(x) sum(x) > 0 , TRUE))))
}
gen_y <- function(ps){
  print(length(taxa_names(core(ps, detection =1, prevalence =0.999))))
}
gen_z <- function(ps){
  print(mean(sample_sums(core(ps, detection =1, prevalence =0.999))) / mean(sample_sums(ps)))
}
gen_zz <- function(ps){
  print(length(taxa_names(core(ps, detection = 1, prevalence = 0))) - 
          length(taxa_names(core(ps, detection = 1, prevalence = 1/length(sample_sums(ps))))))
}

#Create dataframe of length wombatlist
ASVresultsresamp <- data.frame(totalASVs_per_wombat=rep(0,length(Microchip_list)),
                               total_core_ASVs=rep(0,length(Microchip_list)),
                               mean_rel_abun_core_ASVs=rep(0,length(Microchip_list)),
                               num_singleton_ASVs_per_wombat=rep(0, length(Microchip_list)))
row.names(ASVresultsresamp) <- sort(unique(ps_resamples@sam_data$Microchip))

#Loop over each subsetted phyloseq object, outputting results into the dataframe!
for (ps in 1:length(Microchip_list)){
  x <- gen_x(Microchip_list[[ps]])
  y <- gen_y(Microchip_list[[ps]])
  z <- gen_z(Microchip_list[[ps]])
  zz <- gen_zz(Microchip_list[[ps]])
  ASVresultsresamp[ps, 1] <- x
  ASVresultsresamp[ps, 2] <- y
  ASVresultsresamp[ps, 3] <- z
  ASVresultsresamp[ps, 4] <- zz
}

                                #####Figure 2
#Create theme for PCoA

theme_pcoa <- theme(axis.text.x = element_text(face="bold", size=16), 
                    axis.text.y = element_text(face="bold", size=16),
                    axis.title.x = element_text(size=20, face="bold"),
                    axis.title.y = element_text(size=20, face="bold"),
                    axis.line = element_line(colour = "black"),
                    #Background panel
                    panel.background = element_rect(fill = "White"),
                    panel.grid.major = element_line(colour = "white"), 
                    panel.grid.minor = element_line(colour = "white"),
                    #Legend
                    legend.title = element_blank(),
                    legend.text = element_text(size=16),
                    legend.key = element_rect(fill = "white", color = NA),
                    legend.key.size = unit(2.5, "line"))

#Distance calculations for plotting ordination
#Unweighted UniFrac
unweighted_unifrac <- ordinate(ps_resamples, method = "PCoA", distance = "unifrac", weighted=F)

(ord.unwght.resamp <-plot_ordination(physeq = ps_resamples, 
                                   ordination = unweighted_unifrac, 
                                   color = "Microchip",
                                   axes = c(1, 2)) +
                                   theme_minimal() +
                                   geom_line())

(ord.uw.1.2 <- ord.unwght.resamp +
    scale_x_reverse() +
    geom_point(size=7.5) +
    geom_path() +
    geom_text(aes(label = Sample_Order), color = "black", fontface = "bold") +
    theme(legend.position = "none") +
    ggtitle("Unweighted") +
    theme_pcoa)

#Weighted UniFrac
weighted_unifrac <- ordinate(ps_resamples, method = "PCoA", distance = "unifrac", weighted=T)

(ord.wght.resamp <-plot_ordination(physeq = ps_resamples, 
                ordination = weighted_unifrac, 
                color = "Microchip",
                axes = c(1, 2)) +
  theme_minimal() +
  geom_line())

(ord.w.1.2 <- ord.wght.resamp +
  scale_x_reverse() +
  geom_point(size=7.5) +
  geom_path() +
  geom_text(aes(label = Sample_Order), color = "black", fontface = "bold") +
    scale_fill_discrete(labels = c("W1", "W2", "W3","W4", "W5", "W6", "W7")) +
  theme(legend.position = "right") +
    ggtitle("Weighted")+
    theme_pcoa)

#Combine them
(resamp.ord <- plot_grid(ord.uw.1.2, ord.w.1.2 + theme(legend.position="none")))

# create some space to the left of the legend
legend <- get_legend(
  ord.uw.1.2 + theme(legend.box.margin = margin(0, 0, 0, 12)))

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
plot_grid(resamp.ord, legend, rel_widths = c(3, .4))

                  #PERMANOVA Individual Identity and Season

metadata_wombat_resamp <- as(sample_data(ps_resamples), "data.frame")

adonis2(distance(ps_resamples, method="unweighted_unifrac") ~
         Season, 
        data = metadata_wombat_resamp, 
        permutations = 100000)

adonis2(distance(ps_resamples, method="weighted_unifrac") ~ 
          Season, 
        data = metadata_wombat_resamp, 
        permutations = 100000)

                              #Alpha Diversity Individual Identity

#Wilcoxon rank-sum test (2 means)
wilcox.test(evenness_simpson~ Season, data = alpha_diversity_metadata_resamp)
#Also used for variables 'Sex' and 'Age_class'
#Change 'evenness_simpson' to 'evenness_pielou', 'observed', 'diversity_shannon', 'diversity_fisher'

#Kruskal-Wallis test (Typically 3 or more means)
kruskal.test(evenness_simpson~ Microchip, data = alpha_diversity_metadata_resamp)
#Change 'evenness_simpson' to 'evenness_pielou', 'observed', 'diversity_shannon', 'diversity_fisher'


                               #Dominant Phyla
#View Phylum present in individual samples
data_resamp <- ps_resamples %>%
  tax_glom("Phylum") %>%
  phyloseq::transform_sample_counts(function(x) { x/sum(x) }) %>%
  psmelt() %>%
  as_tibble

                               #Supplementary Figure 1
#Microshades relative abundance figure

mdf_prep <- ps_resamples %>%
         tax_glom("Genus") %>%
        phyloseq::transform_sample_counts(function(x) { x/sum(x) }) %>%
         psmelt() %>%
        filter(Abundance > 0)

color_objs_resamp <- create_color_dfs(mdf_prep,selected_groups = c('Proteobacteria', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes') , cvd = TRUE)

mdf_resmp <- color_objs_resamp$mdf
cdf_resamp <- color_objs_resamp$cdf

plot <- plot_microshades(mdf_resamp, cdf_resamp)

resamp_legend <-custom_legend(mdf_resamp, cdf_resamp)

plot_diff <- plot + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(size= 6)) +
  facet_wrap(~Microchip, scales = "free_x", nrow = 2) +
  theme(axis.text.x = element_text(size= 6)) + 
  theme(plot.margin = margin(6,20,6,6)) 

plot_grid(plot_diff, resamp_legend,  rel_widths = c(1, .25))

#Note sample ID labels modified outside of R

                              #Remove samples from 783872B - outlier

ps_resamples2 <- subset_samples(ps_resamples, Microchip %in% c("783872B", "783DCF3", "7ABD26D", "7ABE063", "7AC5251", "94320324122"))
alpha_diversityresamp2 <- alpha(ps_resamples2, index = "all")

#Bring in metadata
metadataresamp2 <- meta(ps_resamples2)
metadataresamp2$name <- rownames(metadataresamp2)
alpha_diversityresamp2$name <- rownames(alpha_diversityresamp2)
alpha_diversity_metadataresamp2 <- merge(alpha_diversityresamp2, metadataresamp2, by = "name")

#Repeat above code for Beta and Alpha diversity with 783872B removed


                                    #LOCATION SAMPLES #  

#Isolate samples from each location: exclude unwanted recapture samples
ps_seasonal <- subset_samples(ps_filt, Seasonal == "Y")
alpha_diversity_seasonal <- alpha(ps_seasonal, index = "all")

metadata_seasonal <- meta(ps_seasonal)
metadata_seasonal$name <- rownames(metadata_seasonal)
alpha_diversity_seasonal$name <- rownames(alpha_diversity_seasonal)
alpha_diversity_metadata_seasonal <- merge(alpha_diversity_seasonal, metadata_seasonal, by = "name")

#Find Core ASVs, relative abundance of Core ASVs etc. for each location
options(digits=3)

##How many ASVs are always present in a given location through time?
#Create a list of each subset (Location)
Location_list <- c(
  ps_Kooloola <- subset_samples(ps_filt, Location=="Kooloola"),
  ps_Wonga <- subset_samples(ps_filt, Location=="Wonga"),
  ps_BonBon <- subset_samples(ps_filt, Location=="Bon Bon"),
  ps_EyreP <- subset_samples(ps_filt, Location=="Eyre Peninsula"), 
)

#Create functions to use in loop
gen_x <- function(ps){
  print(length(taxa_names(filter_taxa(ps, function(x) sum(x) > 0 , TRUE))))
}
gen_y <- function(ps){
  print(length(taxa_names(core(ps, detection =1, prevalence =0.999))))
}
gen_z <- function(ps){
  print(mean(sample_sums(core(ps, detection =1, prevalence =0.999))) / mean(sample_sums(ps)))
}
gen_zz <- function(ps){
  print(length(taxa_names(core(ps, detection = 1, prevalence = 0))) - 
          length(taxa_names(core(ps, detection = 1, prevalence = 1/length(sample_sums(ps))))))
}

#Create dataframe of length Location list
ASVresults <- data.frame(totalASVs_per_location=rep(0,length(Location_list)),
                         total_core_ASVs=rep(0,length(Location_list)),
                         mean_rel_abun_core_ASVs=rep(0,length(Location_list)),
                         num_singleton_ASVs_per_location=rep(0, length(Location_list)))
row.names(ASVresults) <- sort(unique(ps_filt@sam_data$Location))

#Loop over each subsetted phyloseq object, outputting results into the dataframe!
for (ps in 1:length(Location_list)){
  x <- gen_x(Location_list[[ps]])
  y <- gen_y(Location_list[[ps]])
  z <- gen_z(Location_list[[ps]])
  zz <- gen_zz(Location_list[[ps]])
  ASVresults[ps, 1] <- x
  ASVresults[ps, 2] <- y
  ASVresults[ps, 3] <- z
  ASVresults[ps, 4] <- zz
}


                            #Figure 3
#PCoA
#Create PCoA Theme
theme_pcoa <- theme(axis.text.x = element_text(face="bold", size=16), 
                    axis.text.y = element_text(face="bold", size=16),
                    axis.title.x = element_text(size=20, face="bold"),
                    axis.title.y = element_text(size=20, face="bold"),
                    axis.line = element_line(colour = "black"),
                    #Background panel
                    panel.background = element_rect(fill = "White"),
                    panel.grid.major = element_line(colour = "white"), 
                    panel.grid.minor = element_line(colour = "white"),
                    #Legend
                    legend.title = element_blank(),
                    legend.text = element_text(size=16),
                    legend.key = element_rect(fill = "white", color = NA),
                    legend.key.size = unit(2.5, "line"))


#PCOA ordination of unweighted UniFrac distances 
unweighted_unifrac <- ordinate(ps_seasonal, method = "PCoA", distance = "unifrac", weighted=F)

weighted_unifrac <- ordinate(ps_seasonal, method = "PCoA", distance = "unifrac", weighted=T)

#Axes 1/2
Location.uw <-plot_ordination(physeq = ps_seasonal, 
                              ordination = unweighted_unifrac, 
                              color = "Location",
                              axes = c(1, 2)) +
  theme_pcoa +
  theme(legend.position = "top") + 
  geom_point(size = 3, alpha = 0.6)

Location.w <-plot_ordination(physeq = ps_seasonal, 
                             ordination = weighted_unifrac, 
                             color = "Location",
                             axes = c(1, 2)) +
  theme_pcoa +
  geom_point(size = 3, alpha = 0.6)

#Axes 2/3
Location.uw <-plot_ordination(physeq = ps_seasonal, 
                              ordination = unweighted_unifrac, 
                              color = "Location",
                              axes = c(2, 3)) +
  theme_pcoa +
  geom_point(size = 3, alpha = 0.6)

Location.w <-plot_ordination(physeq = ps_seasonal, 
                             ordination = weighted_unifrac, 
                             color = "Location",
                             axes = c(2, 3)) +
  theme_pcoa +
  geom_point(size = 3, alpha = 0.6)

#Create PCoAs for season within each location

# WONGA

ps_locationW <- subset_samples(ps_seasonal, Location == "Wonga")
alpha_diversity <- alpha(ps_locationW, index = "all")

metadata <- meta(ps_filt)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadataW <- merge(alpha_diversity, metadata, by = "name")

#Distance Calculation for PCoA
unweighted_unifrac <- ordinate(ps_locationW, method = "PCoA", distance = "unifrac", weighted=F)

#Axes 1/2
(W <- plot_ordination(physeq = ps_locationW, 
                      ordination = unweighted_unifrac, 
                      title = "Wonga" ,
                      color = "Season",
                      axes = c(1, 2)) +
    theme_minimal() +
    theme_pcoa +
    geom_point(size = 3, alpha = 1, shape = 8) +
    scale_color_manual(values=c("chocolate2", "skyblue1")))

# KOOLOOLA
ps_locationK <- subset_samples(ps_seasonal, Location =="Kooloola")
alpha_diversityK <- alpha(ps_locationK, index = "all")

metadata <- meta(ps_seasonal)
metadata$name <- rownames(metadata)
alpha_diversityK$name <- rownames(alpha_diversityK)
alpha_diversity_metadataK <- merge(alpha_diversityK, metadata, by = "name")

#Distance calculation for PCoA
unweighted_unifrac <- ordinate(ps_locationK, method = "PCoA", distance = "unifrac", weighted=F)

#PCoA
#Axes 1/2
(K <- plot_ordination(physeq = ps_locationK, 
                      ordination = unweighted_unifrac, 
                      title = "Kooloola" ,
                      color = "Season",
                      axes = c(1, 3)) +
    theme_minimal() +
    theme_pcoa +
    geom_point(size = 3, alpha = 1, shape = 7) +
    scale_color_manual(values=c("chocolate2", "skyblue1")))

# BON BON

ps_locationB <- subset_samples(ps_seasonal, Location == "Bon Bon")
alpha_diversity <- alpha(ps_locationB, index = "all")

metadata <- meta(ps_filt)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadataB <- merge(alpha_diversity, metadata, by = "name")

#PCoA
#Axes 1/2
(B <- plot_ordination(physeq = ps_locationB, 
                      ordination = unweighted_unifrac, 
                      title = "Bon Bon" ,
                      color = "Season",
                      shape = "Site" ,
                      axes = c(1, 2)) +
    theme_minimal() +
    theme_pcoa +
    geom_point(size = 3, alpha = 1) +
    scale_color_manual(values=c("chocolate2", "skyblue1")))

# EYRE PENINSULA

ps_locationE <- subset_samples(ps_filt, Location == "Eyre Peninsula")
alpha_diversity <- alpha(ps_locationE, index = "all")

metadata <- meta(ps_filt)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadataE <- merge(alpha_diversity, metadata, by = "name")

#Distance Calculation for PCoA
unweighted_unifrac <- ordinate(ps_locationE, method = "PCoA", distance = "unifrac", weighted=F)

#PCoA
#Axes 1/2
(E <- plot_ordination(physeq = ps_locationE, 
                      ordination = unweighted_unifrac, 
                      title = "Eyre Peninsula" ,
                      color = "Season",
                      shape = "Site",
                      axes = c(1, 2)) +
    theme_minimal() +
    theme_pcoa +
    geom_point(size = 3, alpha = 1) +
    scale_shape_manual(values=c(15, 3)) +
    scale_color_manual(values=c("chocolate2", "skyblue1")))

# Combine main location PCoA with individual location PCoAs
main <- plot_grid(K + theme(legend.position = "none"), W + theme(legend.position = "none"), B + theme(legend.position = "none"), E + theme(legend.position = "none"), 
                  ncol = 2, nrow = 2)
main
legend_b <- get_legend(K + theme(legend.margin = margin(t =0, unit ='cm')))
(main1 <- plot_grid(main, get_legend(K), ncol=2 + theme(legend.position = "bottom", legend.margin = margin(t=0, unit = 'cm'))))
(plot_grid(main, legend_b))
(main2 <- plot_grid(Location, main1, ncol = 1, nrow = 2))

#Location PERMANOVA

metadata_wombat <- as(sample_data(ps_seasonal), "data.frame")

adonis2(distance(ps_seasonal, method="weighted_unifrac") ~ 
          Location, 
        data = metadata_wombat, 
        permutations = 100000)

adonis2(distance(ps_seasonal, method="unweighted_unifrac") ~ 
          Location, 
        data = metadata_wombat, 
        permutations = 100000)

#Repeat PERMANOVA as above for variable 'Season'


                                 
                                  #### Figure 4
#Alpha Diversity of Location and Season

#Observed ASVs boxplot
(Obs.bp <- ggplot(alpha_diversity_metadata_seasonal, aes(x = Site, y = observed, colour = Season)) +
geom_boxplot() +
    theme_classic() +
    ylab("Observed ASVs") +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size=9)) +
        scale_colour_manual(values=c("chocolate2", "skyblue1")))
#Simpson Evenness boxplot
(Evs.bp <- ggplot(alpha_diversity_metadata_seasonal, aes(x = Site, y = evenness_simpson, colour = Season)) +
    geom_boxplot() +
    theme_classic() +
    ylab("Simpson evenness") +
    scale_colour_manual(values=c("chocolate2", "skyblue1")) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size=9)))
#Shannon Diversity boxplot
(Dis.bp <- ggplot(alpha_diversity_metadata_seasonal, aes(x = Site, y = diversity_shannon, colour = Season)) +
    geom_boxplot() +
    theme_classic() +
    ylab("Shannon diversity") +
    scale_colour_manual(values=c("chocolate2", "skyblue1")) +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size=9)))

#Combine all boxplots
(Comb.alpha.bp <- plot_grid(Obs.bp + theme(legend.position = "none"), Evs.bp + theme(legend.position = "none"), Dis.bp + theme(legend.position = "none"), ncol = 1))


#Alpha Diversity Tests

#Repeat all tests with alpha diversity measures 'observed', 'evenness_pielou', 'diversity_shannon', 'diversity_fisher'
#Determine alpha diversity differences between sites, between seasons, and between sites and seasons

#ANOVA and Tukeys Post Hoc Test 
ANOVA1 <-aov(evenness_simpson ~ Site*Season, data = alpha_diversity_metadata_seasonal)
summary(ANOVA1)
TukeyHSD(ANOVA1)



                                           ####Figure 5
#Microshades Relative Abundance Figure
mdf_prep1 <- ps_seasonal %>%
  tax_glom("Genus") %>%
  phyloseq::transform_sample_counts(function(x) { x/sum(x) }) %>%
  psmelt() %>%
  filter(Abundance > 0)

color_objs_snl <- create_color_dfs(mdf_prep1,selected_groups = c('Proteobacteria', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes') , cvd = TRUE)

mdf_snl <- color_objs_snl$mdf
cdf_snl <- color_objs_snl$cdf

plot <- plot_microshades(mdf_snl, cdf_snl)

GP_legend1 <-custom_legend(mdf_snl, cdf_snl)

plot_diff <- plot + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(size= 6)) +
  facet_wrap(~Season + Location, scales = "free_x", nrow = 2) +
  theme(axis.text.x = element_blank()) + 
  theme(plot.margin = margin(6,20,6,6)) +
  theme (strip.text.x = element_text(size = 10))

plot_grid(plot_diff, GP_legend1,  rel_widths = c(1, .25))



                    #Seasonal Variation within each Location
                    # WONGA

ps_locationW <- subset_samples(ps_seasonal, Location == "Wonga")
alpha_diversity <- alpha(ps_locationW, index = "all")

metadata <- meta(ps_filt)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadataW <- merge(alpha_diversity, metadata, by = "name")

#Seasonal Alpha Diversity comparison for Wonga samples
#Wilcoxon rank-sum test 
wilcox.test(evenness_simpson~ Season, data = alpha_diversity_metadataW)
#Repeat with 'observed', 'evenness_pielou', 'diversity_shannon', and 'diversity_fisher'

#Beta Diversity PERMANOVA
metadataW <- as(sample_data(ps_locationW), "data.frame")

adonis2(distance(ps_locationW, method="unifrac") ~ Season, data = metadataW, permutations = 100000)

adonis2(distance(ps_locationW, method="weighted_unifrac") ~ Season, data = metadataW, permutations = 100000)

#View Phylum present in Wonga Samples
data_seasonW <- ps_locationW %>%
  tax_glom("Phylum") %>%
  phyloseq::transform_sample_counts(function(x) { x/sum(x) }) %>%
  psmelt() %>%
  as_tibble

                                   # KOOLOOLA

ps_locationK <- subset_samples(ps_seasonal, Location =="Kooloola")
alpha_diversityK <- alpha(ps_locationK, index = "all")

metadata <- meta(ps_seasonal)
metadata$name <- rownames(metadata)
alpha_diversityK$name <- rownames(alpha_diversityK)
alpha_diversity_metadataK <- merge(alpha_diversityK, metadata, by = "name")

#Seasonal Alpha Diversity comparison for Kooloola samples
#Wilcoxon rank-sum test 
wilcox.test(evenness_simpson~ Season, data = alpha_diversity_metadataK)
#Repeat with 'observed', 'evenness_pielou', 'diversity_shannon', and 'diversity_fisher'

#Beta Diversity PERMANOVA
metadataK <- as(sample_data(ps_locationK), "data.frame")

adonis2(distance(ps_locationK, method="unifrac") ~ Season, data = metadataK, permutations = 100000)

adonis2(distance(ps_locationK, method="weighted_unifrac") ~ Season, data = metadataK, permutations = 100000)

#View Phylum present in Kooloola samples
data_seasonK <- ps_locationK %>%
  tax_glom("Phylum") %>%
  phyloseq::transform_sample_counts(function(x) { x/sum(x) }) %>%
  psmelt() %>%
  as_tibble

                                    # BON BON

ps_locationB <- subset_samples(ps_seasonal, Location == "Bon Bon")
alpha_diversity <- alpha(ps_locationB, index = "all")

metadata <- meta(ps_filt)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadataB <- merge(alpha_diversity, metadata, by = "name")

#Seasonal Alpha Diversity comparison for Bon Bon samples
#ANOVA and Tukeys Post Hoc
ANOVA_B <-aov(evenness_simpson~Site*Season, data = alpha_diversity_metadataB)
TukeyHSD(ANOVA_B)
#Repeat with 'observed', 'evenness_pielou', 'diversity_shannon', and 'diversity_fisher'

#Diversity Calculation  for PCoA
unweighted_unifrac <- ordinate(ps_locationB, method = "PCoA", distance = "unifrac", weighted=F)


#Beta Diversity PERMANOVA
metadataB <- as(sample_data(ps_locationB), "data.frame")

adonis2(distance(ps_locationB, method="unifrac") ~ Season, data = metadataB, permutations = 100000)

adonis2(distance(ps_locationB, method="weighted_unifrac") ~ Season, data = metadataB, permutations = 100000)

#View Phylum present in Bon Bon Samples
data_seasonB <- ps_locationB %>%
  tax_glom("Phylum") %>%
  phyloseq::transform_sample_counts(function(x) { x/sum(x) }) %>%
  psmelt() %>%
  as_tibble

                                       # EYRE PENINSULA

ps_locationE <- subset_samples(ps_filt, Location == "Eyre Peninsula")
alpha_diversity <- alpha(ps_locationE, index = "all")

metadata <- meta(ps_filt)
metadata$name <- rownames(metadata)
alpha_diversity$name <- rownames(alpha_diversity)
alpha_diversity_metadataE <- merge(alpha_diversity, metadata, by = "name")

#Seasonal Alpha Diversity comparison for Eyre Peninsula samples
#ANOVA and Tukeys Post Hoc
ANOVA_E <-aov(evenness_simpson~Site*Season, data = alpha_diversity_metadataE)
TukeyHSD(ANOVA_E)
#Repeat with 'observed', 'evenness_pielou', 'diversity_shannon', and 'diversity_fisher'


#Beta Diversity PERMANOVA
metadataE <- as(sample_data(ps_locationE), "data.frame")

adonis2(distance(ps_locationE, method="unifrac") ~ Season, data = metadataE, permutations = 100000)
adonis2(distance(ps_locationE, method="weighted_unifrac") ~ Season, data = metadataE, permutations = 100000)

#View Phylum present in Eyre Peninsula samples
data_seasonE <- ps_locationE %>%
  tax_glom("Phylum") %>%
  phyloseq::transform_sample_counts(function(x) { x/sum(x) }) %>%
  psmelt() %>%
  as_tibble


##########################################################################################
### ANCOMBC-2 analyses
### Raphael Eisenhofer Jan 2024
##########################################################################################

## ANCOM-BC (Raph)
library(ANCOMBC)
library(knitr)

#Run ancombc2 for wombat location and season
ancombc2_location_season <- ancombc2(data = ps_seasonal, 
                                     assay_name = "counts", 
                                     tax_level = "Genus", #change to agglomerate analysis to a higher taxonomic range
                                     fix_formula = "Location + Season", #fixed variable(s)
                                     lib_cut = 8600, 
                                     group = "Location",
                                     pairwise = TRUE,
                                     n_cl = 5)

#View results for paired locations
ancom_location_season <- ancombc2_location_season$res_pair
write_xlsx(ancom_location_season,"/ancom_location_season_paired.xlsx")

#View results for season
ancom_location_season <- ancombc2_location_season$res
write_xlsx(ancom_location_season,"/ancom_location_season.xlsx")

#Supplementary Information Figure 3
#ANCOMBC Season Figure

#from https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html

res <- ancombc2_location_season$res

df_fig_pair1 = res %>%
  dplyr::filter(diff_SeasonWinter_Spring == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_SeasonWinter_Spring == 1, round(lfc_SeasonWinter_Spring, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_pair2 = res %>%
  dplyr::filter(diff_SeasonWinter_Spring == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_SeasonWinter_Spring == 1, 
                              "aquamarine3")) %>%
  tidyr::pivot_longer(cols = lfc1, 
                      names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_pair = df_fig_pair1 %>%
  dplyr::left_join(df_fig_pair2, by = c("taxon", "group"))

df_fig_pair$group = recode(df_fig_pair$group, 
                           `lfc1` = "Summer_Autumn - Winter_Spring",)

df_fig_pair$group = factor(df_fig_pair$group, 
                           levels = c("Summer_Autumn - Winter_Spring"))

lo = floor(min(df_fig_pair$value))
up = ceiling(max(df_fig_pair$value))
mid = (lo + up)/2
fig_pair = df_fig_pair %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = "black"), size = 3) +
  scale_color_identity(guide = FALSE) +
  facet.by(~Location)+
  labs(x = NULL, y = NULL, title = "Log fold changes by season") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.margin = margin(10, 10, 10, 30))
fig_pair

