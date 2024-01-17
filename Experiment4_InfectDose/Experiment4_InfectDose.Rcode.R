# E5 analy sis 
library(tidyverse)
library(phyloseq)
library(broom)
library(ggrepel)
library(ggplot2)
library(microbiome)
library(qiime2R)
library(xlsx)
library(RColorBrewer)
library(vegan)
library(DESeq2)
library(ANCOMBC) 
library(Maaslin2)
library(mixOmics)
library(metagMisc)
library(eulerr)

setwd("/Users/liztown/Documents/KalanLab/Papers/Vaginal Microbiome/Rcode_forManuscript/Experiment4_InfectDose")

ps.E5infectDose <- qza_to_phyloseq(features = "./table-E5-dada2.qza",
                                   taxonomy = "./taxonomy-E5.qza",
                                   metadata = "./E5manufest.txt")

E5.meta <- data.frame(sample_data(ps.E5infectDose))
E5.meta$SampleID <- row.names(E5.meta)

ps.Decontam <- subset_taxa(ps.E5infectDose, Family !="mitochondria")
ps.Decontam <- subset_taxa(ps.Decontam, Class !="Chloroplast")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !="Cyanobacteria")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !=" ")
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylobacterium") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingomonas") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ralstonia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Facklamia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Chryseobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curvibacter") 
ps.E5 <- subset_samples(ps.Decontam, !sample_names(ps.Decontam) %in% c("LK16S009-371", "LK16S009-372")) # removing negative controls 
ps.E5 <- phyloseq_filter_prevalence(ps.E5, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR") 
OutCV <- psmelt(ps.E5)
# write.csv(OutCV, "./E5_TotalAbundanceOut.csv")
ps.E5REl <- microbiome::transform(ps.E5, "compositional")
E5REl <- psmelt(ps.E5REl)
#write.csv(E5REl, "./E5_RELATIVEAbundanceOut.csv")
E5Rel.df<-as.data.frame(E5REl)

# Figure 3B - Experiment 4 Relative abundance 
OtherE5 <- E5Rel.df # the dataframe of the relative abundance within each sample (from above)
OtherE5 <- subset.data.frame(OtherE5, Abundance > 0.001)
OtherE5$Genus[OtherE5$Genus == '[Prevotella]'] <- 'Prevotella'
OtherE5$Genus[OtherE5$Genus == 'Propionibacterium'] <- 'Cutibacterium'
OtherE5 <- OtherE5 %>% mutate(Genus = ifelse(Abundance < 0.05 , "Other", Genus)) # replacing the genus catagory with "Other" if the abundance is less than 1% 
OtherE5 <- OtherE5 %>% mutate(Phylum = ifelse(Abundance < 0.05 , "Other", Phylum)) 

OtherE5$ViralInfectGroup <- factor(OtherE5$ViralInfectGroup, levels = c("Mock", "10^4","10^6", "10^8"))
OtherE5$Timepoint <- factor(OtherE5$Timepoint, levels = c("T1", "T2"))

OtherE5$Phylum <- factor(OtherE5$Phylum, levels = c("Other", "Actinobacteria",  "Firmicutes", "Bacteroidetes","Fusobacteria","[Thermi]",  "Proteobacteria"))
OtherE5$Genus<- factor(OtherE5$Genus, levels = c("Other",
                                                 "Arthrobacter","Bifidobacterium","Corynebacterium",'Cutibacterium',"Gordonia","Micrococcus","Nesterenkonia","Nocardioides",
                                                 "Propionibacterium","Streptomyces", ##Actinobacteria / actinomycetota (9)
                                                 "Anoxybacillus","Bacillus", "Brochothrix","Clostridium","Enterococcus", "Lactobacillus","Lactococcus","Macrococcus",
                                                 "Oscillospira", "Ruminococcus", "Staphylococcus", "Streptococcus", "Veillonella",# Firmicutes/ Bacillota (13)
                                                 "Cloacibacterium","[Prevotella]","Porphyromonas","Prevotella","Sediminibacterium", # Bacteroidetes / Bacteroidota (5)
                                                 "Fusobacterium","Leptotrichia", # fusobacteria
                                                 "Deinococcus", # Thermi
                                                 "Acinetobacter","Actinobacillus", "Afipia", "Agrobacterium","Brevundimonas","Caulobacter", "Enhydrobacter","Herbaspirillum",
                                                 "Massilia", "Pseudomonas","Psychromonas", "Schlegelella", "Stenotrophomonas","Sutterella","Zoogloea")) #Proteobacteria /pseudomonadota (15)

BigMicro2 <-c("#C1C1C1", #Gray
              colorRampPalette(c("#031A34", "#053061", "#21777C", "#2A9D8F", "#8AB17D", "#C0C27B","#D2C77A", "#F5D279","#EFB366","#F4A261","#DE7350","#C8443F","#913C44", "#733036"))(44))

ggplot(OtherE5, aes(x= Timepoint, y =Abundance, fill = Genus))+
  geom_bar(stat="identity", position="fill") +
  facet_wrap("ViralInfectGroup", scales = "free_x", nrow = 1)+
  ylab("Relative Abundance")+
  theme_light()+
  scale_fill_manual(values = BigMicro2)+
  scale_color_manual(values = BigMicro2)+
  theme(axis.text.x = element_text(size = 4)) +
  ggtitle("Experiment 4 Genus present > 0.5% of reads in a sample")

# Alpha Abundance (Not shown in paper)
ps.E5@sam_data$ViralInfectGroup <- factor(ps.E5@sam_data$ViralInfectGroup, levels = c("Mock", "10^4","10^6", "10^8"))
ps.E5@sam_data$Timepoint <- factor(ps.E5@sam_data$Timepoint, levels = c("T1", "T2"))
plot_richness(ps.E5, x = "Timepoint", measures=c("Shannon"), color = "ViralInfectGroup") + 
  geom_point(size = 4)+
  #stat_summary(fun.y = median,  geom = "point", aes(group = MouseSampleName),shape = 95, size = 12, show.legend = F)+ 
  facet_wrap("ViralInfectGroup", scales = "free_x", ncol = 4) + 
  theme_light()

tab <-microbiome::alpha(ps.E5, index = c("observed", "chao1", "diversity_inverse_simpson", "diversity_gini_simpson","diversity_shannon",	"diversity_coverage", "evenness_camargo",	"evenness_pielou",
                                                     "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp", "dominance_dmn", "dominance_absolute", "dominance_relative", "dominance_simpson", 
                                                     "dominance_core_abundance", "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_rare_abundance"))

AlphaE <- tab
AlphaE$SampleID<- row.names(AlphaE)
AlphaE5 <- merge(AlphaE, E5.meta, by = "SampleID", all = F) # the alpha table above but with the metadata added back in 
AlphaE5$ViralInfectGroup <- factor(AlphaE5$ViralInfectGroup, levels = c("Mock", "10^4","10^6", "10^8"))
AlphaE5$Timepoint <- factor(AlphaE5$Timepoint, levels = c("T1", "T2"))

ggplot(AlphaE5, aes(x = Timepoint, y =diversity_shannon, color = ViralInfectGroup))+
  geom_point(size = 5)+
  stat_summary(fun.y = mean,  geom = "point",shape = 95, size = 13, show.legend = F)+ 
  scale_color_manual(values = c("#21777C","#8AB17D","#EFB366","#DE7350"))+
  facet_wrap("ViralInfectGroup", scales = "free_x", ncol = 4) + 
  theme_light()+
  ggtitle("Shannon Alpha Diversity Experiment 4")

# Figure 3 C Beta Diversity
min(sample_sums(ps.E5))# minimum sample read is 69
median(sample_sums(ps.E5infectDose)) # 7921.5
median(sample_sums(ps.E5)) # 1833
max(sample_sums(ps.E5)) #83188
rarecurve(t(otu_table(ps.E5)), step=100, ylim =c(0,50), xlim=c(0,600))  ## considering using a read cut off of 5000 for beta diversity metrics 

ps.rarefiedE5 = rarefy_even_depth(ps.E5, rngseed=1, sample.size=300, replace=F) # THIS IS THE ONE TO GO WITH

GP = ps.rarefiedE5
GP.ord <- ordinate(GP, "NMDS",  "bray") # For Bray Curtis of the rareified dataset
plot_ordination(GP, GP.ord, type="samples", color="ViralInfectGroup")+# shape = "EdgeCenter") + 
  geom_point(size=3) + ggtitle("Bray Curtis")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(aes(shape = Timepoint, fill = ViralInfectGroup), size = 6)+
  scale_shape_manual(values = c(19,17))+
  scale_color_manual(values = c("#21777C","#8AB17D","#EFB366","#DE7350"))+
  theme_bw()
# additional Beta diversity plots (not in paper)
plot_ordination(GP, GP.ord, type="samples", color="ViralCopyNumberRanked")+# shape = "EdgeCenter") + 
  geom_point(size=3) + ggtitle("Bray Curtis")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(aes(shape = Timepoint, fill = ViralCopyNumberRanked), size = 6)+
  scale_shape_manual(values = c(19,17))+
  scale_color_manual(values = c( "#053061","#21777C","#CCCCCC"))+
  theme_bw()
plot_ordination(GP, GP.ord, type="samples", color="Cleared.Persist.Endpoint")+# shape = "EdgeCenter") + 
  geom_point(size=3) + ggtitle("Bray Curtis")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(aes(shape = ViralInfectGroup, fill = Cleared.Persist.Endpoint), size = 6)+
  scale_shape_manual(values = c(19,17,15,18))+
  scale_color_manual(values = c("#8AB17D", "#CCCCCC", "#21777C"))+
  theme_bw()

## Table 2: Type-2 PERMANOVAs 
sampledf <- data.frame(sample_data(ps.rarefiedE5)) 
E5_bray <- phyloseq::distance(ps.rarefiedE5, method = "bray")
Bray_clust <-hclust(E5_bray)
adonis2(E5_bray ~ ViralInfectGroup, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(E5_bray ~ ViralCopyNumberRanked, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(E5_bray ~ Timepoint, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(E5_bray ~ Cleared.Persist.Endpoint, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(E5_bray ~ DiseaseCINEndpoint, by= "margin", data = sampledf, permutations = 9999) # 

ps.rareE5.INFECTED <- subset_samples(ps.rarefiedE5, ViralInfectGroup != "Mock")
sampledf <- data.frame(sample_data(ps.rareE5.INFECTED)) 
E5_bray <- phyloseq::distance(ps.rareE5.INFECTED, method = "bray")
Bray_clust <-hclust(E5_bray)
adonis2(E5_bray ~ ViralInfectGroup, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(E5_bray ~ ViralCopyNumberRanked, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(E5_bray ~ Timepoint, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(E5_bray ~ Cleared.Persist.Endpoint, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(E5_bray ~ DiseaseCINEndpoint, by= "margin", data = sampledf, permutations = 9999) # 


## Running Maaslin2 for determining differential abundance
# https://huttenhower.sph.harvard.edu/maaslin/
input_genus_data <- as.data.frame(read_csv("./E5_Genus_RELATIVEAbundance.csv"))
rownames(input_genus_data) <-input_genus_data[,1]

input_meta_file <- as.data.frame(E5.meta)
input_meta_file <- subset.data.frame(input_meta_file, MouseID != "PBS")
input_meta_file <- subset.data.frame(input_meta_file, MouseID != "NegEX")
#rownames(input_meta_file)<-input_meta_file[,1]
input_meta_file$ViralInfectGroup = factor(input_meta_file$ViralInfectGroup, levels = c("Mock", "10^4","10^6", "10^8"))

#Supplemental Figure 3B (figure and associated stats)
fit_data = Maaslin2( 
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2.FINAL/E5_ViralInfectDose_Genus_maslin", 
  fixed_effects = c("ViralInfectGroup"),
  reference = c("ViralInfectGroup", "Mock"),
  random_effects = c("MouseID"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2.FINAL/E5_ViralInfectDose_Genus_10.4_maslin", 
  fixed_effects = c("ViralInfectGroup"),
  reference = c("ViralInfectGroup", "10^4"),
  random_effects = c("MouseID"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2.FINAL/E5_ViralInfectDose_Genus_10.6_maslin", 
  fixed_effects = c("ViralInfectGroup"),
  reference = c("ViralInfectGroup", "10^6"),
  random_effects = c("MouseID"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2.FINAL/E5_ViralInfectDose_Genus_10.8_maslin", 
  fixed_effects = c("ViralInfectGroup"),
  reference = c("ViralInfectGroup", "10^8"),
  random_effects = c("MouseID"))

# Supplemental Figure 3C
Infected_meta <- subset.data.frame(input_meta_file, ViralInfectGroup != "Mock")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Infected_meta, 
  output = "./MAASLIN2.FINAL/E5_Infected_Only_ViralHL_genus_maslin", 
  fixed_effects = c("ViralCopyNumberRanked"),
  reference = c("ViralCopyNumberRanked", "Low"),
  random_effects = c("MouseID"))

#Additional MAASLIN plots (not in paper)
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Infected_meta, 
  output = "./MAASLIN2.FINAL/E5__Infected_only_ClearedPersist_Genus_maslin", 
  fixed_effects = c("Cleared.Persist.Endpoint"),
  reference = c("Cleared.Persist.Endpoint", "Cleared"),
  random_effects = c("MouseID"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Infected_meta, 
  output = "./MAASLIN2.FINAL/E5_Infected_Only_Disease_CIN1_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "CIN1"),
  random_effects = c("MouseID"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Infected_meta, 
  output = "./MAASLIN2.FINAL/E5_Infected_Only_Disease_CIN2_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "CIN2"),
  random_effects = c("MouseID"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Infected_meta, 
  output = "./MAASLIN2.FINAL/E5_Infected_Only_Disease_CIN3_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "CIN3"),
  random_effects = c("MouseID"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Infected_meta, 
  output = "./MAASLIN2.FINAL/E5_Infected_Only_Disease_SCC_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "SCC"),
  random_effects = c("MouseID"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Infected_meta, 
  output = "./MAASLIN2.FINAL/E5_Infected_Only_Disease_SCC_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "Normal"),
  random_effects = c("MouseID"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2.FINAL/E5_ALLmice_Disease_CIN1_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "CIN1"),
  random_effects = c("MouseID"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2.FINAL/E5_ALLmice_Disease_CIN2_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "CIN2"),
  random_effects = c("MouseID"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2.FINAL/E5_AllMice_Disease_CIN3_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "CIN3"),
  random_effects = c("MouseID"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2.FINAL/E5_AllMice_Disease_SCC_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "SCC"),
  random_effects = c("MouseID"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2.FINAL/E5_AllMice_Disease_SCC_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "Normal"),
  random_effects = c("MouseID"))

## Mix omics for PLS-DA analysis - GENUS level 
# Figure 3 and Supplemental Figure 3 
ps.E5clean <- tax_glom(ps.E5, "Genus")
ps.E5clean <- phyloseq_filter_prevalence(ps.E5clean, prev.trh = 0.05, abund.trh = 10, threshold_condition = "AND") # using AND here b/c need to ensure that the ASV is at least in 2 (>0.04 % of samples)
relE5.2 <- transform(ps.E5clean, "compositional") 
relE5.2.df <-psmelt(relE5.2)
relE5.2.df<- subset.data.frame(relE5.2.df, Abundance > 0.001)

write.csv(relE5.2.df, "./Mixomics_Genus/E5relativeAbundance_Genus_Mixomics.csv")
E5.relativeAbundanceOTHER_pivot <- read_csv("./Mixomics_Genus/E5relativeAbundance_Genus_Mixomics.Pivot.csv")

Micro.df <- as.data.frame(E5.relativeAbundanceOTHER_pivot) # 77 ASVs included 
rownames(Micro.df) <- Micro.df$SampleID
Micro.df <- Micro.df[,-1]

E5.meta.df <- data.frame(sample_data(ps.E5))

# PLSDA for ALL infected vs. Mock
E5.Meta<- subset.data.frame(E5.meta.df , ViralInfectGroup %in% c("Mock", "10^4", "10^6", "10^8"))
Micro.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E5.Meta))

X <- Micro.df
Y <- E5.Meta
Y.I <- Y$ViralInfectGroup
Y.I<-factor(Y.I, levels = c("Mock", "10^4", "10^6", "10^8"))

plsda.I <- plsda(X, Y.I, ncomp = 3)
plotIndiv(plsda.I, 
          ind.names = F, 
          ellipse = TRUE,
          comp = c(1,2),
          legend=TRUE, cex=6,
          col =c("#21777C","#8AB17D","#EFB366","#DE7350"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 4 Infection v. none')

cord1_Load<-plotLoadings(plsda.I, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Mock== "TRUE" & Contrib.10.4 == "TRUE" &Contrib.10.6 == "TRUE" &Contrib.10.8 == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Mock", "X10.4","X10.6","X10.8", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.I, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings,  Contrib.Mock== "TRUE" & Contrib.10.4 == "TRUE" &Contrib.10.6 == "TRUE" &Contrib.10.8 == "TRUE"  ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Mock", "X10.4","X10.6","X10.8", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))
Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.1)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.4)+
  scale_color_manual(values = c( "#21777C","#8AB17D","#EFB366","#DE7350"))+
  scale_fill_manual(values = c("#21777C","#8AB17D","#EFB366","#DE7350"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Mock", "10^4", "10^6", "10^8")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 4 -Viral Infection Inoculation Dose groups top driving Taxa')

#disease severity ALL mice 
E5.CIN.Meta<- subset.data.frame(E5.meta.df , DiseaseCINEndpoint %in% c("Normal", "CIN1", "CIN2", "CIN3", "SCC"))
Micro.CIN.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E5.CIN.Meta))

X.CIN <- Micro.CIN.df
Y <- E5.CIN.Meta
Y.CIN <- Y$DiseaseCINEndpoint
Y.CIN<-factor(Y.CIN, levels = c("Normal", "CIN1", "CIN2", "CIN3", "SCC"))

plsda.CIN<- plsda(X.CIN, Y.CIN, ncomp = 3)
plotIndiv(plsda.CIN, 
          ind.names = F, 
          ellipse = TRUE,
          comp = c(2,3),
          legend=TRUE, cex=6,
          col =c("#BBBBBB", "#21777C","#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 4 All mice Disease Severity')

cord1_Load<-plotLoadings(plsda.CIN, comp = 2, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Normal== "TRUE" & Contrib.CIN1 == "TRUE" &Contrib.CIN2 == "TRUE" &Contrib.CIN3 == "TRUE" &Contrib.SCC =="TRUE") # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Normal", "CIN1","CIN2","CIN3","SCC", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.CIN, comp = 3, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings,  Contrib.Normal== "TRUE" & Contrib.CIN1 == "TRUE" &Contrib.CIN2 == "TRUE" &Contrib.CIN3 == "TRUE" &Contrib.SCC== "TRUE"  ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Normal", "CIN1","CIN2","CIN3","SCC", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))
Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.19)
TopCord<-subset.data.frame(TopCord, Variable %in% c("Acinetobacter", "Cutibacterium","Pseudomonas","Herbaspirillum","Corynebacterium","Actinobacillus","Brochothrix","Rhodococcus", "Kocuria", "Staphylococcus", "Stenotrophomonas"))
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.4)+
  scale_color_manual(values = c("#21777C","#8AB17D","#F6D27A","#F5A261"))+
  scale_fill_manual(values = c( "#21777C","#8AB17D","#F6D27A","#F5A261"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Normal","CIN1","CIN2", "CIN3","SCC")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 4 -All mice Disease Severity')


# INFECTED mice only viral load H L 
E5.I.Meta<- subset.data.frame(E5.Meta , ViralInfectGroup != "Mock")

E5.I.VHL.Meta<- subset.data.frame(E5.I.Meta , ViralCopyNumberRanked %in% c("Low", "High"))
Micro.I.VHL.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E5.I.VHL.Meta))

X.VHL <-Micro.I.VHL.df
Y.VHL <- E5.I.VHL.Meta$ViralCopyNumberRanked
Y.VHL <- factor(Y.VHL, levels = c("Low", "High" ))

E5.I.PC.Meta<- subset.data.frame(E5.I.Meta, Cleared.Persist.Endpoint %in% c("Cleared", "Persistent"))
Micro.I.PC.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E5.I.PC.Meta))

X.PC<-Micro.I.PC.df
Y.PC <- E5.I.PC.Meta$Cleared.Persist.Endpoint
Y.PC <- factor(Y.PC, levels = c("Cleared", "Persistent" ))

plsda.VHL <- plsda(X.VHL, Y.VHL, ncomp = 3)
plotIndiv(plsda.VHL, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp = c(1,2),
          col =c(  "#2A9D8F", "#053061"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 4 infected mice only Viral Load')

cord1_Load<-plotLoadings(plsda.VHL, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Low== "TRUE" & Contrib.High == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Low", "High", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.VHL, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Low== "TRUE" & Contrib.High == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Low", "High", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.001)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.5)+
  scale_color_manual(values = c("#053061","#2A9D8F"))+
  scale_fill_manual(values = c("#053061", "#2A9D8F"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Low", "High")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 4 infected mice - Viral Load, top driving Taxa')


plsda.PC <- plsda(X.PC, Y.PC, ncomp = 3)
plotIndiv(plsda.PC, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp = c(1,2),
          col =c("#8AB17D",  "#2A9D8F"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 4 infected mice only Viral persist and clearance')

cord1_Load<-plotLoadings(plsda.PC, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Cleared== "TRUE" & Contrib.Persistent == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Cleared", "Persistent", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.PC, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Cleared== "TRUE" & Contrib.Persistent == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Cleared", "Persistent", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.001)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .02)+
  xlim(-0,0.4)+
  scale_color_manual(values = c("#8AB17D","#2A9D8F"))+
  scale_fill_manual(values = c("#8AB17D", "#2A9D8F"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Cleared", "Persistent")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 4 infected mice - Viral persistance and clearance, top driving Taxa')

# disease severity Infected Mice 
E5.CIN.Meta<- subset.data.frame(E5.I.Meta , DiseaseCINEndpoint %in% c("Normal", "CIN1", "CIN2", "CIN3", "SCC"))
Micro.CIN.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E5.CIN.Meta))

X.CIN <- Micro.CIN.df
Y <- E5.CIN.Meta
Y.CIN <- Y$DiseaseCINEndpoint
Y.CIN<-factor(Y.CIN, levels = c("Normal", "CIN1", "CIN2", "CIN3", "SCC"))

plsda.CIN<- plsda(X.CIN, Y.CIN, ncomp = 3)
plotIndiv(plsda.CIN, 
          ind.names = F, 
          ellipse = TRUE,
          comp = c(1,2),
          legend=TRUE, cex=6,
          col =c("#BBBBBB", "#21777C","#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 4 Infected Mice Disease Severity')

cord1_Load<-plotLoadings(plsda.CIN, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Normal== "TRUE" & Contrib.CIN1 == "TRUE" &Contrib.CIN2 == "TRUE" &Contrib.CIN3 == "TRUE" &Contrib.SCC =="TRUE") # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Normal", "CIN1","CIN2","CIN3","SCC", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.CIN, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings,  Contrib.Normal== "TRUE" & Contrib.CIN1 == "TRUE" &Contrib.CIN2 == "TRUE" &Contrib.CIN3 == "TRUE" &Contrib.SCC== "TRUE"  ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Normal", "CIN1","CIN2","CIN3","SCC", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))
Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.01)
#TopCord<-subset.data.frame(TopCord, Variable %in% c("Acinetobacter", "Cutibacterium","Pseudomonas","Herbaspirillum","Corynebacterium","Actinobacillus","Brochothrix","Rhodococcus", "Kocuria", "Staphylococcus", "Stenotrophomonas"))
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.4)+
  scale_color_manual(values = c("#21777C","#8AB17D","#F6D27A","#F5A261"))+
  scale_fill_manual(values = c( "#21777C","#8AB17D","#F6D27A","#F5A261"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Normal","CIN1","CIN2", "CIN3","SCC")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 4-Infected Mice isease Severity')


# Additional ploots - CORE TAXA for infection dose groups in Experiment 4 
# https://microbiome.github.io/tutorials/core_venn.html
ps.E5REl # the relative abundance phyloseq file 
ps.E5REl_Genera <- tax_glom(ps.E5REl, taxrank = "Genus") # merging ASVs at the genera level (the OTU numbers indicated are for ONE genus only - not for unique reads)
# NOTE the OTU numbers are just ONE of the OTUs for each of the genera. these numbers will not be all the OTUs in the origional dataset

table(meta(ps.E5REl_Genera)$ViralInfectGroup) # 62 E2 samples and 34 E3 samples

VIG <- unique(as.character(meta(ps.E5REl_Genera)$ViralInfectGroup))

# fore loop to go through the experiments and combine identified taxa into a core taxa list 
list_core <- c() # an empty object to store information
for (n in VIG){ # for each variable n in Ex23 (the ex 2 and ex 3)
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps.E5REl_Genera, VIG == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.0001 in at least 30% samples 
                         prevalence = 0.5)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core)

mycols <- c("#21777C","#8AB17D","#EFB366","#DE7350") # colors for the ven diagram
plot(venn(list_core),
     fills = mycols)
plot(euler(list_core),
     fills = mycols)

# format names
ps.E5REl_Genera.f <- format_to_besthit(ps.E5REl_Genera)
# check names
taxa_names(ps.E5REl_Genera.f)[1:5]

list_core <- c() # an empty object to store information
for (n in VIG){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps.E5REl_Genera.f, VIG == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.00001, # 0.0001 in at least 50% samples 
                         prevalence = 0.5)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core)


