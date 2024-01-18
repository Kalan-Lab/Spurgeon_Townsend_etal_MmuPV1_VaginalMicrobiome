#  NATURAL Vaginal Microbiome -- Experiments 1 (pilot) and experiment 2 (Baseline)

library(tidyverse)
library(phyloseq)
library(microbiomeutilities) # for the core microbiome analysis 
library(eulerr)
library(ggplot2)
library(microbiome)
library(metagMisc)
library(qiime2R)
library(xlsx)
library(RColorBrewer)
library(microshades)
library(decontam)
library(vegan)
library(DESeq2)
library(ANCOMBC) 
library(Maaslin2)
library(pheatmap)
library(mixOmics)

setwd("/Users/liztown/Documents/KalanLab/Papers/Vaginal Microbiome/Rcode_forManuscript/Experiment1_Pilot")

ps.E2pilot <- qza_to_phyloseq(features = "./table-Mouse-dada2-20-160.qza",
                              taxonomy = "./taxonomy-Mouse.qza",
                              metadata = "./E2PilotMouseManufest.txt",
                              tree = "./MousePilot-rooted-tree.qza")

setwd("/Users/liztown/Documents/KalanLab/Papers/Vaginal Microbiome/Rcode_forManuscript/Experiment2_Baseline")

ps.E3Base <- qza_to_phyloseq(features = "./table-Baseline-dada2.qza",
                             taxonomy = "./taxonomy-Baseline.qza",
                             metadata = "./E3BaselineManufest.txt",
                             tree = "./Baseline-rooted-tree.qza")


setwd("/Users/liztown/Documents/KalanLab/Papers/Vaginal Microbiome/Rcode_forManuscript/NaturalMicrobiome_Experiments1and2")

### decontaminating all the files and making csv for the TOTAL abundance in each of the data sets
ps.Decontam <- subset_taxa(ps.E2pilot, Family !="mitochondria")
ps.Decontam <- subset_taxa(ps.Decontam, Class !="Chloroplast")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !="Cyanobacteria")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !=" ")
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylobacterium") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingomonas") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ralstonia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Facklamia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Chryseobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curvibacter") 
ps.E2pilotDe <- phyloseq_filter_prevalence(ps.Decontam, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR")
ps.E2pilotDe <- subset_samples(ps.E2pilotDe, SampleOrigin == "Mouse") 
ps.E2pilotDe# THE CLEANED FILE 

ps.Decontam <- subset_taxa(ps.E3Base, Family !="mitochondria")
ps.Decontam <- subset_taxa(ps.Decontam, Class !="Chloroplast")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !="Cyanobacteria")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !=" ")
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylobacterium") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingomonas") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ralstonia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Facklamia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Chryseobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curvibacter") 
ps.E3BaseDe <- phyloseq_filter_prevalence(ps.Decontam, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR")
ps.E3BaseDe <- subset_samples(ps.E3BaseDe, SampleOrigin == "Mouse") # removing the negative control
ps.E3BaseDe# THE CLEANED FILE 

# PART 1 - The NATURAL mouse vaginal microbiome
### this portion combines data from the untouched mice (so no mock infection [mock infection has depo], no depo, no nothing) from the E2 (pilot) and E3 (baseline experiments)
ps.E2natural <- subset_samples(ps.E2pilotDe, NaturalMicro == "Natural")
ps.E3natural <- subset_samples(ps.E3BaseDe, NaturalMicro == "Natural")

ps.Natural <- merge_phyloseq(ps.E2natural, ps.E3natural)

OutCV <- psmelt(ps.Natural)
write.csv(OutCV, "/Users/liztown/Documents/KalanLab/VaginalMicrobiome/Combined16Sresults/NaturalVaginalMicrobiome/E2E3NaturalMicrobiomeAbundanceTotalsOut.csv")

ps.NaturalREl <- transform(ps.Natural, "compositional")
NatRelout <- psmelt(ps.NaturalREl)
write.csv(NatRelout, "/Users/liztown/Documents/KalanLab/VaginalMicrobiome/Combined16Sresults/NaturalVaginalMicrobiome/E2E3NaturalMicrobiomeRELATIVEAbundanceTotalsOut.csv")

NaturalRel.df<-as.data.frame(NatRelout)
Natural.meta <- data.frame(sample_data(ps.Natural))
Natural.meta$SampleID <- row.names(Natural.meta)

# Relative abundance # figure 1
ps.NaturalREl <- microbiome::transform(ps.Natural, "compositional")
ps.NaturalREl_Merge = tax_glom(ps.NaturalREl, "Genus")
NatRelout <- psmelt(ps.NaturalREl_Merge)

OtherNat <- as.data.frame(NatRelout) # the dataframe of the relative abundance within each sample
OtherNat <- subset.data.frame(OtherNat, Abundance >=0.0001)
OtherNat <- OtherNat %>% mutate(Genus = ifelse(Abundance < 0.05, "Other", Genus)) # replacing the genus catagory with "Other" if the abundance is less than 1% 
OtherNat <- OtherNat %>% mutate(Phylum = ifelse(Abundance < 0.05, "Other", Phylum)) 
OtherNat<- OtherNat %>% mutate(Genus = ifelse(Genus == "[Prevotella]","Prevotella", Genus))
table(OtherNat$Genus, OtherNat$Phylum)

OtherNat$Phylum <- factor(OtherNat$Phylum, levels = c("Other", "Actinobacteria",  "Firmicutes", "Bacteroidetes", "Proteobacteria"))
OtherNat$Genus<- factor(OtherNat$Genus, levels = c("Other",
                                                   "Corynebacterium","Dermacoccus", "Micrococcus", 
                                                   "Propionibacterium","Rhodococcus","Rubrobacter", "Streptomyces", ##Actinobacteria / actinomycetota (7)
                                                   
                                                   "Aerococcus",  "Alloiococcus", "Anaerobacillus", "Bacillus","Clostridium","Coprococcus", "Enterococcus","Gemella", "Lactobacillus", "Lactococcus",
                                                   "Leuconostoc", "Oscillospira",  "Sporosarcina","Staphylococcus", "Streptococcus","Turicibacter", # Firmicutes/ Bacillota (16)
                                                   
                                                   "Bacteroides","Flavisolibacter","Odoribacter",  "Porphyromonas",  "Prevotella",# Bacteroidetes / Bacteroidota (5)
                                                   
                                                   "Acinetobacter","Aminobacter","Caulobacter","Cellvibrio","Halomonas", "Massilia","Mycoplana", "Neisseria", 
                                                   "Paracoccus","Perlucidibaca", "Pseudomonas","Roseomonas", "Sphingopyxis", "Stenotrophomonas", "Xylella",#Proteobacteria /pseudomonadota (15)
                                                   "SargSea-WGS")) #SAR406 

ggplot(OtherNat, aes(x= Ex, y =Abundance, fill = Genus, color = Genus))+
  geom_bar(stat="identity", position="fill") +
  facet_wrap("Ex", scales = "free_x", nrow = 1)+
  ylab("Relative Abundance")+
  theme_light()+
  scale_color_manual(values = c("#C1C1C1", colorRampPalette(c("#053061", "#175169","#21777C", "#2A9D8F", "#8AB17D","#C0C27B",
                                                              "#F5D279","#EFB366","#F4A261","#DE7350","#C8443F", "#913C44","#6E2C33" ))(21)))+
  scale_fill_manual(values = c("#C1C1C1", colorRampPalette(c("#053061", "#175169","#21777C", "#2A9D8F", "#8AB17D","#C0C27B",
                                                             "#F5D279","#EFB366","#F4A261","#DE7350","#C8443F", "#913C44","#6E2C33" ))(21)))+
  theme(axis.text.x = element_text(size = 4)) +
  ggtitle("Genus present > 0.5% of reads in a sample")


# core microbiome analysis - between hte E2 and E3 experiment Figure 1D 
# https://microbiome.github.io/tutorials/core_venn.html
ps.NaturalREl # the relative abundance phyloseq file 
ps.NatRelGenera <- tax_glom(ps.NaturalREl, taxrank = "Genus") # merging ASVs at the genera level (the OTU numbers indicated are for ONE genus only - not for unique reads)
table(meta(ps.NatRelGenera)$Ex) # 62 E2 samples and 34 E3 samples

Ex23 <- unique(as.character(meta(ps.NatRelGenera)$Ex))
print(Ex23)

# fore loop to go through the experiments and combine identified taxa into a core taxa list 
list_core <- c() # an empty object to store information
for (n in Ex23){ # for each variable n in Ex23 (the ex 2 and ex 3)
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps.NatRelGenera, Ex == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.0001 in at least 50% samples 
                         prevalence = 0.5)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core)

mycols <- c(E2="#F5D279", E3="#F4A261") # colors for the ven diagram
plot(venn(list_core),
     fills = mycols)

# format names
ps.NatRelGenera.f <- format_to_besthit(ps.NatRelGenera)
# check names
taxa_names(ps.NatRelGenera.f)[1:5]

list_core <- c() # an empty object to store information

for (n in Ex23){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps.NatRelGenera.f, Ex == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.0001 in at least 50% samples 
                         prevalence = 0.5)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core)

#Alpha abundance of the natural vaginal microbiome (not shown in paper )
AlphaTab <-microbiome::alpha(ps.Natural, index = c("observed", "chao1", "diversity_inverse_simpson", "diversity_gini_simpson","diversity_shannon",	"diversity_coverage", "evenness_camargo",	"evenness_pielou",
                                                   "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp", "dominance_dmn", "dominance_absolute", "dominance_relative", "dominance_simpson", 
                                                   "dominance_core_abundance", "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_rare_abundance"))
AlphaTable<-as.data.frame(AlphaTab)
AlphaTable$SampleID <- row.names(AlphaTable)
rownames(AlphaTable) <- NULL

NatAlphaTable <- merge(Natural.meta, AlphaTable, by = "SampleID", all = T )

ggplot(NatAlphaTable, aes(x = Ex, y =diversity_shannon, color = Ex))+
  geom_jitter(size = 4)+
  stat_summary(fun.y = mean,  geom = "point",shape = 95, size = 15, show.legend = F)+ 
  scale_color_manual(values = c("#EFB366","#DE7350"))+
  scale_shape_manual(values = c(19,17))+
  ylim(0,6)+
  theme_light()+
  ggtitle("Shannon Diversity of the Natural Mouse Vaginal Microbiome")

# beta diversity of the natural microbiome 
min(sample_sums(ps.Natural))# minimum sample read is 11
median(sample_sums(ps.Natural)) # 6835
max(sample_sums(ps.Natural)) #63016
table(sample_sums(ps.Natural))

rarecurve(t(otu_table(ps.Natural)), step=100, ylim =c(0,100), xlim=c(0,1100))  ## considering using a read cut off of 5000 for beta diversity metrics 

ps.NatRarefied = rarefy_even_depth(ps.Natural, rngseed=1, sample.size=700, replace=F) # 9 samples removed due to too few reads 

GP = ps.NatRarefied
GP.ord <- ordinate(GP, "NMDS",  "bray") # For Bray Curtis of the rareified dataset
plot_ordination(GP, GP.ord, type="samples", color="Ex")+# shape = "EdgeCenter") + 
  geom_point(size=3) + ggtitle("Bray Curtis")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(size = 6)+
  scale_shape_manual(values = c(3,4,19,17))+
  scale_color_manual(values = c("#EFB366","#DE7350"))+
  theme_bw()

sampledf <- data.frame(sample_data(ps.NatRarefied)) 
Nat_bray <- phyloseq::distance(ps.NatRarefied, method = "bray")
Bray_clust <-hclust(Nat_bray)
adonis2(Nat_bray ~ Ex, by= "margin", data = sampledf, permutations = 9999) # p-value < 1e-4 ***
adonis2(Nat_bray ~ Ex+ NaturalMicro, by= "margin", data = sampledf, permutations = 9999) #ex P-Value = 0.0001**, Natural Micro P.value = 0.0069 **


# core microbiome analysis - between hte experiment 1 and experiment2 NATURAL MICROBIOME ONLY Figure 1D
# https://microbiome.github.io/tutorials/core_venn.html
ps.NaturalREl # the relative abundance phyloseq file 
ps.NaturalRElNat <- subset_samples(ps.NaturalREl, NaturalMicro == "Natural")
ps.NatRelGenera <- tax_glom(ps.NaturalRElNat, taxrank = "Genus") # merging ASVs at the genera level (the OTU numbers indicated are for ONE genus only - not for unique reads)
# NOTE the OTU numbers are just ONE of the OTUs for each of the genera. these numbers will not be all the OTUs in the origional dataset

table(meta(ps.NatRelGenera)$Ex) # 62 E2 samples and 34 E3 samples

Ex23 <- unique(as.character(meta(ps.NatRelGenera)$Ex))
print(Ex23)

# fore loop to go through the experiments and combine identified taxa into a core taxa list 
list_core <- c() # an empty object to store information
for (n in Ex23){ # for each variable n in Ex23 (the ex 2 and ex 3)
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps.NatRelGenera, Ex == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.0001 in at least 50% samples 
                         prevalence = 0.3)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core)

mycols <- c(E2="#F5D279", E3="#F4A261") # colors for the ven diagram
plot(venn(list_core),
     fills = mycols)

# format names
ps.NatRelGenera.f <- format_to_besthit(ps.NatRelGenera)
# check names
taxa_names(ps.NatRelGenera.f)[1:5]

list_core <- c() # an empty object to store information

for (n in Ex23){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps.NatRelGenera.f, Ex == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.0001 in at least 50% samples 
                         prevalence = 0.3)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core)

# core microbiome analysis - between hte experiment 1 and 2 experiment POST-DEPO only Figure 1E
# https://microbiome.github.io/tutorials/core_venn.html
ps.E2PD <- subset_samples(ps.E2pilotDe, NaturalMicro == "Post-Depo")
ps.E2PD <- subset_samples(ps.E2PD, AntibioticTx %in% c("None"))
ps.E3PD <- subset_samples(ps.E3BaseDe, AntibioticTx %in% c("None")) # untouxhed mice in pre-depo timrpoints = "Natural", Untouched mice at post -depo timepoimt = "Untouched'
ps.E3PD <- subset_samples(ps.E3PD, NaturalMicro == "Post-Depo")
ps.PD <- merge_phyloseq(ps.E2PD , ps.E3PD)
ps.PDREl <- transform(ps.PD, "compositional")

ps.PDRelGenera <- tax_glom(ps.PDREl, taxrank = "Genus") # merging ASVs at the genera level (the OTU numbers indicated are for ONE genus only - not for unique reads)
# NOTE the OTU numbers are just ONE of the OTUs for each of the genera. these numbers will not be all the OTUs in the origional dataset

table(meta(ps.PDRelGenera)$Ex) # 62 E2 samples and 34 E3 samples

Ex23 <- unique(as.character(meta(ps.PDRelGenera)$Ex))
print(Ex23)

# fore loop to go through the experiments and combine identified taxa into a core taxa list 
list_core <- c() # an empty object to store information
for (n in Ex23){ # for each variable n in Ex23 (the ex 2 and ex 3)
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps.PDRelGenera, Ex == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.0001 in at least 50% samples 
                         prevalence = 0.4)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core)

mycols <- c(E2="#F5D279", E3="#F4A261") # colors for the ven diagram
plot(venn(list_core),
     fills = mycols)

# format names
ps.PDRelGenera.f <- format_to_besthit(ps.PDRelGenera)
# check names
taxa_names(ps.PDRelGenera.f)[1:5]

list_core <- c() # an empty object to store information

for (n in Ex23){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps.PDRelGenera.f, Ex == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.0001 in at least 50% samples 
                         prevalence = 0.4)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core)


## MAASLIN # Difference between E2 and E2 in the NATURAL microbiome 
input_genus_data <- as.data.frame(read_csv("./Natural_genus_RELATIVEAbundanceTotalsOut.csv"))
rownames(input_genus_data) <-input_genus_data[,1]

input_meta_file <-Natural.meta
Nat_input_meta_file <-subset.data.frame(input_meta_file, NaturalMicro == "Natural")
Nat_input_meta_file$Ex <- factor(Nat_input_meta_file$Ex, levels = c("E2", "E3"))

Nat_PDN_meta_file <- subset.data.frame(Nat_meta_file, NaturalMicro != "Untouched")
E2_meta_file <- subset.data.frame(input_meta_file, Ex == "E2")
E3_meta_file <- subset.data.frame(input_meta_file, Ex == "E3")
E3_PDN_meta_file <- subset.data.frame(E3_meta_file, NaturalMicro != "Untouched")
E3_PDN_NoUNTOUCHED_MF <- subset.data.frame(E3_meta_file, Infected != "Untouched")

# Evaluating differences in the natural microbiome between E2 and E3 
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Nat_input_meta_file, 
  output = "./MAASLIN/Natural_Microbiome_E2vE3", 
  fixed_effects = c("Ex"),
  reference = c("Ex", "E2"),
  random_effects = c("MP"))

# Evaluating differences in the Post-Depo and Natural in E2 and E3 TOGETHER (NO Associations)
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Nat_input_meta_file, 
  output = "./MAASLIN/E2E3_PostDepo_v_Natural_Microbiome", 
  fixed_effects = c("NaturalMicro"),
  reference = c("NaturalMicro", "Natural"),
  random_effects = c("Ex"))

# Evaluating differences in the Post-Depo and Natural in E2 and E3 separatly 
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = E2_meta_file, 
  output = "./MAASLIN/E2_PostDepo_v_Natural_Microbiome", 
  fixed_effects = c("NaturalMicro"),
  reference = c("NaturalMicro", "Natural"),
  random_effects = c("MP"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = E3_meta_file, 
  output = "./MAASLIN/E3_PostDepo_v_Natural_v_untouched_Microbiome", 
  fixed_effects = c("NaturalMicro"),
  reference = c("NaturalMicro", "Natural"),
  random_effects = c("MP"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = E3_PDN_meta_file, 
  output = "./MAASLIN/E3_PostDepo_v_Natural_Microbiome", 
  fixed_effects = c("NaturalMicro"),
  reference = c("NaturalMicro", "Natural"),
  random_effects = c("MP"))

fit_data = Maaslin2( # no untouched even in the natural group
  input_data = input_genus_data, 
  input_metadata = E3_PDN_NoUNTOUCHED_MF, 
  output = "./MAASLIN/E3_PostDepo_v_Natural_NOUNTOUCHEDatBaseline_Microbiome", 
  fixed_effects = c("NaturalMicro"),
  reference = c("NaturalMicro", "Natural"),
  random_effects = c("MP"))


