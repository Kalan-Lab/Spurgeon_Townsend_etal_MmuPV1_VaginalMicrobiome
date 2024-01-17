# Core microbiome analysis for INFECTED mice from Experiments 1-4 (vaginal levages of the microbiome)
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

setwd("/Users/liztown/Documents/KalanLab/Papers/Vaginal Microbiome/Rcode_forManuscript")
#Importing EACH experiment / QIIME object 
ps.E2pilot <- qza_to_phyloseq(features = "./Experiment1_Pilot/table-Mouse-dada2-20-160.qza",
                              taxonomy = "./Experiment1_Pilot/taxonomy-Mouse.qza",
                              metadata = "./Experiment1_Pilot/E2PilotMouseManufest.txt")# ,
ps.E3Base <- qza_to_phyloseq(features = ".Experiment2_Baseline/table-Baseline-dada2.qza",
                             taxonomy = "./Experiment2_Baseline/taxonomy-Baseline.qza",
                             metadata = "./Experiment2_Baseline/E3BaselineManufest.txt") #,
ps.E4infectFirst <-qza_to_phyloseq(features = "./Experiment3_InfectFirst/table-E4-dada2.qza",
                                   taxonomy = "./Experiment3_InfectFirst/taxonomy-E4.qza",
                                   metadata = "./Experiment3_InfectFirst/E4InfectFirstManufest.txt")
ps.E5infectDose <- qza_to_phyloseq(features = "./Experiment4_InfectDose/table-E5-dada2.qza",
                                   taxonomy = "./Experiment4_InfectDose/taxonomy-E5.qza",
                                   metadata = "./Experiment4_InfectDose/E5manufest.txt")


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

ps.Decontam <- subset_taxa(ps.E4infectFirst, Family !="mitochondria")
ps.Decontam <- subset_taxa(ps.Decontam, Class !="Chloroplast")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !="Cyanobacteria")
ps.Decontam <- subset_taxa(ps.Decontam, Phylum !=" ")
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Methylobacterium") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Sphingomonas") # Highly previlent contaminants
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Ralstonia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Facklamia") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Chryseobacterium") 
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Curvibacter") 
ps.E4 <- subset_samples(ps.Decontam, !sample_names(ps.Decontam) %in% c("LK16S009-241", "LK16S009-267", "LK16S009-267", "LK16S009-293", "LK16S009-318", "LK16S009-343", "LK16S009-368", "LK16S009-369", "LK16S009-370")) # removing negative controls 
ps.E4 <- phyloseq_filter_prevalence(ps.E4, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR") #This function will remove taxa (OTUs) with low prevalence, where prevalence is the fraction of total samples in which an OTU is observed.

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


### CORE Microbiome for INFECTED MICE Figure 2 I
# PART 1 - The NATURAL mouse vaginal microbiome  and POST Depo 
### this portion combines data from the untouched mice (so no mock infection [mock infection has depo], no depo, no nothing) from the E2 (pilot) and E3 (baseline experiments)
ps.E2PostInfect <- subset_samples(ps.E2pilotDe, InfectStage2 %in% c("Infected"))
ps.E2PostInfect <- subset_samples(ps.E2PostInfect, AntibioticTx %in% c("None"))
ps.E2PostInfect <- subset_samples(ps.E2PostInfect, InfectStage2== "Infected")

ps.E3PostInfect <- subset_samples(ps.E3BaseDe, AntibioticTx %in% c("None")) # untouxhed mice in pre-depo timrpoints = "Natural", Untouched mice at post -depo timepoimt = "Untouched'
ps.E3PostInfect <- subset_samples(ps.E3PostInfect, Infected %in% c("Infected"))
ps.E3PostInfect <- subset_samples(ps.E3PostInfect, InfectStage2== "Infected")

ps.E4PostInfect <- subset_samples(ps.E4, ABX %in% c("None")) # untouxhed mice in pre-depo timrpoints = "Natural", Untouched mice at post -depo timepoimt = "Untouched'
ps.E4PostInfect <- subset_samples(ps.E4PostInfect, Infected %in% c("Infected"))

ps.E5PostInfect <- subset_samples(ps.E5, ViralInfectGroup %in% c("10^4", "10^6", "10^8")) # untouxhed mice in pre-depo timrpoints = "Natural", Untouched mice at post -depo timepoimt = "Untouched'


ps.PostInfect <- merge_phyloseq(ps.E2PostInfect , ps.E3PostInfect, ps.E4PostInfect, ps.E5PostInfect)
ps.PostInfectREl <- transform(ps.PostInfect, "compositional")

ps.PostInfectREl # the relative abundance phyloseq file 
ps.PostInfectRElGenera <- tax_glom(ps.PostInfectREl, taxrank = "Genus") # merging ASVs at the genera level (the OTU numbers indicated are for ONE genus only - not for unique reads)
# NOTE the OTU numbers are just ONE of the OTUs for each of the genera. these numbers will not be all the OTUs in the origional dataset

table(meta(ps.PostInfectRElGenera)$Ex) # 62 E2 samples and 34 E3 samples

Ex23 <- unique(as.character(meta(ps.PostInfectRElGenera)$Ex))
print(Ex23)

# fore loop to go through the experiments and combine identified taxa into a core taxa list 
list_core <- c() # an empty object to store information
for (n in Ex23){ # for each variable n in Ex23 (the ex 2 and ex 3)
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps.PostInfectRElGenera, Ex == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.0001, # 0.0001 in at least 50% samples 
                         prevalence = 0.3)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core)

mycols <- c(E2="#F5D279", E3="#F4A261", E4 ="#DE7350", E5 ="#21777C" ) # colors for the ven diagram
plot(venn(list_core),
     fills = mycols)

mycols2 <- c(E2="#F5D279", E3="#21777C", E4 ="#DE7350", E5 ="#F4A261" ) # colors for the ven diagram
plot(euler(list_core, shape = "ellipse"),
     fills = mycols2)

# format names
ps.PostInfectRElGenera.f <- format_to_besthit(ps.PostInfectRElGenera)
# check names
taxa_names(ps.PostInfectRElGenera.f)[1:5]

list_core <- c() # an empty object to store information

for (n in Ex23){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps.PostInfectRElGenera.f, Ex == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.05, # 0.0001 in at least 50% samples 
                         prevalence = 0.3)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
print(list_core)


