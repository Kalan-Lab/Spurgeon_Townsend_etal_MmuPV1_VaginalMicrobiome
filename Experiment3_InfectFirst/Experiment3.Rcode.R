# E4 analysis --- NO antibiotic tx only - looking only at the effects of the infection 
# NOTE there is NO pre-infection group - noo baseline. all timepoints are 4.5 wpi + to 24.5 wkp 
library(tidyverse)
library(phyloseq)
library(metagMisc)
library(ggplot2)
library(microbiome)
library(qiime2R)
library(xlsx)
library(RColorBrewer)
library(microshades)
library(decontam)
library(vegan)
library(DESeq2)
library(ANCOMBC) 
library(Maaslin2)
library(mixOmics)

setwd("/Users/liztown/Documents/KalanLab/Papers/Vaginal Microbiome/Rcode_forManuscript/Experiment3_InfectFirst")

ps.E4infectFirst <-qza_to_phyloseq(features = "./table-E4-dada2.qza",
                                   taxonomy = "./taxonomy-E4.qza",
                                   metadata = "./E4InfectFirstManufest.txt")


E4.meta <- data.frame(sample_data(ps.E4infectFirst))
E4.meta$SampleID <- row.names(E4.meta)
E4.NoABX.meta <- subset.data.frame(E4.meta, ABX != "ABX")

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

OutCV <- psmelt(ps.E4)
#write.csv(OutCV, "./E4_TotalAbundanceOut.csv")
ps.E4REl <- transform(ps.E4, "compositional")
E4REl <- psmelt(ps.E4REl)
#write.csv(E4REl, "./E4_RELATIVEAbundanceOut.csv")
E4Rel.df<-as.data.frame(E4REl)

# Relative abundance Experiment 3 - supplemental figure 2F 
OtherE4 <- E4Rel.df # the dataframe of the relative abundance within each sample
OtherE4 <- subset.data.frame(OtherE4, ABX == "None")
OtherE4 <- subset.data.frame(OtherE4, Genus != "NA")
OtherE4 <- subset.data.frame(OtherE4, Abundance > 0.01)
OtherE4 <- OtherE4 %>% mutate(Genus = ifelse(Abundance < 0.05, "Other", Genus)) # replacing the genus catagory with "Other" if the abundance is less than 1% 
OtherE4 <- OtherE4 %>% mutate(Phylum = ifelse(Abundance < 0.05, "Other", Phylum)) 
OtherE4$Infected <- factor(OtherE4$Infected, levels = c("Mock","Infected"))
table(OtherE4$Genus, OtherE4$Phylum)

OtherE4$Phylum <- factor(OtherE4$Phylum, levels = c("Other", "Actinobacteria",  "Firmicutes", "Bacteroidetes", "Fusobacteria","[Thermi]",  "Proteobacteria"))
OtherE4$Genus<- factor(OtherE4$Genus, levels = c("Other",
                                                 "Adlercreutzia","Bifidobacterium", "Brachybacterium","Brevibacterium", "Corynebacterium","Dietzia","Geodermatophilus", "Microbacterium", "Micrococcus",
                                                 "Propionibacterium","Rhodococcus","Rothia", "Streptomyces", ##Actinobacteria / actinomycetota (13)
                                                 
                                                 "Bacillus","Coprococcus", "Gemella","Lactobacillus","Leuconostoc", "Oscillospira",
                                                 "Ruminococcus","Sporosarcina", "Staphylococcus", "Streptococcus","Tetragenococcus","Turicibacter","Weissella", # Firmicutes/ Bacillota (13)
                                                 
                                                 "Parabacteroides", "Porphyromonas","Sediminibacterium", "Sphingobacterium", # Bacteroidetes / Bacteroidota 
                                                 "Fusobacterium", # fusobacteria
                                                 
                                                 "Achromobacter",  "Acinetobacter","Actinobacillus","Aeromonas", "Afipia","Brevundimonas","Burkholderia", "Caulobacter",
                                                 "Comamonas", "Delftia", "Enhydrobacter","Erwinia","Ewingella",  "Herbaspirillum", "Lysobacter", "Ochrobactrum", 
                                                 "Paracoccus", "Pseudomonas","Rahnella", "Rubellimicrobium","Serratia",  
                                                 "Stenotrophomonas")) #Proteobacteria /pseudomonadota (21)

BigMicro2 <-c("#C1C1C1", #Gray
              colorRampPalette(c("#031A34", "#053061", "#21777C", "#2A9D8F", "#8AB17D", "#C0C27B","#D2C77A", "#F5D279","#EFB366","#F4A261","#DE7350","#C8443F","#913C44", "#733036"))(53))

ggplot(OtherE4, aes(x= Timepoint, y =Abundance, fill = Genus, color = Genus))+
  geom_bar(stat="identity", position="fill") +
  facet_wrap("Infected", scales = "free_x", nrow = 1)+
  ylab("Relative Abundance")+
  theme_light()+
  scale_fill_manual(values = BigMicro2)+
  scale_color_manual(values = BigMicro2)+
  theme(axis.text.x = element_text(size = 4)) +
  ggtitle("Genus present > 0.5% of reads in a sample")


# Beta Diversity (not shown in paper)
min(sample_sums(ps.E4))# minimum sample read is 101
median(sample_sums(ps.E4)) # 1611
max(sample_sums(ps.E4)) #84870
ps.E4.Not <- subset_samples(ps.E4infectFirst, ABX != "ABX")
median(sample_sums(ps.E4.Not))
ps.E4.NoABX <- subset_samples(ps.E4, ABX != "ABX")
min(sample_sums(ps.E4.NoABX))# minimum sample read is 188
median(sample_sums(ps.E4.NoABX)) # 2373
max(sample_sums(ps.E4.NoABX)) #84870

rarecurve(t(otu_table(ps.E4.NoABX)), step=100, ylim =c(0,100), xlim=c(0,1000))  ## considering using a read cut off of 5000 for beta diversity metrics 

ps.E4ra.NoABX  = rarefy_even_depth(ps.E4.NoABX, rngseed=1, sample.size=600, replace=F) # THIS IS THE ONE TO GO WITH

GP = ps.E4ra.NoABX
GP.ord <- ordinate(GP, "NMDS",  "bray") # For Bray Curtis of the rareified dataset
plot_ordination(GP, GP.ord, type="samples", color="Group")+# shape = "EdgeCenter") + 
  geom_point(size=3) + ggtitle("Bray Curtis Experiment 3 - infected v. mock")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(aes(shape = Group), size = 8)+
  #scale_shape_manual(values = c(8,19,17,25,15,23))+
  scale_color_manual(values = c("#053061", "#8AB17D")) + #"#21777C","#8AB17D","#EFB366","#DE7350"))+
  scale_fill_manual(values = c("#053061", "#8AB17D")) + #values = c("#21777C","#8AB17D","#EFB366","#DE7350"))+
  theme_bw()

E4.inf.rar <- subset_samples(ps.E4ra.NoABX, Infected == "Infected")
GP = E4.inf.rar
GP.ord <- ordinate(GP, "NMDS",  "bray")
plot_ordination(GP, GP.ord, type="samples", color="InfectStageML")+# shape = "EdgeCenter") + 
  geom_point(size=3) + ggtitle("Bray Curtis Experiment 3 - infected v. mock")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(aes(shape = InfectStageML), size = 8)+
  scale_shape_manual(values = c(19,15,17))+
  scale_color_manual(values = c("#186275","#021327", "#053061")) + #"#21777C","#8AB17D","#EFB366","#DE7350"))+
  scale_fill_manual(values = c("#186275","#021327","#053061")) + #values = c("#21777C","#8AB17D","#EFB366","#DE7350"))+
  theme_bw()

#permanova For Table 2 
sampledf <- data.frame(sample_data(ps.E4ra.NoABX)) 
LCM_bray <- phyloseq::distance(ps.E4ra.NoABX, method = "bray")
Bray_clust <-hclust(LCM_bray)
plot(Bray_clust)

adonis2(LCM_bray ~ Infected  , by= "margin", data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~  ViralCopyNumberHighLow, data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~  Cleared.PersistentAtEndpoint, data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~  DiseaseCINEndpoint , data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~  InfectStageML , data = sampledf, permutations = 9999) # 

Ps.rare.LCMinf <-subset_samples(ps.E4ra.NoABX, Infected == "Infected")
sampledf <- data.frame(sample_data(Ps.rare.LCMinf )) 
LCM_bray <- phyloseq::distance(Ps.rare.LCMinf , method = "bray")
Bray_clust <-hclust(LCM_bray)
adonis2(LCM_bray ~  ViralCopyNumberHighLow, data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~  Cleared.PersistentAtEndpoint, data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~  DiseaseCINEndpoint, data = sampledf, permutations = 9999)
adonis2(LCM_bray ~  InfectStageML , data = sampledf, permutations = 9999) # # 

Ps.rare.MOCK <-subset_samples(ps.E4ra.NoABX, Infected == "Mock")
sampledf <- data.frame(sample_data(Ps.rare.MOCK)) 
LCM_bray <- phyloseq::distance(Ps.rare.MOCK , method = "bray")
Bray_clust <-hclust(LCM_bray)
adonis2(LCM_bray ~  InfectStageML , data = sampledf, permutations = 9999) # # 


## Running Maaslin2 for determining differential abundance Supplemental Figures 2 and 3 
E4.relativeAbundanceOTHER_pivot <- read_csv("./E4relativeAbundanceMixomics_Genus.Pivot.csv")
input_genus_data <- as.data.frame(E4.relativeAbundanceOTHER_pivot)
rownames(input_genus_data) <-input_genus_data[,1]

input_meta_file <- E4.meta
input_meta_file <- subset.data.frame(input_meta_file, MouseID != "PBS")
input_meta_file <- subset.data.frame(input_meta_file, MouseID != "NegEX")
input_meta_file <- subset.data.frame(input_meta_file, ABX != "ABX")
input_meta_file$Infected <- factor(input_meta_file$Infected, levels = c("Mock", "Infected"))
VCHF.MF <- subset.data.frame(input_meta_file, ViralCopyNumberHighLow %in% c("Low", "High"))
VCP.MF <- subset.data.frame(input_meta_file, Cleared.PersistentAtEndpoint %in% c("Cleared", "Persistent"))
Inf.Stage.MF <- subset.data.frame(input_meta_file, Infected == "Infected")
Dis.CIN.MF <- subset.data.frame(input_meta_file, DiseaseCINEndpoint %in% c("Normal/Hyperplasia", "CIN3", "At least CIN3"))
Dis.CIN.MF2 <- subset.data.frame(input_meta_file, DiseaseCINEndpoint %in% c( "CIN3", "At least CIN3"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN/Mock.v.Infected", 
  fixed_effects = c("Infected"),
  reference = c("Infected", "Mock"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = VCHF.MF , 
  output = "./MAASLIN/ViralCopy.H.vs.L-2", 
  fixed_effects = c("ViralCopyNumberHighLow"),
  reference = c("ViralCopyNumberHighLow", "Low"),
  random_effects = c("MouseID"))

fit_data = Maaslin2(
  input_data = input_genus_data , 
  input_metadata = VCP.MF, 
  output = "./MAASLIN/Viral.Persist.v.Clear", 
  fixed_effects = c("Cleared.PersistentAtEndpoint"),
  reference = c("Cleared.PersistentAtEndpoint", "Cleared"),
  random_effects = c("MouseID"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Dis.CIN.MF, 
  output = "./MAASLIN/DiseaseCINEndpoint.v.CIN3", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "CIN3"),
  random_effects = c("MouseID"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Dis.CIN.MF2, 
  output = "./MAASLIN/DiseaseCINEndpoint.v.CIN3.vplus", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "CIN3"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Dis.CIN.MF, 
  output = "./MAASLIN/DiseaseCINEndpoint.v.Hyperplasia", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "Normal/Hyperplasia"),
  random_effects = c("MouseID"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Dis.CIN.MF, 
  output = "./MAASLIN/DiseaseCINEndpoint.v.CIN3plus", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "At Least CIN3"),
  random_effects = c("MouseID"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata =  Inf.Stage.MF, 
  output = "./MAASLIN/InfectedMice.InfectionStage_EarlyInfection.MAASLIN ", 
  fixed_effects = c("InfectStageML"),
  reference = c("InfectStageML", "Early Infection 5 & 7 wpi"),
  random_effects = c("MouseID"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata =  Inf.Stage.MF, 
  output = "./MAASLIN/InfectedMice.InfectionStage_MidInfection.MAASLIN ", 
  fixed_effects = c("InfectStageML"),
  reference = c("InfectStageML", "Mid Infection 9 & 14 wpi"),
  random_effects = c("MouseID"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata =  Inf.Stage.MF, 
  output = "./MAASLIN/InfectedMice.InfectionStage_LateInfection.MAASLIN ", 
  fixed_effects = c("InfectStageML"),
  reference = c("InfectStageML", "Late Infection 18 & 25 wpi"),
  random_effects = c("MouseID"))

## Mix omics for PLS-DA analysis --- at the GENUS level 
## the microbiome data - they highly recommend normalizing - 
#     i chose to import the relative abundance table (everything normalized to out of 1) 
#      of the whole data set so all samples were represented) then filtering out anything with an abundance
#      < 0.01 (1%) and calling it "other"-
ps.E4clean<- tax_glom(ps.E4 , "Genus")
ps.E4clean <- phyloseq_filter_prevalence(ps.E4clean, prev.trh = 0.03, abund.trh = 10, threshold_condition = "AND") # using AND here b/c need to ensure that the ASV is at least in 2 (>0.04 % of samples)
relE4.2 <- transform(ps.E4clean, "compositional") 
relE4.2.df <-psmelt(relE4.2)
relE4.2.df<- subset.data.frame(relE4.2.df, Abundance > 0.001)

write.csv(relE4.2.df, "./E4relativeAbundanceMixomics_Genus.csv")
E4.relativeAbundanceOTHER_pivot <- read_csv("./E4relativeAbundanceMixomics_Genus.Pivot.csv")
Micro.df <- as.data.frame(E4.relativeAbundanceOTHER_pivot) # 77 ASVs included 
rownames(Micro.df) <- Micro.df$SampleID
Micro.df <- Micro.df[,-1]
E4.meta.df <- data.frame(sample_data(ps.E4))

# PLSDA for ALL infected vs. Mock
E4.Meta<- subset.data.frame(E4.meta.df , Infected %in% c("Infected", "Mock"))
Micro.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E4.Meta))

X <- Micro.df
Y <- E4.Meta
Y.I <- Y$Infected
Y.I<-factor(Y.I, levels = c( "Mock","Infected"))

plsda.I <- plsda(X, Y.I, ncomp = 3)
plotIndiv(plsda.I, 
          ind.names = F, 
          ellipse = TRUE,
          comp = c(1,2),
          legend=TRUE, cex=6,
          col =c("#8AB17D","#053061"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'E4 Infection v. none')

cord1_Load<-plotLoadings(plsda.I, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Infected== "TRUE" & Contrib.Mock == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Mock", "Infected", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.I, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Mock== "TRUE" & Contrib.Infected == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Mock", "Infected", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))
Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.001)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.3)+
  scale_color_manual(values = c( "#053061","#8AB17D"))+
  scale_fill_manual(values = c("#053061","#8AB17D"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Mock", "Infected")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 3-Invection. v. mock. top driving Taxa')

# INFECTED mice only viral load H L 
E4.I.Meta<- subset.data.frame(E4.meta.df , Infected %in% c("Infected"))
Micro.I.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E4.I.Meta))

E4.I.dis.Meta<- subset.data.frame(E4.I.Meta , DiseaseCINEndpoint %in% c("Normal/Hyperplasia", "CIN3", "At least CIN3"))
Micro.I.dis.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E4.I.dis.Meta))

X.dis <-Micro.I.dis.df
Y.dis <- E4.I.dis.Meta$DiseaseCINEndpoint
Y.dis<-factor(Y.dis, levels = c("Normal/Hyperplasia", "CIN3", "At least CIN3", "SCC"))

E4.I.VHL.Meta<- subset.data.frame(E4.I.Meta , ViralCopyNumberHighLow %in% c("Low", "High"))
Micro.I.VHL.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E4.I.VHL.Meta))

X.VHL <-Micro.I.VHL.df
Y.VHL <- E4.I.VHL.Meta$ViralCopyNumberHighLow
Y.VHL <- factor(Y.VHL, levels = c("Low", "High" ))

E4.I.PC.Meta<- subset.data.frame(E4.I.Meta, Cleared.PersistentAtEndpoint %in% c("Cleared", "Persistent"))
Micro.I.PC.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E4.I.PC.Meta))

X.PC<-Micro.I.PC.df
Y.PC <- E4.I.PC.Meta$Cleared.PersistentAtEndpoint
Y.PC <- factor(Y.PC, levels = c("Cleared", "Persistent" ))

E4.I.ML.Meta<- subset.data.frame(E4.I.Meta, InfectStageML%in% c("Early Infection 5 & 7 wpi", "Mid Infection 9 & 14 wpi","Late Infection 18 & 25 wpi"))
Micro.I.ML.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E4.I.ML.Meta))

X.ML<-Micro.I.ML.df
Y.ML <- E4.I.ML.Meta$InfectStageML
Y.ML <- factor(Y.ML, levels = c("Early Infection 5 & 7 wpi", "Mid Infection 9 & 14 wpi","Late Infection 18 & 25 wpi"))

plsda.dis <- plsda(X.dis, Y.dis, ncomp = 3) # note no SCC in this group
plotIndiv(plsda.dis, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp= c(1,2),
          col =c("#8AB17D","#F6D27A","#F5A261"), 
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 3 infected mice Disease severity')

cord1_Load<-plotLoadings(plsda.dis, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Normal.Hyperplasia == "TRUE" & Contrib.CIN3 == "TRUE" & Contrib.At.least.CIN3 == "TRUE") # removing variables that contributed to ALL the groups
# "Normal/Hyperplasia", "CIN 3", "At least CIN3", "SCC"
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Normal.Hyperplasia", "CIN3", "At.least.CIN3",  "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.dis, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Normal.Hyperplasia == "TRUE" & Contrib.CIN3 == "TRUE" & Contrib.At.least.CIN3=="TRUE") # removing variables that contributed to ALL the groups
# removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Normal.Hyperplasia", "CIN3", "At.least.CIN3","GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.005)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(0,0.35)+
  scale_color_manual(values = c( "#F5A261","#F6D27A","#8AB17D"))+
  scale_fill_manual(values = c("#F5A261","#F6D27A","#8AB17D"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Normal/Hyperplasia", "CIN3", "At least CIN3")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 3 infected mice Disease severity, top driving Taxa')

E4.I.dis.Meta<- subset.data.frame(E4.I.Meta , DiseaseCINEndpoint %in% c("CIN3", "At least CIN3"))
Micro.I.dis.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E4.I.dis.Meta))

X.dis2 <-Micro.I.dis.df
Y.dis2 <- E4.I.dis.Meta$DiseaseCINEndpoint
Y.dis2<-factor(Y.dis2, levels = c("CIN3", "At least CIN3"))

plsda.dis2 <- plsda(X.dis2, Y.dis2, ncomp = 3) # note no SCC in this group
plotIndiv(plsda.dis2, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp= c(1,2),
          col =c("#F6D27A","#F5A261"), 
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'experiment 3 infected mice ONLY Disease severity')

cord1_Load<-plotLoadings(plsda.dis2, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.CIN3 == "TRUE" & Contrib.At.least.CIN3 == "TRUE") # removing variables that contributed to ALL the groups
# "Normal/Hyperplasia", "CIN 3", "At least CIN3", "SCC"
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable", "CIN3", "At.least.CIN3",  "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.dis2, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings,  Contrib.CIN3 == "TRUE" & Contrib.At.least.CIN3=="TRUE") # removing variables that contributed to ALL the groups
# removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable", "CIN3", "At.least.CIN3","GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.005)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(0,0.35)+
  scale_color_manual(values = c( "#F5A261","#F6D27A"))+
  scale_fill_manual(values = c("#F5A261","#F6D27A"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c( "CIN3", "At least CIN3")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 3 infected mice ONLY Disease severity, top driving Taxa')

plsda.VHL <- plsda(X.VHL, Y.VHL, ncomp = 3)
plotIndiv(plsda.VHL, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp = c(1,2),
          col =c(  "#2A9D8F", "#053061"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'E4 infected mice only Viral Load')

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

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.12)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.3)+
  scale_color_manual(values = c("#2A9D8F", "#053061"))+
  scale_fill_manual(values = c( "#2A9D8F", "#053061"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Low", "High")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 3 infected mice - Viral Load, top driving Taxa')


plsda.PC <- plsda(X.PC, Y.PC, ncomp = 3)
plotIndiv(plsda.PC, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp = c(1,2),
          col =c(  "#2A9D8F", "#8AB17D"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 3 infected mice only Viral persist and clearance')

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

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.1)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .03)+
  xlim(-0,0.5)+
  scale_color_manual(values = c("#8AB17D","#2A9D8F"))+
  scale_fill_manual(values = c("#8AB17D", "#2A9D8F"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Cleared", "Persistent")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 3 infected mice - Viral persistance and clearance, top driving Taxa')


plsda.ML <- plsda(X.ML, Y.ML, ncomp = 3) # note no SCC in this group
plotIndiv(plsda.ML, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp= c(1,3),
          #pch = c(3,4,5,2,1),
          col =c("#186275","#021327", "#053061" ), #"#DE7350"), 
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 3 - early, mid, late infection stage - infectedd mice only')

cord1_Load <-plotLoadings(plsda.ML, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings,   Contrib.Early.Infection.5...7.wpi=="TRUE"& Contrib.Mid.Infection.9...14.wpi == "TRUE" & Contrib.Late.Infection.18...25.wpi=="TRUE") # removing variables that contributed to ALL the groups
# "Normal/Hyperplasia", "CIN 3", "At least CIN3", "SCC"
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable", "Contrib.Early.Infection.5...7.wpi", "Contrib.Mid.Infection.9...14.wpi","Contrib.Mid.Infection.9...14.wpi",  "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.ML, comp = 3, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings,  Contrib.Early.Infection.5...7.wpi=="TRUE"& Contrib.Mid.Infection.9...14.wpi == "TRUE" & Contrib.Mid.Infection.9...14.wpi=="TRUE") # removing variables that contributed to ALL the groups
# removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Contrib.Early.Infection.5...7.wpi", "Contrib.Mid.Infection.9...14.wpi","Contrib.Mid.Infection.9...14.wpi",  "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.07)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(0,0.4)+
  scale_color_manual(values = c("#186275","#021327", "#053061" ))+
  scale_fill_manual(values = c("#186275","#021327", "#053061" ))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Early Infection 5 & 7 wpi", "Mid Infection 9 & 14 wpi","Late Infection 18 & 25 wpi")), scales = "free_y", nrow =1)+
  theme_bw()+
  ggtitle('Experiment 3 Infection Stage- top driving Taxa')
