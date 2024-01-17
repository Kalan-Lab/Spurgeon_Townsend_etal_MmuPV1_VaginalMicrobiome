# Vaginal Microbiome 
# Experiment 5 - laser capture microdisection sample Microbiome Analysis 

library(ggplot2)
library(ggrepel)
library(microbiome)
library(tidyverse)
library(dplyr)
library(broom)
library(qiime2R)
library(xlsx)
library(RColorBrewer)
library(vegan)
library(phyloseq)
library(magrittr)
library(Maaslin2)
library(mixOmics)
library(metagMisc)

setwd("/Users/liztown/Documents/KalanLab/Papers/Vaginal Microbiome/Rcode_forManuscript/Experiment5_LCM")

ps.LCMre <- qza_to_phyloseq(features = "./table-LCMre-dada2.qza",
                            taxonomy = "./taxonomy-LCMre.qza",
                            metadata = "./LCMreSeqManufest.tsv",
                            tree = "./LCMre-rooted-tree.qza")
LCM.meta <- data.frame(sample_data(ps.LCMre))
LCM.meta$SampleID <- row.names(LCM.meta)

# cleaning the data set 
ps.LCMclean<- subset_taxa(ps.LCMre, Family !="mitochondria")
ps.LCMclean<- subset_taxa(ps.LCMclean, Phylum !="Chloroflexi")
ps.LCMclean<- subset_taxa(ps.LCMclean, Class !="Chloroplast")
ps.LCMclean<- subset_taxa(ps.LCMclean, Phylum !="")
ps.LCMclean<- subset_taxa(ps.LCMclean, Genus !="Methylobacterium")
ps.LCMclean<- subset_taxa(ps.LCMclean, Genus !="Sphingomonas")
ps.LCMclean <- subset_taxa(ps.LCMclean, Genus !="Ralstonia") 
ps.LCMclean <- subset_taxa(ps.LCMclean, Genus !="Facklamia") 
ps.LCMclean <- subset_taxa(ps.LCMclean, Genus !="Chryseobacterium") 
ps.LCMclean <- subset_taxa(ps.LCMclean, Genus !="Curvibacter") 
ps.LCMclean <- subset_taxa(ps.LCMclean, taxa(ps.LCMclean@otu_table) !="836ac518d86fdb0dd587761ae39b8845") # removing a corynebacterium stratium sequence that is ~70 percent of the reads in 2 of the samples. probable contaminant 
ps.LCMclean <- subset_samples(ps.LCMclean, sample_names(ps.LCMclean) != "LK16S004-305") # removing negative controls 
ps.LCMclean <- subset_samples(ps.LCMclean, sample_names(ps.LCMclean) != "LK16S004-306")
ps.LCMclean <- phyloseq_filter_prevalence(ps.LCMclean, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR") 
ps.LCMclean

#reordering by treatment 
ps.LCMclean@sam_data$Treatment <- factor(ps.LCMclean@sam_data$Treatment, levels = c("Mock", "MmuPV1 only", "Mock + 6m E2", "MmuPV1 + 6m E2"))

# Output files
relLCM <- transform(ps.LCMclean, "compositional") # transform is the function to make it relative abundance
relLCM.df <-psmelt(relLCM)
LCMout <- psmelt(ps.LCMclean)
#write.csv(LCMout, "./LCMtotalAbundanceOut.csv")
LCMrelOut <-psmelt(relLCM)
#write.csv(LCMrelOut, "./LCMrelativeAbundanceOut.csv")

### Experiment 5  RELATIVE ABUNDANCE PLOTS - Figure 4
relLCM.df <-psmelt(relLCM)
OtherLCM <- relLCM.df# the dataframe of the relative abundance within each sample
OtherLCM <- subset.data.frame(OtherLCM, Abundance > 0.02)
OtherLCM <- OtherLCM %>% mutate(Genus = ifelse(Abundance < 0.05, "Other", Genus)) # replacing the genus catagory with "Other" if the abundance is less than 1% 
OtherLCM <- OtherLCM %>% mutate(Phylum = ifelse(Abundance < 0.05, "Other", Phylum)) 
OtherLCM$Treatment <- factor(OtherLCM$Treatment, levels = c("Mock", "MmuPV1 only","Mock + 6m E2", "MmuPV1 + 6m E2"))
table(OtherLCM$Genus, OtherLCM$Phylum)

OtherLCM$Phylum <- factor(OtherLCM$Phylum, levels = c("Other", "Actinobacteria",  "Firmicutes", "Bacteroidetes","Fusobacteria", "Nitrospirae", "Proteobacteria", "Tenericutes", "Verrucomicrobia"))
OtherLCM$Genus<- factor(OtherLCM$Genus, levels = c("Other",
                                                   "Bifidobacterium","Brachybacterium","Corynebacterium","Propionibacterium", "Dermacoccus", "Friedmanniella", "Gardnerella", "Kocuria",
                                                   "Microbacterium", "Micrococcus","Mycobacterium","Pseudonocardia", "Rubrobacter","Solirubrobacter", 
                                                   "Streptomyces", ##Actinobacteria / actinomycetota (15)
                                                   
                                                   "Anaerobacillus",  "Carnobacterium", "Clostridium" , "Enterococcus", "Faecalibacterium", "Gemella","Lactobacillus", 
                                                   "Roseburia", "Staphylococcus", "Streptococcus","Vagococcus", # Firmicutes/ Bacillota (11)
                                                   
                                                   "Alistipes","Bacteroides",  "Cloacibacterium","Cytophaga","Flavobacterium", "Porphyromonas","Prevotella", "Sediminibacterium",
                                                   "Sphingobacterium","Spirosoma", "Wautersiella", # Bacteroidetes / Bacteroidota (11)
                                                   "Leptotrichia", # fusobacteria
                                                   "Nitrospira", # Nitrospirae
                                                   
                                                   "Acinetobacter","Actinobacillus","Aeromonas", "Aggregatibacter","Brachymonas","Burkholderia", "Cellvibrio", "Enhydrobacter",
                                                   "Massilia","Novosphingobium", "Ochrobactrum","Paracoccus","Pedomicrobium","Petrobacter",  "Pseudomonas", "Rubellimicrobium","Skermanella",  
                                                   "Sphingobium", "Stenotrophomonas",#Proteobacteria /pseudomonadota (19)
                                                   
                                                   "Anaeroplasma", #Tenericutes
                                                   "Akkermansia","DA101")) #Verrucomicrobia
BigMicro2 <-c("#C1C1C1", #Gray
              colorRampPalette(c("#031A34", "#053061", "#21777C", "#2A9D8F", "#8AB17D", "#C0C27B","#D2C77A", "#F5D279","#EFB366","#F4A261","#DE7350","#C8443F","#913C44", "#733036"))(61))

ggplot(OtherLCM, aes(x= Treatment, y =Abundance, fill = Genus, color = Genus))+
  geom_bar(stat="identity", position="fill") +
  facet_wrap("Treatment", scales = "free_x", nrow = 1)+
  ylab("Relative Abundance")+
  theme_light()+
  scale_fill_manual(values = BigMicro2)+
  scale_color_manual(values = BigMicro2)+
  theme(axis.text.x = element_text(size = 4)) +
  ggtitle("Genus present > 0.5% of reads in a sample")


# Alpha diversity (not in paper)
plot_richness(ps.LCMclean, x = "Treatment", measures=c("Shannon")) + 
  geom_boxplot(fill = "gray")  + 
  theme(text = element_text(size = 18))+
  #facet_wrap("Treatment", scales = "free_x") + 
  theme_bw()
ps.LCMclean@sam_data$Disease2 <- factor(ps.LCMclean@sam_data$Disease2, levels = c("No disease", "Hyperplasia (No Dysplasia)", "Low CIN", "High CIN", "SCC")) 
ps.LCMclean@sam_data$Estrus.Stage <- factor(ps.LCMclean@sam_data$Estrus.Stage, levels = c("Proestrus", "Estrus", "Estrus/Metestrus", "Metestrus ", "Metestrus/Diestrus", "Diestrus", "Diestrus/Proestrus")) 
tab <-microbiome::alpha(ps.LCMclean, index = "all")
#write.csv(tab, "./LCMoutAlphaTable.csv")
AlphaLCM <- tab
AlphaLCM$SampleID<- row.names(AlphaLCM)
AlphaLCM <- merge(AlphaLCM, LCM.meta, by = "SampleID", all = F) # the alpha table above but with the metadata added back in 
AlphaLCM$Treatment <- factor(AlphaLCM$Treatment, levels = c("Mock", "Mock + 6m E2","MmuPV1 only", "MmuPV1 + 6m E2"))
AlphaLCM$Disease2 <- factor(AlphaLCM$Disease2, levels = c("No disease", "Hyperplasia (No Dysplasia)", "Low CIN", "High CIN", "SCC")) 
AlphaLCM$Estrus.Stage <- factor(AlphaLCM$Estrus.Stage, levels = c("Proestrus", "Estrus", "Estrus/Metestrus", "Metestrus ", "Metestrus/Diestrus", "Diestrus", "Diestrus/Proestrus")) 

ggplot(AlphaLCM, aes(x = Treatment, y =diversity_shannon, color = Treatment))+
  geom_point(size = 5)+
  stat_summary(fun.y = mean,  geom = "point",shape = 95, size = 13, show.legend = F)+ 
  scale_color_manual(values = c("#21777C","#8AB17D","#EFB366","#DE7350"))+
  facet_wrap("Treatment", scales = "free_x", ncol = 4) + 
  #ylim(0,6)+
  theme_light()+
  ggtitle("Shannon Alpha Diversity LCM")


# Beta Diversity Figure 4 C
min(sample_sums(ps.LCMclean)) #383
median(sample_sums(ps.LCMclean)) #3805
max(sample_sums(ps.LCMclean)) #33632

rarecurve(t(otu_table(ps.LCMclean)), step=20, xlim= c(0,1000)) ## considering using a read cut off od 500 for beta diversity metrics 

ps.rarefied = rarefy_even_depth(ps.LCMclean, rngseed=1, sample.size=800, replace=F) # THIS IS THE ONE TO GO WITH
LCMreRareified <-psmelt(ps.rarefied)

GP = ps.rarefied 
GP.ord <- ordinate(GP, "NMDS", "bray") # For Bray Curtis of the rareified dataset
plot_ordination(GP, GP.ord, type="samples", color="Treatment") + 
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(size=5, ) + ggtitle("LCM - Bray Curtis ")+
  scale_color_manual(values = c("#DE7350","#EFB366","#8AB17D", "#21777C"))+
  theme_bw()



# Type 2 Permanovas for Table 2
sampledf <- data.frame(sample_data(ps.rarefied)) 
LCM_bray <- phyloseq::distance(ps.rarefied, method = "bray")
Bray_clust <-hclust(LCM_bray)

adonis2(LCM_bray ~ Infected + E2 + Estrus.Stage  , by= "margin", data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~ Infected   , by= "margin", data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~  E2  , by= "margin", data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~  Estrus.Stage  , by= "margin", data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~  ViralCopyNumberHighLow, data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~  Disease2, data = sampledf, permutations = 9999) # 

Ps.rare.LCMinf <-subset_samples(ps.rarefied, Infected == "MmuPV")
sampledf <- data.frame(sample_data(Ps.rare.LCMinf )) 
LCM_bray <- phyloseq::distance(Ps.rare.LCMinf , method = "bray")
Bray_clust <-hclust(LCM_bray)
adonis2(LCM_bray ~ E2 + Estrus.Stage  , by= "margin", data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~  ViralCopyNumberHighLow, data = sampledf, permutations = 9999) # 
adonis2(LCM_bray ~  Disease2, data = sampledf, permutations = 9999) # 




## Running Maaslin2 for determining differential abundance - Supplemental figure 5 
# (and a few analysises that did not make it into the final figures)
# https://huttenhower.sph.harvard.edu/maaslin/
LCMrelativeAbundance_Genus <- read.csv("./MAASLIN2/LCMrelativeAbundance_Genus.csv")
input_genus_data <- as.data.frame(LCMrelativeAbundance_Genus)
rownames(input_genus_data) <-input_genus_data[,1]

input_meta_file <- LCM.meta
input_meta_file$Disease_Levels = factor(input_meta_file$Disease2, levels = c("No disease", "Hyperplasia (No Dysplasia)", "Low CIN", "High CIN", "SCC"))
input_meta_file$Treatment = factor(input_meta_file$Treatment, levels = levels = c("Mock", "MmuPV1 only","Mock + 6m E2", "MmuPV1 + 6m E2"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2/LCM_TreatmentGroup_Genus", 
  fixed_effects = c("Treatment"),
  reference = c("Treatment", "Mock"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2/LCM_DiseaseSeverity_allGroups_Genus", 
  fixed_effects = c("Disease2"),
  reference = c("Disease2", "No Disease"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2/LCM_DiseaseSeveritySCC_allGroups_Genus", 
  fixed_effects = c("Disease2"),
  reference = c("Disease2", "SCC"))

NoEstrogen<-subset.data.frame(input_meta_file, E2 == "none")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = NoEstrogen, 
  output = "./MAASLIN2/LCM_Infection_inNoEstrogenGroup_genus", 
  fixed_effects = c("Infected"),
  reference = c("Infected","Mock"))

Estrogen<-subset.data.frame(input_meta_file, E2 == "E2")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Estrogen, 
  output = "./MAASLIN2/LCM_Infection_inEstrogenGroup_genus", 
  fixed_effects = c("Infected"),
  reference = c("Infected","Mock"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = input_meta_file, 
  output = "./MAASLIN2/Allmice_EstrogenV.none_try2_genus", 
  fixed_effects = c("E2"),
  reference = c("E2", "None"))

Mock<-subset.data.frame(input_meta_file, Infected == "Mock")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Mock, 
  output = "./MAASLIN2/MockMice_EstrogenV.none_try2_genus", 
  fixed_effects = c("E2"),
  reference = c("E2", "None"))

Infected<-subset.data.frame(input_meta_file, Infected == "MmuPV")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Infected, 
  output = "./MAASLIN2/InfectedMice_EstrogenV.none_try2_genus", 
  fixed_effects = c("E2"),
  reference = c("E2", "None"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Infected, 
  output = "./MAASLIN2/LCM_DiseaseSeverity_inInfectionGroup_genus", 
  fixed_effects = c("Disease2"),
  reference = c("Disease2","Low CIN"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Infected, 
  output = "./MAASLIN2/LCM_DiseaseSeveritySCC_inInfectionGroup_genus", 
  fixed_effects = c("Disease2"),
  reference = c("Disease2","SCC"))


## MixOmics for PLS-DA analysis - at the genera level -for Figure 4 and Supplemental Figure 5
## the microbiome data - they highly recommend normalizing - 
#     i chose to import the relative abundance table (everything normalized to out of 1) 
#      of the whole data set so all samples were represented) then filtering out anything with an abundance
#      < 0.01 (1%) and calling it "other"-
## Note - all the dataframes need the same number of rows. (ie 10 rows - one per subjecr in subject order)
ps.ps.LCMclean2<- tax_glom(ps.LCMclean, "Genus")
ps.LCMclean2 <-metagMisc:: phyloseq_filter_prevalence(ps.LCMclean2, prev.trh = 0.05, abund.trh = 10, threshold_condition = "AND") # using AND here b/c need to ensure that the ASV is at least in 2 (>0.04 % of samples)
relLCM2 <- microbiome::transform(ps.LCMclean2, "compositional") 
relLCM2.df <-psmelt(relLCM2)
OtherLCM2 <- relLCM2.df# the dataframe of the relative abundance within each sample
OtherLCM2 <- subset.data.frame(OtherLCM2, Abundance > 0.001)

#write.csv(OtherLCM2, "./Mixomics.Genus.PLSDA/LCMrelativeAbundanceOTHER.Genus.csv")
LCMrelativeAbundanceOTHER_pivot <- read_csv("./Mixomics.Genus.PLSDA/LCMrelativeAbundanceOTHER.Genus.pivot.csv")

LCM.Micro.df <- as.data.frame(LCMrelativeAbundanceOTHER_pivot) # 77 ASVs included 
rownames(LCM.Micro.df) <- LCM.Micro.df$SampleID
LCM.Micro.df <- LCM.Micro.df[,-1]

LCM.meta <- data.frame(sample_data(ps.LCMre))
LCM.meta$SampleID <- row.names(LCM.meta)

# PLSDA for ALL infected 
LCM.Meta.Infect<- subset.data.frame(LCM.meta, Infected == "MmuPV")
LCM.inf.Micro.df <- subset.data.frame(LCM.Micro.df, row.names(LCM.Micro.df) %in% row.names(LCM.Meta.Infect))

X <- LCM.inf.Micro.df
Y <- LCM.Meta.Infect
Y.dis <- Y$Disease2
Y.dis<-factor(Y.dis, levels = c("Low CIN", "High CIN", "SCC"))
Y.VHL <- Y$ViralCopyNumberHighLow
Y.VHL <- factor(Y.VHL, levels = c("Low", "High" ))
Y.PC <- Y$Cleared.Persistent # not assessed here - all mice had persistent infection 
Y.E2 <- Y$E2
Y.E2 <- factor(Y.E2, levels = c("none", "E2"))

plsda.dis <- plsda(X, Y.dis, ncomp = 2)
plotIndiv(plsda.dis, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          col =c("#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'LCM - infected mice only Disease severity')

cord1_Load<-plotLoadings(plsda.dis, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% mutate(VariateAxis = "Variate Axis 1")
cord2_Load<-plotLoadings(plsda.dis, comp = 2, contrib = 'max', method = 'median')
cord2_Load <-cord2_Load %>% mutate(VariateAxis = "Variate Axis 2")
Cord_Loadings_dataframe <- rbind(cord1_Load , cord2_Load)

Cord_Loadings_dataframe <- tibble::rownames_to_column(Cord_Loadings_dataframe, "Variable")
TiesX3 <- filter(Cord_Loadings_dataframe, Contrib.Low.CIN == "TRUE" & Contrib.High.CIN == "TRUE" & Contrib.SCC == "TRUE") # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_Loadings_dataframe_filt<-Cord_Loadings_dataframe[ !Cord_Loadings_dataframe$Variable %in% Tie3Var,]

Low.CIN <- filter(Cord_Loadings_dataframe_filt, Contrib.Low.CIN == "TRUE")
Low.CIN$GroupContrib[Low.CIN$GroupContrib == "tie"] <- "LowCIN"
High.CIN <- filter(Cord_Loadings_dataframe_filt, Contrib.High.CIN == "TRUE")
High.CIN$GroupContrib[High.CIN$GroupContrib == "tie"] <- "HighCIN"
SCC<- filter(Cord_Loadings_dataframe_filt, Contrib.SCC == "TRUE")
SCC$GroupContrib[SCC$GroupContrib == "tie"] <- "SCC"

Cord_Filtered_Loadings.DF <- rbind(Low.CIN,High.CIN,SCC)

TopCord <-subset.data.frame(Cord_Filtered_Loadings.DF, abs(importance) > 0.1)
ggplot(TopCord, aes(x = importance, y = reorder(Variable, abs(importance)), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*abs(importance)))+
  geom_col(width = .04)+
  xlim(-.3,0.3)+
  scale_color_manual(values = c("#F6D27A", "#8AB17D","#F5A261"))+
  scale_fill_manual(values = c("#F6D27A", "#8AB17D","#F5A261"))+
  facet_wrap(~VariateAxis, scales = "free_y")+
  #facet_wrap(~factor(GroupContrib, levels = c("Low CIN", "High CIN", "SCC")), scales = "free_y")+
  theme_bw()+
  ggtitle('LCM - infected mice only Disease severity, top driving Taxa')

ggplot(TopCord, aes(x = abs(importance), y = reorder(Variable, abs(importance)), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*abs(importance)))+
  geom_col(width = .04)+
  xlim(0,0.3)+
  scale_color_manual(values = c("#F6D27A", "#8AB17D","#F5A261"))+
  scale_fill_manual(values = c("#F6D27A", "#8AB17D","#F5A261"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Low CIN", "High CIN", "SCC")), scales = "free_y")+
  theme_bw()+
  ggtitle('LCM - infected mice only Disease severity, top driving Taxa')

cord1_Load<-plotLoadings(plsda.dis, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Low.CIN == "TRUE" & Contrib.High.CIN == "TRUE" & Contrib.SCC == "TRUE") # removing variables that contributed to ALL the groups
# removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Low.CIN", "High.CIN", "SCC", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.dis, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Low.CIN == "TRUE" & Contrib.High.CIN == "TRUE" & Contrib.SCC == "TRUE") # removing variables that contributed to ALL the groups
# removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Low.CIN", "High.CIN", "SCC", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.05)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(0,0.3)+
  scale_color_manual(values = c("#F6D27A", "#8AB17D","#F5A261"))+
  scale_fill_manual(values = c("#F6D27A", "#8AB17D","#F5A261"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Low CIN", "High CIN", "SCC")), scales = "free_y")+
  theme_bw()+
  ggtitle('LCM - infected mice only Disease severity, top driving Taxa')


plsda.VHL <- plsda(X, Y.VHL, ncomp = 2)
plotIndiv(plsda.VHL, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          col =c(  "#2A9D8F", "#053061"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'LCM - infected mice only Viral Load')


cord1_Load<-plotLoadings(plsda.VHL, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% mutate(VariateAxis = "Variate Axis 1")
cord2_Load<-plotLoadings(plsda.VHL, comp = 2, contrib = 'max', method = 'median')
cord2_Load <-cord2_Load %>% mutate(VariateAxis = "Variate Axis 2")
Cord_Loadings_dataframe <- rbind(cord1_Load , cord2_Load)

Cord_Loadings_dataframe <- tibble::rownames_to_column(Cord_Loadings_dataframe, "Variable")
TiesX3 <- filter(Cord_Loadings_dataframe, Contrib.Low== "TRUE" & Contrib.High == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_Loadings_dataframe_filt.df<-Cord_Loadings_dataframe[ !Cord_Loadings_dataframe$Variable %in% Tie3Var,]

TopCord <-subset.data.frame(Cord_Loadings_dataframe_filt.df, abs(importance) > 0.05)
ggplot(TopCord, aes(x = importance, y = reorder(Variable, abs(importance)), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*abs(importance)))+
  geom_col(width = .04)+
  xlim(-.3,0.3)+
  scale_color_manual(values = c("#2A9D8F", "#053061"))+
  scale_fill_manual(values = c("#2A9D8F", "#053061"))+
  facet_wrap(~VariateAxis, scales = "free_y")+
  #facet_wrap(~factor(GroupContrib, levels = c("Low CIN", "High CIN", "SCC")), scales = "free_y")+
  theme_bw()+
  ggtitle('LCM - infected mice only Viral Load, top driving Taxa')


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

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.01)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.3)+
  scale_color_manual(values = c("#053061", "#2A9D8F"))+
  scale_fill_manual(values = c("#053061", "#2A9D8F"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Low", "High")), scales = "free_y")+
  theme_bw()+
  ggtitle('LCM - infected mice only Viral Load, top driving Taxa')

# 
plsda.In.E2 <- plsda(X, Y.E2, ncomp = 3)
plotIndiv(plsda.In.E2, 
          ind.names = F, 
          comp = c(1,2),
          ellipse = TRUE,
          legend=TRUE, cex=6,
          col =c("#F6D27A","#F5A261"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'LCM - infected mice only Estrogen')


cord1_Load<-plotLoadings(plsda.In.E2, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% mutate(VariateAxis = "Variate Axis 1")
cord2_Load<-plotLoadings(plsda.In.E2, comp = 2, contrib = 'max', method = 'median')
cord2_Load <-cord2_Load %>% mutate(VariateAxis = "Variate Axis 2")
Cord_Loadings_dataframe <- rbind(cord1_Load , cord2_Load)

Cord_Loadings_dataframe <- tibble::rownames_to_column(Cord_Loadings_dataframe, "Variable")
TiesX3 <- filter(Cord_Loadings_dataframe, Contrib.none== "TRUE" & Contrib.E2 == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_Loadings_dataframe_filt.df<-Cord_Loadings_dataframe[ !Cord_Loadings_dataframe$Variable %in% Tie3Var,]

TopCord <-subset.data.frame(Cord_Loadings_dataframe_filt.df, abs(importance) > 0.05)
ggplot(TopCord, aes(x = importance, y = reorder(Variable, abs(importance)), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*abs(importance)))+
  geom_col(width = .04)+
  xlim(-.3,0.3)+
  scale_color_manual(values = c("#F6D27A","#F5A261"))+
  scale_fill_manual(values = c("#F6D27A","#F5A261"))+
  facet_wrap(~VariateAxis, scales = "free_y")+
  #facet_wrap(~factor(GroupContrib, levels = c("Low CIN", "High CIN", "SCC")), scales = "free_y")+
  theme_bw()+
  ggtitle('LCM - infected mice only Viral Load, top driving Taxa')


cord1_Load<-plotLoadings(plsda.In.E2, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("importance" ="Variate.1.importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.none== "TRUE" & Contrib.E2 == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","none", "E2", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.In.E2, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename( "importance" = "Variate.2.importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.none== "TRUE" & Contrib.E2 == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","none", "E2", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.2)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.42)+
  scale_color_manual(values = c("#F6D27A","#F5A261"))+
  scale_fill_manual(values = c("#F6D27A","#F5A261"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("none", "E2")), scales = "free_y")+
  theme_bw()+
  ggtitle('LCM - infected mice only Estrogen v.none, top driving Taxa')

# PLSDA for ALL infected W/ estrogen
LCM.Meta.Inf.E2 <- subset.data.frame(LCM.Meta.Infect, E2 == "E2")
LCM.E2.inf.Micro.df <- subset.data.frame(LCM.Micro.df, row.names(LCM.Micro.df) %in% row.names(LCM.Meta.Inf.E2))

X <- LCM.E2.inf.Micro.df
Y <- LCM.Meta.Inf.E2
Y.dis <- Y$Disease2
Y.dis<-factor(Y.dis, levels = c("Low CIN", "High CIN", "SCC"))
Y.VHL <- Y$ViralCopyNumberHighLow # Only one had Low viral load - so unable to assess
Y.VHL <- factor(Y.VHL, levels = c("Low", "High" ))
Y.PC <- Y$Cleared.Persistent # not assessed here - all mice had persistent infection 

plsda.dis <- plsda(X, Y.dis, ncomp = 2)
plotIndiv(plsda.dis, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          col =c("#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'LCM - infected mice + Estrogen only Disease severity')

cord1_Load<-plotLoadings(plsda.dis, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% mutate(VariateAxis = "Variate Axis 1")
cord2_Load<-plotLoadings(plsda.dis, comp = 2, contrib = 'max', method = 'median')
cord2_Load <-cord2_Load %>% mutate(VariateAxis = "Variate Axis 2")
Cord_Loadings_dataframe <- rbind(cord1_Load , cord2_Load)

Cord_Loadings_dataframe <- tibble::rownames_to_column(Cord_Loadings_dataframe, "Variable")
TiesX3 <- filter(Cord_Loadings_dataframe, Contrib.Low.CIN == "TRUE" & Contrib.High.CIN == "TRUE" & Contrib.SCC == "TRUE") # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_Loadings_dataframe_filt<-Cord_Loadings_dataframe[ !Cord_Loadings_dataframe$Variable %in% Tie3Var,]

Low.CIN <- filter(Cord_Loadings_dataframe_filt, Contrib.Low.CIN == "TRUE")
Low.CIN$GroupContrib[Low.CIN$GroupContrib == "tie"] <- "LowCIN"
High.CIN <- filter(Cord_Loadings_dataframe_filt, Contrib.High.CIN == "TRUE")
High.CIN$GroupContrib[High.CIN$GroupContrib == "tie"] <- "HighCIN"
SCC<- filter(Cord_Loadings_dataframe_filt, Contrib.SCC == "TRUE")
SCC$GroupContrib[SCC$GroupContrib == "tie"] <- "SCC"

Cord_Filtered_Loadings.DF <- rbind(Low.CIN,High.CIN,SCC)

TopCord <-subset.data.frame(Cord_Filtered_Loadings.DF, abs(importance) > 0.1)
ggplot(TopCord, aes(x = importance, y = reorder(Variable, abs(importance)), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*abs(importance)))+
  geom_col(width = .04)+
  xlim(-.3,0.3)+
  scale_color_manual(values = c("#F6D27A", "#8AB17D","#F5A261"))+
  scale_fill_manual(values = c("#F6D27A", "#8AB17D","#F5A261"))+
  facet_wrap(~VariateAxis, scales = "free_y")+
  #facet_wrap(~factor(GroupContrib, levels = c("Low CIN", "High CIN", "SCC")), scales = "free_y")+
  theme_bw()+
  ggtitle('LCM - infected mice +E2 only Disease severity, top driving Taxa')

ggplot(TopCord, aes(x = abs(importance), y = reorder(Variable, abs(importance)), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*abs(importance)))+
  geom_col(width = .04)+
  xlim(0,0.3)+
  scale_color_manual(values = c("#F6D27A", "#8AB17D","#F5A261"))+
  scale_fill_manual(values = c("#F6D27A", "#8AB17D","#F5A261"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Low CIN", "High CIN", "SCC")), scales = "free_y")+
  theme_bw()+
  ggtitle('LCM - infected mice only Disease severity, top driving Taxa')

cord1_Load<-plotLoadings(plsda.dis, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Low.CIN == "TRUE" & Contrib.High.CIN == "TRUE" & Contrib.SCC == "TRUE") # removing variables that contributed to ALL the groups
# removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Low.CIN", "High.CIN", "SCC", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.dis, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Low.CIN == "TRUE" & Contrib.High.CIN == "TRUE" & Contrib.SCC == "TRUE") # removing variables that contributed to ALL the groups
# removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Low.CIN", "High.CIN", "SCC", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.09)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(0,0.3)+
  scale_color_manual(values = c("#F6D27A", "#8AB17D","#F5A261"))+
  scale_fill_manual(values = c("#F6D27A", "#8AB17D","#F5A261"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Low CIN", "High CIN", "SCC")), scales = "free_y")+
  theme_bw()+
  ggtitle('LCM - infected mice with Estrogen only Disease severity, top driving Taxa')

### for MOCK mice (estrogen v. none )
LCM.Meta.Mock<- subset.data.frame(LCM.meta, Infected == "Mock")
LCM.Mock.Micro.df <- subset.data.frame(LCM.Micro.df, row.names(LCM.Micro.df) %in% row.names(LCM.Meta.Mock))

X <- LCM.Mock.Micro.df
Y <- LCM.Meta.Mock
Y.E2 <- Y$E2
Y.E2 <- factor(Y.E2, levels = c("none", "E2"))


plsda.M.E2 <- plsda(X, Y.E2, ncomp = 3)
plotIndiv(plsda.M.E2, 
          ind.names = F, 
          comp = c(1,2),
          ellipse = TRUE,
          legend=TRUE, cex=6,
          col =c("#F6D27A","#F5A261"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'LCM - Mock mice only Estrogen')


cord1_Load<-plotLoadings(plsda.M.E2, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% mutate(VariateAxis = "Variate Axis 1")
cord2_Load<-plotLoadings(plsda.M.E2, comp = 2, contrib = 'max', method = 'median')
cord2_Load <-cord2_Load %>% mutate(VariateAxis = "Variate Axis 2")
Cord_Loadings_dataframe <- rbind(cord1_Load , cord2_Load)

Cord_Loadings_dataframe <- tibble::rownames_to_column(Cord_Loadings_dataframe, "Variable")
TiesX3 <- filter(Cord_Loadings_dataframe, Contrib.none== "TRUE" & Contrib.E2 == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_Loadings_dataframe_filt.df<-Cord_Loadings_dataframe[ !Cord_Loadings_dataframe$Variable %in% Tie3Var,]

TopCord <-subset.data.frame(Cord_Loadings_dataframe_filt.df, abs(importance) > 0.05)
ggplot(TopCord, aes(x = importance, y = reorder(Variable, abs(importance)), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*abs(importance)))+
  geom_col(width = .04)+
  xlim(-.3,0.3)+
  scale_color_manual(values = c("#F6D27A","#F5A261"))+
  scale_fill_manual(values = c("#F6D27A","#F5A261"))+
  facet_wrap(~VariateAxis, scales = "free_y")+
  #facet_wrap(~factor(GroupContrib, levels = c("Low CIN", "High CIN", "SCC")), scales = "free_y")+
  theme_bw()+
  ggtitle('LCM - Mock mice only Viral Load, top driving Taxa')


cord1_Load<-plotLoadings(plsda.M.E2, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("importance" ="Variate.1.importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.none== "TRUE" & Contrib.E2 == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","none", "E2", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.M.E2, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename( "importance" = "Variate.2.importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.none== "TRUE" & Contrib.E2 == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","none", "E2", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.01)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.4)+
  scale_color_manual(values = c("#F6D27A","#F5A261"))+
  scale_fill_manual(values = c("#F6D27A","#F5A261"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("none", "E2")), scales = "free_y")+
  theme_bw()+
  ggtitle('LCM - Mock mice only Estrogen v.none, top driving Taxa')


