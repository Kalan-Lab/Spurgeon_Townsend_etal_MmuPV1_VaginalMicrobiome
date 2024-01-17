# Pilot Experiment E2 - cleaned 

library(tidyverse)
library(broom)
library(ggrepel)
library(phyloseq)
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
library(pheatmap)
library(metagMisc)

setwd("/Users/liztown/Documents/KalanLab/Papers/Vaginal Microbiome/Rcode_forManuscript/Experiment1_Pilot")

ps.E2pilot <- qza_to_phyloseq(features = "./table-Mouse-dada2-20-160.qza",
                              taxonomy = "./taxonomy-Mouse.qza",
                              metadata = "./E2PilotMouseManufest.txt",
                              tree = "./MousePilot-rooted-tree.qza")

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
ps.E2pilotDe <- subset_samples(ps.E2pilotDe, SampleOrigin == "Mouse")  # removing the negative control
ps.E2pilotDe# THE CLEANED FILE 

OutCV <- psmelt(ps.E2pilotDe)
#write.csv(OutCV, "./E2PilotAbundanceTotalsOut.csv")
ps.PilotREl <- transform(ps.E2pilotDe, "compositional")
PilotRelout <- psmelt(ps.PilotREl)
#write.csv(PilotRelout, "./E2PilotRelativeAbundanceTotalsOut.csv")
PilotRel.df<-as.data.frame(PilotRelout)

Pilot.meta <- data.frame(sample_data(ps.E2pilotDe))
Pilot.meta$SampleID <- row.names(Pilot.meta)

# color pallet 
WTB_pallete = c("#053061", "#175169","#21777C", "#2A9D8F", "#8AB17D","#F5D279","#EFB366","#F4A261","#DE7350","#C8443F" )

# Relative Abundance plots - Experiment 1 supplemental Figure 2A 
ps.PilotREl <- transform(ps.E2pilotDe, "compositional")
ps.PilotREl<-taxglom(ps.PilotREl, taxrank = "Genus")
PilotRelout <- psmelt(ps.PilotREl)
PilotRel.df<-as.data.frame(PilotRelout)
OtherE2 <- PilotRel.df # the dataframe of the relative abundance within each sample
OtherE2 <- subset.data.frame(OtherE2, AntibioticTx =="None" )
OtherE2 <- subsel.data.frame(OtherE2, Abundance > 0.001)
OtherE2 <- OtherE2 %>% mutate(Genus = ifelse(Abundance < 0.05, "Other", Genus)) # replacing the genus catagory with "Other" if the abundance is less than 1% 
OtherE2 <- OtherE2 %>% mutate(Phylum = ifelse(Abundance < 0.05, "Other", Phylum)) 
OtherE2$DiseaseCINEndpoint <- factor(OtherE2$DiseaseCINEndpoint, levels = c("Hyperplasia", "CIN 2","CIN 3", "SCC", "No Tissue"))
table(OtherE2$Genus, OtherE2$Phylum)

OtherE2$Phylum <- factor(OtherE2$Phylum, levels = c("Other", "Actinobacteria",  "Firmicutes", "Bacteroidetes","Crenarchaeota","Fusobacteria", "Proteobacteria", "[Thermi]"))
OtherE2$Genus<- factor(OtherE2$Genus, levels = c("Other",
                                                 "Brachybacterium", "Corynebacterium","Microbacterium", "Micrococcus","Rothia","Saccharopolyspora",  "Streptomyces", ##Actinobacteria / actinomycetota (6)
                                                 
                                                 "Alloiococcus", "Bacillus","Enterococcus","Faecalibacterium",  "Lactobacillus","Roseburia",
                                                 "Staphylococcus", "Streptococcus",  "Turicibacter", # Firmicutes/ Bacillota 8)
                                                 
                                                 "Pedobacter", "Porphyromonas", "Sphingobacterium", # Bacteroidetes / Bacteroidota (3)
                                                 "Candidatus Nitrososphaera", # Crenarchaeota" (1)
                                                 
                                                 "Acidovorax", "Acinetobacter","Brevundimonas", "Enhydrobacter","Janthinobacterium","Neisseria", "Paracoccus", "Pseudomonas","Roseomonas",
                                                 "Stenotrophomonas")) #Proteobacteria /pseudomonadota (9)

BigMicro2 <-c("#C1C1C1", #Gray
              colorRampPalette(c("#031A34", "#053061", "#21777C", "#2A9D8F", "#8AB17D", "#C0C27B","#D2C77A", "#F5D279","#EFB366","#F4A261","#DE7350","#C8443F","#913C44", "#63292E"))(21))

ggplot(OtherE2, aes(x= Timepoint, y =Abundance, fill = Genus))+
  geom_bar(stat="identity", position="fill") +
  facet_wrap("Infected", scales = "free_x", nrow = 1)+
  ylab("Relative Abundance")+
  theme_light()+
  scale_fill_manual(values = BigMicro2)+
  scale_color_manual(values = BigMicro2)+
  theme(axis.text.x = element_text(size = 4, angle = 45)) +
  ggtitle("Genus present > 0.5% of reads in a sample")


# Alpha diversity (not shown in paper)- for non-antibitic treated mice in experiment 1 
AlphaTab <-microbiome::alpha(ps.E2pilotDe, index = c("observed", "chao1", "diversity_inverse_simpson", "diversity_gini_simpson","diversity_shannon",	"diversity_coverage", "evenness_camargo",	"evenness_pielou",
                                                    "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp", "dominance_dmn", "dominance_absolute", "dominance_relative", "dominance_simpson", 
                                                    "dominance_core_abundance", "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_rare_abundance"))
# write.csv(AlphaTab, "./Pilot_FINAL_AlphaTable.csv")

E2.NoABX <- subset_samples(ps.E2pilotDe, AntibioticTx == "None")
plot_richness(E2.NoABX, x = "TimeP", measures=c("Shannon"), color = "Cleared.PersistentAtEndpoint", shape = "DiseaseCINEndpoint ") + 
  geom_point(size = 6, width = 0.24)+
  #stat_summary(fun.y = mean,  geom = "point",shape = 95, size = 5, show.legend = F)+ 
  #facet_wrap(~ AntibioticTx, nrow = 1)+
  scale_color_manual(values = c("#053061", "#21777C",  "#8AB17D","#CCCCCC"))+
  theme_bw()+
  ggtitle("Shannon alpha Diversity for Experiment 1 Pilot Ex - NO antibiotic tx group")


# Beta Diversity (not shown in paper)
E2.NoABX <- subset_samples(ps.E2pilotDe, AntibioticTx == "None")
E2.not <-subset_samples(ps.E2pilot, AntibioticTx == "None")
E2.NoABX <- subset_samples(E2.NoABX, LK16Sid != "LKMB005_014")
min(sample_sums(E2.NoABX))# minimum sample read is 1343
median(sample_sums(E2.NoABX)) # 11156
median(sample_sums(E2.not))
max(sample_sums(E2.NoABX )) #46861

rarecurve(t(otu_table(E2.NoABX)), step=100, ylim =c(0,150), xlim=c(0,2000))  ## considering using a read cut off of 5000 for beta diversity metrics 
ps.rarefiedE2 = rarefy_even_depth(E2.NoABX, rngseed=2, sample.size=2000, replace=F) # THIS IS THE ONE TO GO WITH

GP = ps.rarefiedE2
GP.ord <- ordinate(GP, "NMDS",  "bray") # For Bray Curtis of the rareified dataset
plot_ordination(GP, GP.ord, type="samples", color="InfectStage2")+# shape = "EdgeCenter") + 
  geom_point(size=3) + 
  #geom_text(aes(label = LK16Sid), size = 3, hjust = 1.1)+
  ggtitle("Bray Curtis")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point( size = 6)+
  #scale_shape_manual(values = c(8,19,17,25,15,23))+
  scale_color_manual(values = c("#053061","#BBBBBB",  "#8AB17D","#EFB366","#DE7350","#C8443F" ))+
  #scale_fill_manual(values = c("#053061", "#21777C",  "#8AB17D","#EFB366","#DE7350","#C8443F" ))+
  theme_bw()
plot_ordination(GP, GP.ord, type="samples", color="Timepoint")+# shape = "EdgeCenter") + 
  geom_point(size=3) + ggtitle("Bray Curtis for the mice w/o antibiotics")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point(aes(shape = DiseaseCINEndpoint), size = 6)+
  scale_shape_manual(values = c(19,17,23,8, 15))+
  scale_color_manual(values = c("#053061", "#21777C","#2A9D8F",  "#8AB17D","#EFB366","#DE7350","#C8443F" ))+
  theme_bw()


#PERMANOVAS for Table 2 - Experiment 1  NO ANTIBIPTIOC treatment mice only 
sampledf.E2 <- data.frame(sample_data(ps.rarefiedE2)) 
E2_bray <- phyloseq::distance(ps.rarefiedE2, method = "bray")

adonis2(E2_bray ~ InfectStage2 , by= "margin", data = sampledf.E2, permutations = 9999) # 
adonis2(E2_bray ~ InfectStage , by= "margin", data = sampledf.E2, permutations = 9999) # 
adonis2(E2_bray ~ InfectStageDEM , by= "margin", data = sampledf.E2, permutations = 9999) # 
adonis2(E2_bray ~ DiseaseCINEndpoint  , by= "margin", data = sampledf.E2, permutations = 9999) # 
adonis2(E2_bray ~ DiseaseEndpointRanked , by= "margin", data = sampledf.E2, permutations = 9999) # 
adonis2(E2_bray ~ ViralCopyNumberHighLow  , by= "margin", data = sampledf.E2, permutations = 9999) # 
adonis2(E2_bray ~ Cleared.PersistentAtEndpoint  , by= "margin", data = sampledf.E2, permutations = 9999) # 

# Just the POST infection groups 
post.inf.rare <- subset_samples(ps.rarefiedE2, InfectStage2 == "Infected")
sampledf.E2 <- data.frame(sample_data(post.inf.rare)) 
E2_bray <- phyloseq::distance(post.inf.rare, method = "bray")

adonis2(E2_bray ~ InfectStage2 , by= "margin", data = sampledf.E2, permutations = 9999) # 
adonis2(E2_bray ~ InfectStage , by= "margin", data = sampledf.E2, permutations = 9999) # 
adonis2(E2_bray ~ InfectStageDEM , by= "margin", data = sampledf.E2, permutations = 9999) # 
adonis2(E2_bray ~ DiseaseCINEndpoint  , by= "margin", data = sampledf.E2, permutations = 9999) # 
adonis2(E2_bray ~ DiseaseEndpointRanked , by= "margin", data = sampledf.E2, permutations = 9999) # 
adonis2(E2_bray ~ ViralCopyNumberHighLow  , by= "margin", data = sampledf.E2, permutations = 9999) # 
adonis2(E2_bray ~ Cleared.PersistentAtEndpoint  , by= "margin", data = sampledf.E2, permutations = 9999) # 


## MIXOMICS -GENUS level - Experiment 1 - Figure 1 and Supplemental figure 2 
## the microbiome data - they highly recommend normalizing - 
#     i chose to import the relative abundance table (everything normalized to out of 1) 
#      of the whole data set so all samples were represented) then filtering out anything with an abundance
#      < 0.01 (1%) and calling it "other"-
## Note - all the dataframes need the same number of rows. (ie 10 rows - one per subjecr in subject order)
E2.NoABX <- subset_samples(ps.E2pilotDe, AntibioticTx == "None")
E2.NoABX <- tax_glom(E2.NoABX , "Genus")
ps.E2clean <- phyloseq_filter_prevalence(E2.NoABX, prev.trh = 0.05, abund.trh = 10, threshold_condition = "AND") # using AND here b/c need to ensure that the ASV is at least in 2 (>0.04 % of samples)
relE2.2 <- transform(ps.E2clean, "compositional") 
relE2.2.df <-psmelt(relE2.2)
relE2.2.df<- subset.data.frame(relE2.2.df, Abundance > 0.001)
# write.csv(relE2.2.df, "./Mixomics_Genus/E2relativeAbundance_Genus_Mixomics.csv")
E2.relativeAbundanceOTHER_pivot <- read_csv("./Mixomics_Genus/E2relativeAbundance_Genus_Mixomics.Pivot.csv")

Micro.df <- as.data.frame(E2.relativeAbundanceOTHER_pivot) # 77 ASVs included 
rownames(Micro.df) <- Micro.df$SampleID
Micro.df <- Micro.df[,-1]

E2.meta.df <- data.frame(sample_data(ps.E2clean))

#Natural microbiome - pre v. post depo - Experiment 1 Figure 1F 
E2.Nat.Meta<- subset.data.frame(E2.meta.df , NaturalMicro %in% c("Natural", "Post-Depo"))
Micro.Nat.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E2.Nat.Meta))
X.Nat <- Micro.Nat.df
Y <- E2.Nat.Meta
Y.Nat <- Y$NaturalMicro
Y.Nat<-factor(Y.Nat, levels = c("Natural", "Post-Depo"))
plsda.Nat <- plsda(X.Nat, Y.Nat, ncomp = 3)
plotIndiv(plsda.Nat, 
          ind.names = F, 
          ellipse = TRUE,
          comp = c(1,2),
          legend=TRUE, cex=6,
          col =c("#BBBBBB","#F4A261"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = "Eexperiment 1 NATURAL Pre-infect, Depo")

cord1_Load <-plotLoadings(plsda.Nat, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Natural == "TRUE"& Contrib.Post.Depo == "TRUE") # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Natural", "Post.Depo",  "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.Nat, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Natural == "TRUE"& Contrib.Post.Depo =="TRUE"  ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Natural", "Post.Depo",  "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))
Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.1)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.4)+
  scale_color_manual(values = c("#BBBBBB", "#F4A261"))+
  scale_fill_manual(values = c( "#BBBBBB", "#F4A261"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Natural", "Post-Depo", "Untouched","Mock", "Infected")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 1 Natural v Depo top driving Taxa')


# PLSDA for experiment 1 (no antibitoic treatment) pre infected v infected Supplemental figure 2
E2.Meta<- subset.data.frame(E2.meta.df , Infected %in% c("Infected"))
Micro.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E2.Meta))

X <- Micro.df
Y <- E2.Meta
Y.I <- Y$InfectStage2
Y.I<-factor(Y.I, levels = c( "Pre-Infection","Infected"))

plsda.I <- plsda(X, Y.I, ncomp = 3)
plotIndiv(plsda.I, 
          ind.names = F, 
          ellipse = TRUE,
          comp = c(1,2),
          legend=TRUE, cex=6,
          col =c("#BBBBBB","#053061"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 1  pre-Infection v. Infection')

cord1_Load<-plotLoadings(plsda.I, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Infected== "TRUE" & Contrib.Pre.Infection == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Pre.Infection", "Infected", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.I, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Pre.Infection== "TRUE" & Contrib.Infected == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Pre.Infection", "Infected", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))
Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.1)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.4)+
  scale_color_manual(values = c("#053061","#BBBBBB"))+
  scale_fill_manual(values = c("#053061","#BBBBBB"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Pre-Infection", "Infected")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 1- pre-infection v. Invection.top driving Taxa')

# infection stage (pre, depo, est, early, mid) Experiment 1 supplementalfigure 2 
E2.Ist.Meta<- subset.data.frame(E2.Meta, InfectStageDEM %in% c("Pre-Infection", "Post-Depo", "Infection Establishment (2wpi)", "Early Infection (4 & 6wpi)", "Mid Infection (8 & 10 wpi)"))
Micro.Ist.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E2.Ist.Meta))

X.Ist <-Micro.Ist.df
Y.Ist <- E2.Ist.Meta$InfectStageDEM
Y.Ist<-factor(Y.Ist, levels = c("Pre-Infection", "Post-Depo", "Infection Establishment (2wpi)", "Early Infection (4 & 6wpi)", "Mid Infection (8 & 10 wpi)"))

plsda.Ist <- plsda(X.Ist, Y.Ist, ncomp = 3) # note no SCC in this group
plotIndiv(plsda.Ist, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp= c(1,3),
          #pch = c(3,4,5,2,1),
          col =c("#BBBBBB","#EFB366","#268A86","#13546F","#053061" ), #"#DE7350"), 
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 1 infection stage')

cord1_Load <-plotLoadings(plsda.Ist, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Pre.Infection == "TRUE" & Contrib.Post.Depo == "TRUE" & Contrib.Infection.Establishment..2wpi.=="TRUE"& Contrib.Early.Infection..4...6wpi. == "TRUE" & Contrib.Mid.Infection..8...10.wpi.=="TRUE") # removing variables that contributed to ALL the groups
# "Normal/Hyperplasia", "CIN 3", "At least CIN3", "SCC"
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Contrib.Pre.Infection","Contrib.Post.Depo", "Contrib.Infection.Establishment..2wpi.", "Contrib.Early.Infection..4...6wpi.","Contrib.Mid.Infection..8...10.wpi.",  "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.Ist, comp = 3, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Pre.Infection == "TRUE" & Contrib.Post.Depo == "TRUE" & Contrib.Infection.Establishment..2wpi.=="TRUE"& Contrib.Early.Infection..4...6wpi. == "TRUE" & Contrib.Mid.Infection..8...10.wpi.=="TRUE") # removing variables that contributed to ALL the groups
# removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Contrib.Pre.Infection","Contrib.Post.Depo", "Contrib.Infection.Establishment..2wpi.", "Contrib.Early.Infection..4...6wpi.","Contrib.Mid.Infection..8...10.wpi.",  "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.05)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(0,0.45)+
  scale_color_manual(values = c("#13546F","#268A86","#053061", "#EFB366","#BBBBBB","#DE7350"))+
  scale_fill_manual(values = c("#13546F","#268A86","#053061", "#EFB366","#BBBBBB","#DE7350"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Pre-Infection", "Post-Depo", "Infection Establishment (2wpi)", "Early Infection (4 & 6wpi)", "Mid Infection (8 & 10 wpi)")), scales = "free_y", nrow =1)+
  theme_bw()+
  ggtitle('Experiment 1 Infection Stage- top driving Taxa')

# Infected mice (post infection) Infection stage
E2.I.Meta<- subset.data.frame(E2.Meta , InfectStage2 %in% c("Infected"))
Micro.I.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E2.I.Meta))

X.Ist <-Micro.I.df
Y.Ist <- E2.I.Meta$InfectStageDEM
Y.Ist<-factor(Y.Ist, levels = c("Infection Establishment (2wpi)", "Early Infection (4 & 6wpi)", "Mid Infection (8 & 10 wpi)"))

plsda.Ist <- plsda(X.Ist, Y.Ist, ncomp = 3) # note no SCC in this group
plotIndiv(plsda.Ist, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp= c(1,2),
          #pch = c(3,4,5,2,1),
          col =c("#268A86","#13546F","#053061" ), #"#DE7350"), 
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 1 infection stage')

cord1_Load <-plotLoadings(plsda.Ist, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings,  Contrib.Infection.Establishment..2wpi.=="TRUE"& Contrib.Early.Infection..4...6wpi. == "TRUE" & Contrib.Mid.Infection..8...10.wpi.=="TRUE") # removing variables that contributed to ALL the groups
# "Normal/Hyperplasia", "CIN 3", "At least CIN3", "SCC"
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable", "Contrib.Infection.Establishment..2wpi.", "Contrib.Early.Infection..4...6wpi.","Contrib.Mid.Infection..8...10.wpi.",  "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.Ist, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings,  Contrib.Infection.Establishment..2wpi.=="TRUE"& Contrib.Early.Infection..4...6wpi. == "TRUE" & Contrib.Mid.Infection..8...10.wpi.=="TRUE") # removing variables that contributed to ALL the groups
# removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable", "Contrib.Infection.Establishment..2wpi.", "Contrib.Early.Infection..4...6wpi.","Contrib.Mid.Infection..8...10.wpi.",  "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.15)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(0,0.45)+
  scale_color_manual(values = c("#13546F","#268A86","#053061" ))+
  scale_fill_manual(values = c("#13546F","#268A86","#053061"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Infection Establishment (2wpi)", "Early Infection (4 & 6wpi)", "Mid Infection (8 & 10 wpi)")), scales = "free_y", nrow =1)+
  theme_bw()+
  ggtitle('Experiment 1 Infection Stage- top driving Taxa')

# INFECTED mice only viral load H L 
E2.I.Meta<- subset.data.frame(E2.Meta , InfectStage2 %in% c("Infected"))
Micro.I.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E2.I.Meta))

E2.I.VHL.Meta<- subset.data.frame(E2.I.Meta , ViralCopyNumberHighLow %in% c("Low", "High"))
Micro.I.VHL.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E2.I.VHL.Meta))

X.VHL <-Micro.I.VHL.df
Y.VHL <- E2.I.VHL.Meta$ViralCopyNumberHighLow
Y.VHL <- factor(Y.VHL, levels = c("Low", "High" ))

E2.I.PC.Meta<- subset.data.frame(E2.I.Meta, Cleared.PersistentAtEndpoint %in% c("Cleared", "Persistent"))
Micro.I.PC.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E2.I.PC.Meta))

plsda.VHL <- plsda(X.VHL, Y.VHL, ncomp = 3)
plotIndiv(plsda.VHL, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp = c(1,2),
          col =c(  "#2A9D8F", "#053061"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 1 infected mice only Viral Load')

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

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.15)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.45)+
  scale_color_manual(values = c("#053061","#2A9D8F"))+
  scale_fill_manual(values = c( "#053061","#2A9D8F"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Low", "High")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 1 infected mice - Viral Load, top driving Taxa')

# INFECTED mice only DiseaseSeverity
E2.I.Meta<- subset.data.frame(E2.Meta , InfectStage2 %in% c("Infected"))

E2.I.CIN.Meta<- subset.data.frame(E2.I.Meta , DiseaseCINEndpoint %in% c("CIN 2", "CIN 3", "SCC"))
Micro.I.CIN.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E2.I.CIN.Meta))

X.CIN <-Micro.I.CIN.df
Y.CIN <- E2.I.CIN.Meta$DiseaseCINEndpoint
Y.CIN <- factor(Y.CIN, levels = c("CIN 2", "CIN 3", "SCC" ))

plsda.CIN <- plsda(X.CIN, Y.CIN, ncomp = 3)
plotIndiv(plsda.CIN, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp = c(1,2),
          col =c("#8AB17D","#F6D27A","#F5A261"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'Experiment 1 infected mice only disease severity')

cord1_Load<-plotLoadings(plsda.CIN, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.CIN.2== "TRUE" & Contrib.CIN.3 == "TRUE"  & Contrib.SCC == "TRUE") # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","CIN.2", "CIN.3", "SCC", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.CIN, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.CIN.2== "TRUE" & Contrib.CIN.3 == "TRUE"  & Contrib.SCC == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","CIN.2", "CIN.3", "SCC", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.115)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.3)+
  scale_color_manual(values = c("#8AB17D","#F6D27A","#F5A261"))+
  scale_fill_manual(values = c( "#8AB17D","#F6D27A","#F5A261"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("CIN 2", "CIN 3", "SCC")), scales = "free_y")+
  theme_bw()+
  ggtitle('Experiment 1 infected mice -Disease Severity, top driving Taxa')

# MAASLIN For Experiment 1 - Supplemental figures 1 and 2
E2pilotRElgenus <- read_csv("./E2pilotRElgenus.csv")
input_genus_data <- as.data.frame(E2pilotRElgenus)
rownames(input_genus_data) <-input_genus_data[,1]

input_meta_file <-Pilot.meta
input_meta_file$DiseaseCINEndpoint <- factor(input_meta_file$DiseaseCINEndpoint, levels = c("Pre-Infection", "CIN 2","CIN 3", "SCC", "No Tissue"))

NoABX_input_meta_file <-subset.data.frame(input_meta_file, AntibioticTx == "None")
NoABX_input_meta_file <- subset.data.frame(NoABX_input_meta_file, DiseaseCINEndpoint != "No Tissue")
VC.postInfection <- subset.data.frame(NoABX_input_meta_file, ViralCopyNumberHighLow != "Pre-Infection")
Post_InfectE2.MF <- subset.data.frame(NoABX_input_meta_file, DiseaseCINEndpoint !="Pre-Infection" )

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = NoABX_DISEASE_input_meta_file, 
  output = "./MAASLIN/E2_Post-InfectDiseaseSCC_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "SCC"))

# natural microbiome v. post depo Supplemental Figure 1 
E2.natural.meta <- subset.data.frame(input_meta_file, NaturalMicro %in% c("Natural", "Post-Depo"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = E2.natural.meta, 
  output = "./MAASLIN/E2_Natural.v.PostDepo_genus.", 
  fixed_effects = c("NaturalMicro"),
  reference = c("NaturalMicro", "Natural"),
  random_effects = c("MP"))
                                    

# Evaluating the effects of the virus in the group NOT treated with antibiotics - Supplemental Figure 2 
NoABX.Pre.Post_input_meta_file <-subset.data.frame(input_meta_file, AntibioticTx == "None")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = NoABX.Pre.Post_input_meta_file ,
  output = "./MAASLIN/E2_Pre.v.PostInfection_genus_maslin", 
  fixed_effects = c("InfectStage2"),
  reference = c("InfectStage2", "Pre-Infection"),
  random_effects = c("MP"))

#ViralCopy number - POST infection 
VC.postInfection <- subset.data.frame(NoABX_input_meta_file, ViralCopyNumberHighLow != "Pre-Infection")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = VC.postInfection, 
  output = "./MAASLIN/E2_VIralCopyNumber_post_infection_genus_maslin", 
  fixed_effects = c("ViralCopyNumberHighLow"),
  reference = c("ViralCopyNumberHighLow", "Low"),
  random_effects = c("MP"))

# DIsease severity post infection 
postInfection <- subset.data.frame(NoABX_input_meta_file, InfectStage2!= "Pre-Infection")
postInfection.CIN <- subset.data.frame(NoABX_input_meta_file, DiseaseCINEndpoint %in% c("CIN 2", "CIN 3", "SCC"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = postInfection.CIN, 
  output = "./MAASLIN/E2_DiseaseCIN_postInfection_CIN2_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "CIN 2"),
  random_effects = c("MP"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = postInfection.CIN, 
  output = "./MAASLIN/E2_DiseaseCIN_postInfection_SCC_genus_maslin", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "SCC"),
  random_effects = c("MP"))

# Infection stage 
postInfectionStage <- subset.data.frame(NoABX_input_meta_file, InfectStage2 != "Pre-Infection")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = postInfectionStage, 
  output = "./MAASLIN/E2_InfectStage_Establishment_genus_maslin", 
  fixed_effects = c("InfectStageDEM"),
  reference = c("InfectStageDEM", "Infection Establishment (2wpi)"),
  random_effects = c("MP"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = postInfectionStage, 
  output = "./MAASLIN/E2_InfectStage_Early_genus_maslin", 
  fixed_effects = c("InfectStageDEM"),
  reference = c("InfectStageDEM", "Early Infection (4 & 6wpi)"),
  random_effects = c("MP"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = postInfectionStage, 
  output = "./MAASLIN/E2_InfectStage_MidInf_genus_maslin", 
  fixed_effects = c("InfectStageDEM"),
  reference = c("InfectStageDEM", "Mid Infection (8 & 10 wpi)"),
  random_effects = c("MP"))

