# Cleaned R code for the Baseline Experiment 

library(tidyverse)
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
library(mixOmics)

setwd("/Users/liztown/Documents/KalanLab/Papers/Vaginal Microbiome/Rcode_forManuscript/Experiment2_Baseline")

ps.E3Base <- qza_to_phyloseq(features = "./table-Baseline-dada2.qza",
                             taxonomy = "./taxonomy-Baseline.qza",
                             metadata = "./E3BaselineManufest.txt",
                             tree = "./Baseline-rooted-tree.qza")

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
ps.Decontam <- subset_taxa(ps.Decontam, Genus !="Candidatus Nitrososphaera") 
ps.E3BaseDe <- phyloseq_filter_prevalence(ps.Decontam, prev.trh = 0.05, abund.trh = 10, threshold_condition = "OR")
ps.E3BaseDe <- subset_samples(ps.E3BaseDe, SampleOrigin == "Mouse") # removing the negative control
ps.E3BaseDe# THE CLEANED FILE 

OutCV <- psmelt(ps.E3BaseDe)
#write.csv(OutCV, "/./E3BaseAbundanceTotalsOut.csv")

ps.baseREl <- transform(ps.E3BaseDe, "compositional")
BaseRelout <- psmelt(ps.baseREl)
#write.csv(BaseRelout, "./E3BaseRelativeAbundanceTotalsOut.csv")
BaseRel.df<-as.data.frame(BaseRelout)

Base.meta <- data.frame(sample_data(ps.E3BaseDe))
Base.meta$SampleID <- row.names(Base.meta)

# color pallet 
WTB_pallete = c("#053061", "#175169","#21777C", "#2A9D8F", "#8AB17D","#F5D279","#EFB366","#F4A261","#DE7350","#C8443F" )

# Relative abundance Figure 2
ps.baseREl <- transform(ps.E3BaseDe, "compositional")
ps.baseREl<- tax_glom(ps.baseREl, taxrank = "Genus")
BaseRelout <- psmelt(ps.baseREl)
OtherE3<-as.data.frame(BaseRelout)# the dataframe of the relative abundance within each sample
OtherE3 <- subset.data.frame(OtherE3, AntibioticTx =="None")
OtherE3 <- subset.data.frame(OtherE3, Abundance > 0.001)
OtherE3 <- OtherE3 %>% mutate(Genus = ifelse(Abundance < 0.05, "Other", Genus)) # replacing the genus catagory with "Other" if the abundance is less than 1% 
OtherE3 <- OtherE3 %>% mutate(Phylum = ifelse(Abundance < 0.05, "Other", Phylum)) 
OtherE3 <- subset.data.frame(OtherE3, Genus != "Candidatus Nitrososphaera")
OtherE3[OtherE3 == "[Prevotella]"] <- "Prevotella"
OtherE3[OtherE3 == "[Ruminococcus]"] <- "Ruminococcus"
OtherE3[OtherE3 == "Propionibacterium"] <- "Cutibacterium"
OtherE3$Infected <- factor(OtherE3$Infected, levels = c("Untouched", "Mock","Infected"))
table(OtherE3$Genus, OtherE3$Phylum)

OtherE3$Genus<- factor(OtherE3$Genus, levels = c("Other",
                                                 "Actinomyces","Actinomycetospora","Actinotalea", "Arthrobacter","Bifidobacterium","Blastococcus","Brevibacterium","Collinsella",
                                                 "Corynebacterium","Cryocola","Cutibacterium", "Dermacoccus","Friedmanniella", "Kocuria", "Leucobacter", "Micrococcus","Mycobacterium","N09","Nocardiopsis",
                                                 "Propionibacterium","Rhodococcus","Rothia","Saccharopolyspora", "Sanguibacter","Serinicoccus", "Streptomyces", ##Actinobacteria / actinomycetota (24)
                                                 
                                                 "Aerococcus", "Alloiococcus", "Anaerobacillus", "Anaerococcus", "Anoxybacillus","Bacillus","Blautia", "Brevibacillus","Catenibacterium","Clostridium",
                                                 "Coprococcus","Desemzia","Enterococcus","Faecalibacterium", "Gemella","Gemmiger","Geobacillus", "Lactobacillus","Lactococcus","Megasphaera",
                                                 "Oscillospira","Paenibacillus","Roseburia","Ruminococcus","Sporosarcina",
                                                 "Staphylococcus", "Streptococcus", "Thermicanus", "Turicibacter","Veillonella", # Firmicutes/ Bacillota (30)
                                                 
                                                 "Bacteroides", "Cloacibacterium","Cytophaga", "Dyadobacter","Hymenobacter","Parabacteroides","Pedobacter", "Porphyromonas",
                                                 "Prevotella","[Prevotella]", "Sphingobacterium", # Bacteroidetes / Bacteroidota (11)
                                                 "Mucispirillum", # Deferribacteres/ Deferribacterota
                                                 "Fusobacterium", #Fusobacteria/ Fusobacteriota
                                                 "Gemmata",  #
                                                 
                                                 "Acidovorax", "Acinetobacter","Agrobacterium","Brevundimonas", "Burkholderia","Cellvibrio","Haemophilus","Hydrogenophilus","Lautropia",
                                                 "Massilia","Neisseria","Novosphingobium", "Paracoccus","Providencia", "Pseudomonas","Pseudoxanthomonas","Rhodanobacter", "Rhodoplanes",
                                                 "Rubellimicrobium","Sphingobium", "Stenotrophomonas","Steroidobacter","Xylella", #Proteobacteria /pseudomonadota (24)
                                                 
                                                 "Deinococcus" )) # Thermi

BigMicro2 <-c("#C1C1C1", #Gray
              colorRampPalette(c("#031A34", "#053061", "#21777C", "#2A9D8F", "#8AB17D", "#C0C27B","#D2C77A", "#F5D279","#EFB366","#F4A261","#DE7350","#C8443F","#913C44", "#733036"))(58))

ggplot(OtherE3, aes(x= Timepoint, y =Abundance, fill = Genus))+
  geom_bar(stat="identity", position="fill") +
  facet_wrap("Infected", scales = "free_x", nrow = 1)+
  ylab("Relative Abundance")+
  theme_light()+
  scale_fill_manual(values = BigMicro2)+
  scale_color_manual(values = BigMicro2)+
  theme(axis.text.x = element_text(size = 4, angle = 45)) +
  ggtitle("Experiment 2 Genus present > 0.5% of reads in a sample")

# Alpha abundance (not in paper)
plot_richness(ps.E3BaseDe, x = "TimeP", measures=c("Shannon"), color = "ABX.Infect") + 
  geom_jitter(size = 2.5, width = 0.24)+
  stat_summary(fun.y = mean,  geom = "point",shape = 95, size = 5, show.legend = F)+ 
  facet_wrap(~ ABX.Infect, nrow = 1)+
  scale_color_manual(values = c("#053061", "#21777C",  "#8AB17D","#EFB366","#DE7350","#C8443F" ))+
  theme_bw()+
  ggtitle("Shannon alpha Diversity for Experiment 2 Baseline Ex")

AlphaTab <-microbiome::alpha(ps.E3BaseDe, index = c("observed", "chao1", "diversity_inverse_simpson", "diversity_gini_simpson","diversity_shannon",	"diversity_coverage", "evenness_camargo",	"evenness_pielou",
                                                   "evenness_simpson", "evenness_evar", "evenness_bulla", "dominance_dbp", "dominance_dmn", "dominance_absolute", "dominance_relative", "dominance_simpson", 
                                                   "dominance_core_abundance", "dominance_gini", "rarity_log_modulo_skewness", "rarity_low_abundance", "rarity_rare_abundance"))
ps.NOABXE3BaseDe <- subset_samples(ps.E3BaseDe, AntibioticTx == "None") 
plot_richness(ps.NOABXE3BaseDe, x = "TimeP", measures=c("Shannon"), color = "Infected") + 
  geom_point(size = 5)+
  stat_summary(fun.y = mean,  geom = "point",shape = 95, size = 6, show.legend = F)+ 
  facet_wrap(~ Infected, nrow = 1)+
  scale_color_manual(values = c("#053061",   "#8AB17D","#EFB366" ))+
  theme_bw()+
  ggtitle("Shannon alpha Diversity for Experiment 2- Baseline Ex-- not antibiotic treated mice")

# BEta Diversity Figure 2B
min(sample_sums(ps.E3BaseDe))# minimum sample read is 0
median(sample_sums(ps.E3Base)) # 10589
median(sample_sums(ps.E3BaseDe)) # 1126
max(sample_sums(ps.E3BaseDe)) #33601
ps.NOABXE3BaseDe <- subset_samples(ps.E3BaseDe, AntibioticTx == "None") # removing mice treated with antibitotics
ps.Base.no <-subset_samples(ps.E3Base, AntibioticTx == "None")# removing mice treated with antibitotics
median(sample_sums(ps.Base.no))
min(sample_sums(ps.NOABXE3BaseDe))# minimum sample read is 11
median(sample_sums(ps.NOABXE3BaseDe)) # 1177
max(sample_sums(ps.NOABXE3BaseDe)) #32808
table(sample_sums(ps.NOABXE3BaseDe))

rarecurve(t(otu_table(ps.NOABXE3BaseDe)), step=100, ylim =c(0,300), xlim=c(0,1500))  ## considering using a read cut off of 5000 for beta diversity metrics 

Rare.NoABX = rarefy_even_depth(ps.NOABXE3BaseDe, rngseed=1, sample.size=300, replace=F) # THIS IS THE ONE TO GO WITH
Rare.NoABX@sam_data$InfectStagePreUM <- factor(Rare.NoABX@sam_data$InfectStagePreUM, levels = c("Pre-Infection", "Post-Depo","Untouched", "Mock", "Infected"))

GP = Rare.NoABX
GP.ord <- ordinate(GP, "NMDS",  "bray") # For Bray Curtis of the rareified dataset
plot_ordination(GP, GP.ord, type="samples", color="InfectStagePreUM")+# shape = "EdgeCenter") + 
  geom_point(size=3) + ggtitle("Bray Curtis for the mice w/o antibiotics- Iinfect stage ")+
  stat_ellipse(type = "t", level = 0.9)+
  geom_point( size = 6)+
  #scale_shape_manual(values = c(8,19,17,25,15,23))+
  scale_color_manual(values = c("#BBBBBB","#F4A261","#EFB366",  "#8AB17D","#053061"))+
  #WTB_pallete = c("#053061", "#175169","#21777C", "#2A9D8F", "#8AB17D","#F5D279","#EFB366","#F4A261","#DE7350","#C8443F" )
  #scale_fill_manual(values = c("#053061", "#21777C",  "#8AB17D","#EFB366","#DE7350","#C8443F" ))+
  theme_bw()

Rare.NoABX.Inf <- subset_samples(Rare.NoABX,Infected == "Infected")
Rare.NoABX.Inf <- subset_samples(Rare.NoABX.Inf, !LK16Sid %in% c("LK16S004_057","LK16S004_071","LK16S004_072"))
Rare.NoABX.Inf@sam_data$InfectStageDEM <- factor(Rare.NoABX.Inf@sam_data$InfectStageDEM, levels = c("Pre-Infection", "Post-Depo", "Infection Establishment (2wpi)", "Early Infection (4 & 6wpi)", "Mid Infection (8-10-12wpi)"))

GP = Rare.NoABX.Inf
GP.ord <- ordinate(GP, "NMDS",  "bray") # For Bray Curtis of the rareified dataset
plot_ordination(GP, GP.ord, type="samples", color="InfectStageDEM")+# shape = "EdgeCenter") + 
  geom_point(size=3) + ggtitle("Bray Curtis for the mice w/o antibiotics- Iinfect stage ")+
  stat_ellipse(type = "t", level = 0.9)+
  #geom_text(aes(label = LK16Sid))+
  geom_point( size = 6)+
  #scale_shape_manual(values = c(8,19,17,25,15,23))+
  scale_color_manual(values = c("#BBBBBB","#EFB366","#268A86","#13546F","#053061"))+
  #WTB_pallete = c("#053061", "#175169","#21777C", "#2A9D8F", "#8AB17D","#F5D279","#EFB366","#F4A261","#DE7350","#C8443F" )
  #scale_fill_manual(values = c("#053061", "#21777C",  "#8AB17D","#EFB366","#DE7350","#C8443F" ))+
  theme_bw()


# Permanova - for table 2
sampledf <- data.frame(sample_data(Rare.NoABX)) 
E3_bray <- phyloseq::distance(Rare.NoABX, method = "bray")

adonis2(E3_bray ~ Infected, by= "margin", data = sampledf, permutations = 9999) # Mock v. Untreated v. Infected 
adonis2(E3_bray ~ InfectStage2, by= "margin", data = sampledf, permutations = 9999) # Prev.post infect
adonis2(E3_bray ~ InfectStageDEM, by= "margin", data = sampledf, permutations = 9999) # pre, depo, post infection
adonis2(E3_bray ~ InfectStagePreUM, by= "margin", data = sampledf, permutations = 9999) # pre, depo, post infection
adonis2(E3_bray ~ ViralCopyNumberHighLow,  by= "margin", data = sampledf, permutations = 9999) # 
adonis2(E3_bray ~ DiseaseCINEndpoint, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(E3_bray ~ Cleared.Persistent, by= "margin", data = sampledf, permutations = 9999) # 

post <- subset_samples(Rare.NoABX, InfectStage2 =="Infected")
sampledf <- data.frame(sample_data(post)) 
E3_bray <- phyloseq::distance(post, method = "bray")

adonis2(E3_bray ~ Infected, by= "margin", data = sampledf, permutations = 9999) # Mock v. Untreated v. Infected 

Moc <- subset_samples(Rare.NoABX, Infected == "Mock")
sampledf <- data.frame(sample_data(Moc)) 
E3_bray <- phyloseq::distance(Moc, method = "bray")

adonis2(E3_bray ~ InfectStage2, by= "margin", data = sampledf, permutations = 9999) # Prev.post infect
adonis2(E3_bray ~ InfectStageDEM, by= "margin", data = sampledf, permutations = 9999) # pre, depo, post infection


Untou <- subset_samples(Rare.NoABX, Infected== "Untouched")
sampledf <- data.frame(sample_data(Untou)) 
E3_bray <- phyloseq::distance(Untou, method = "bray")

adonis2(E3_bray ~ InfectStage2, by= "margin", data = sampledf, permutations = 9999) # Prev.post infect
adonis2(E3_bray ~ InfectStageDEM, by= "margin", data = sampledf, permutations = 9999) # pre, depo, post infection


Rare.NoABX.Inf <- subset_samples(Rare.NoABX,Infected == "Infected")
Rare.NoABX.Inf <- subset_samples(Rare.NoABX.Inf, !LK16Sid %in% c("LK16S004_057","LK16S004_071","LK16S004_072"))
sampledf <- data.frame(sample_data(Rare.NoABX.Inf)) 
E3_bray <- phyloseq::distance(Rare.NoABX.Inf, method = "bray")

adonis2(E3_bray ~ InfectStage2, by= "margin", data = sampledf, permutations = 9999) # Prev.post infect
adonis2(E3_bray ~ InfectStageDEM, by= "margin", data = sampledf, permutations = 9999) # pre, depo, post infection
adonis2(E3_bray ~ ViralCopyNumberHighLow, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(E3_bray ~ DiseaseCINEndpoint, by= "margin", data = sampledf, permutations = 9999) # 
# adonis2(E3_bray ~ Cleared.Persistent, by= "margin", data = sampledf, permutations = 9999) # 


inf.Post <-subset_samples(Rare.NoABX.Inf, InfectStage2 == "Infected")
sampledf <- data.frame(sample_data(inf.Post)) 
E3_bray <- phyloseq::distance(inf.Post, method = "bray")

adonis2(E3_bray ~ InfectStageDEM, by= "margin", data = sampledf, permutations = 9999) # pre, depo, post infection
adonis2(E3_bray ~ ViralCopyNumberHighLow + Timepoint, by= "margin", data = sampledf, permutations = 9999) # 
adonis2(E3_bray ~ DiseaseCINEndpoint, by= "margin", data = sampledf, permutations = 9999) # 
# adonis2(E3_bray ~ Cleared.Percistent, by= "margin", data = sampledf, permutations = 9999) # 

 
## MIXOMICS - Genus level - for Figure 2 - Experiemtn 2
## the microbiome data - they highly recommend normalizing - 
#     i chose to import the relative abundance table (everything normalized to out of 1) 
#      of the whole data set so all samples were represented) then filtering out anything with an abundance
#      < 0.01 (1%) and calling it "other"-
## Note - all the dataframes need the same number of rows. (ie 10 rows - one per subjecr in subject order)
E3.NoABX <- subset_samples(ps.E3BaseDe, AntibioticTx == "None")
E3.NoABX <- tax_glom(E3.NoABX , "Genus")
ps.E3clean <- phyloseq_filter_prevalence(E3.NoABX, prev.trh = 0.03, abund.trh = 10, threshold_condition = "AND") # using AND here b/c need to ensure that the ASV is at least in 2 (>0.04 % of samples)
relE3.2 <- transform(ps.E3clean, "compositional") 
relE3.2.df <-psmelt(relE3.2)
relE3.2.df<- subset.data.frame(relE3.2.df, Abundance > 0.001)
#write.csv(relE3.2.df, "./Mixomics_Genus/E3relativeAbundance_Genus_Mixomics.csv")
E3.relativeAbundanceOTHER_pivot <- read_csv("./Mixomics_Genus/E3relativeAbundance_Genus_Mixomics.Pivot.csv")

Micro.df <- as.data.frame(E3.relativeAbundanceOTHER_pivot) # 269 ASVs included 
rownames(Micro.df) <- Micro.df$SampleID
Micro.df <- Micro.df[,-1]

E3.meta.df <- data.frame(sample_data(ps.E3clean))

# PLSDA for ALL infected vs. Mock vs. Untouched
E3.Meta<- subset.data.frame(E3.meta.df , Infected %in% c("Untouched", "Mock", "Infected"))
Micro.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E3.Meta))

X <- Micro.df
Y <- E3.Meta
Y.IUM <- Y$Infected
Y.IUM<-factor(Y.IUM, levels = c( "Untouched","Mock", "Infected"))

plsda.IUM <- plsda(X, Y.IUM, ncomp = 3)
plotIndiv(plsda.IUM, 
          ind.names = F, 
          ellipse = TRUE,
          comp = c(1,2),
          legend=TRUE, cex=6,
          col =c( "#EFB366",  "#8AB17D", "#053061"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = "E3 Untouched v. Mock v. Infection - all timepoints")

cord1_Load<-plotLoadings(plsda.IUM, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Infected== "TRUE" & Contrib.Mock == "TRUE"& Contrib.Untouched== "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Untouched","Mock", "Infected", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.IUM, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Mock== "TRUE" & Contrib.Infected == "TRUE" & Contrib.Untouched == "TRUE"  ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Untouched","Mock", "Infected", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))
Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.001)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.3)+
  scale_color_manual(values = c( "#053061", "#8AB17D", "#EFB366"))+
  scale_fill_manual(values = c("#053061", "#8AB17D", "#EFB366"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Untouched","Mock", "Infected")), scales = "free_y")+
  theme_bw()+
  ggtitle('E3 Untouched v. Mock v. Invection.top driving Taxa')


X <- Micro.df
Y <- E3.Meta
Y.PIUM <- Y$InfectStagePreUM
Y.PIUM<-factor(Y.PIUM, levels = c("Pre-Infection", "Post-Depo", "Untouched","Mock", "Infected"))

plsda.PIUM <- plsda(X, Y.PIUM, ncomp = 3)
plotIndiv(plsda.PIUM, 
          ind.names = F, 
          ellipse = TRUE,
          comp = c(1,2),
          legend=TRUE, cex=6,
          col =c("#BBBBBB","#F4A261","#EFB366", "#8AB17D","#053061"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = "E3 Pre-infect, Depo, Untouched v. Mock v. Infection - all timepoints")

cord1_Load <-plotLoadings(plsda.PIUM, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Pre.Infection == "TRUE"& Contrib.Post.Depo == "TRUE"& Contrib.Infected== "TRUE" & Contrib.Mock == "TRUE"& Contrib.Untouched== "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Pre.Infection", "Post.Depo", "Untouched","Mock", "Infected", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.PIUM, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Pre.Infection == "TRUE"& Contrib.Post.Depo =="TRUE" & Contrib.Mock== "TRUE" & Contrib.Infected == "TRUE" & Contrib.Untouched == "TRUE"  ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Pre.Infection", "Post.Depo", "Untouched","Mock", "Infected", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))
Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.05)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.35)+
  scale_color_manual(values = c("#053061", "#F4A261", "#BBBBBB","#EFB366", "#8AB17D","#053061"))+
  scale_fill_manual(values = c("#053061", "#F4A261", "#BBBBBB","#EFB366", "#8AB17D","#053061"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Pre-Infection", "Post-Depo", "Untouched","Mock", "Infected")), scales = "free_y")+
  theme_bw()+
  ggtitle('E3 Pre-infection, Depo, Untouched v. Mock v. Invection.top driving Taxa')

#Natural microbiome - pre v. post depo 
E3.Nat.Meta<- subset.data.frame(E3.meta.df , NaturalMicro %in% c("Natural", "Post-Depo"))
Micro.Nat.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E3.Nat.Meta))

X.Nat <- Micro.Nat.df
Y <- E3.Nat.Meta
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
          title = "E3 NATURAL Pre-infect, Depo")

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
  ggtitle('E3 Natural v Depo top driving Taxa')


# post infevtion timepoints 
E3.I.Meta<- subset.data.frame(E3.Meta, InfectStage2%in% c("Infected"))
Micro.I.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E3.I.Meta))

X <- Micro.I.df
Y <- E3.I.Meta
Y.IUM <- Y$Infected
Y.IUM<-factor(Y.IUM, levels = c( "Untouched","Mock", "Infected"))

plsda.IUM <- plsda(X, Y.IUM, ncomp = 3)
plotIndiv(plsda.IUM, 
          ind.names = F, 
          ellipse = TRUE,
          comp = c(1,2),
          legend=TRUE, cex=6,
          col =c( "#EFB366",  "#8AB17D", "#053061"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = "E3 Untouched v. Mock v. Infection - Post infection timepoints")

cord1_Load<-plotLoadings(plsda.IUM, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Infected== "TRUE" & Contrib.Mock == "TRUE"& Contrib.Untouched== "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Untouched","Mock", "Infected", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.IUM, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Mock== "TRUE" & Contrib.Infected == "TRUE" & Contrib.Untouched == "TRUE"  ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Untouched","Mock", "Infected", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))
Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.001)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.3)+
  scale_color_manual(values = c("#053061",  "#8AB17D", "#EFB366"))+
  scale_fill_manual(values = c("#053061",  "#8AB17D", "#EFB366"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Untouched","Mock", "Infected")), scales = "free_y")+
  theme_bw()+
  ggtitle('E3 Untouched v. Mock v. Invection.- POST infection timepoints top driving Taxa')

# Infected group ALL timepoints (infection stage: pre, depo, est, early, mid)
E3.Ist.Meta <- subset.data.frame(E3.Meta, Infected %in% c("Infected"))
Micro.Ist.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E3.Ist.Meta))

X.Ist <-Micro.Ist.df
Y.Ist <- E3.Ist.Meta$InfectStageDEM
Y.Ist<-factor(Y.Ist, levels = c("Pre-Infection", "Post-Depo", "Infection Establishment (2wpi)", "Early Infection (4 & 6wpi)", "Mid Infection (8-10-12wpi)"))

plsda.Ist <- plsda(X.Ist, Y.Ist, ncomp = 3) # note no SCC in this group
plotIndiv(plsda.Ist, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp= c(1,2),
          #pch = c(3,4,5,2,1),
          col =c("#BBBBBB","#EFB366","#268A86","#13546F","#053061" ), #"#DE7350"), 
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'E3 infection stage - infected mice only')

cord1_Load <-plotLoadings(plsda.Ist, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.Pre.Infection == "TRUE" & Contrib.Post.Depo == "TRUE" & Contrib.Infection.Establishment..2wpi.=="TRUE"& Contrib.Early.Infection..4...6wpi. == "TRUE" & Contrib.Mid.Infection..8.10.12wpi.=="TRUE") # removing variables that contributed to ALL the groups
# "Normal/Hyperplasia", "CIN 3", "At least CIN3", "SCC"
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Contrib.Pre.Infection","Contrib.Post.Depo", "Contrib.Infection.Establishment..2wpi.", "Contrib.Early.Infection..4...6wpi.","Contrib.Mid.Infection..8.10.12wpi.",  "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.Ist, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Pre.Infection == "TRUE" & Contrib.Post.Depo == "TRUE" & Contrib.Infection.Establishment..2wpi.=="TRUE"& Contrib.Early.Infection..4...6wpi. == "TRUE" & Contrib.Mid.Infection..8.10.12wpi.=="TRUE") # removing variables that contributed to ALL the groups
# removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Contrib.Pre.Infection","Contrib.Post.Depo", "Contrib.Infection.Establishment..2wpi.", "Contrib.Early.Infection..4...6wpi.","Contrib.Mid.Infection..8.10.12wpi.",  "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.05)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(0,0.3)+
  scale_color_manual(values = c("#13546F","#268A86","#053061", "#EFB366","#BBBBBB","#DE7350"))+
  scale_fill_manual(values = c("#13546F","#268A86","#053061", "#EFB366","#BBBBBB","#DE7350"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Pre-Infection", "Post-Depo", "Infection Establishment (2wpi)", "Early Infection (4 & 6wpi)", "Mid Infection (8-10-12wpi)")), 
             scales = "free_y", nrow =1)+
  theme_bw()+
  ggtitle('E3 Infection Stage- top driving Taxa')

X.I <-Micro.Ist.df
Y.I <- E3.Ist.Meta$InfectStage2
Y.I<-factor(Y.I, levels = c( "Pre-Infection","Infected"))

plsda.I <- plsda(X.I, Y.I, ncomp = 3)
plotIndiv(plsda.I, 
          ind.names = F, 
          ellipse = TRUE,
          comp = c(1,2),
          legend=TRUE, cex=6,
          col =c("#BBBBBB","#053061"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'E3 pre-Infection v. Infection')

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
  ggtitle('E3- pre-infection v. Invection.top driving Taxa')


# infected group at post infection timepoints 
E3.Ist.Meta <- subset.data.frame(E3.Meta, Infected %in% c("Infected"))
E3.PostI<- subset.data.frame(E3.Ist.Meta, InfectStage2 == "Infected")
#Micro.PostI.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E3.PostI.Meta))

E3.I.VHL.Meta<- subset.data.frame(E3.PostI, ViralCopyNumberHighLow %in% c("Low", "High"))
Micro.I.VHL.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E3.I.VHL.Meta))

X.VHL <-Micro.I.VHL.df
Y.VHL <- E3.I.VHL.Meta$ViralCopyNumberHighLow
Y.VHL <- factor(Y.VHL, levels = c("Low", "High" ))

E3.I.dis.Meta <- subset.data.frame(E3.PostI, DiseaseCINEndpoint %in% c("CIN2", "CIN3", "CIN2 (considered 3)"))
E3.I.dis.Meta[E3.I.dis.Meta == "CIN2 (considered 3)"] <- "CIN3"
Micro.I.dis.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E3.I.dis.Meta))

X.dis <-Micro.I.dis.df
Y.dis <- E3.I.dis.Meta$DiseaseCINEndpoint
Y.dis <- factor(Y.dis, levels = c("CIN2", "CIN3" ))

plsda.VHL <- plsda(X.VHL, Y.VHL, ncomp = 3)
plotIndiv(plsda.VHL, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp = c(1,2),
          col =c(  "#2A9D8F", "#053061"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'E3 infected mice only Viral Load')

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

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.2)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.35)+
  scale_color_manual(values = c("#053061","#2A9D8F"))+
  scale_fill_manual(values = c( "#053061","#2A9D8F"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Low", "High")), scales = "free_y")+
  theme_bw()+
  ggtitle('E3 infected mice post infection - Viral Load, top driving Taxa')


plsda.dis <- plsda(X.dis, Y.dis, ncomp = 3)
plotIndiv(plsda.dis, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp = c(1,2),
          col =c("#8AB17D","#F6D27A"),
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'E3 infected mice disease severity')

cord1_Load<-plotLoadings(plsda.dis, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings, Contrib.CIN2== "TRUE" & Contrib.CIN3 == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","CIN2", "CIN3", "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.dis, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.CIN2== "TRUE" & Contrib.CIN3 == "TRUE" ) # removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","CIN2", "CIN3", "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.05)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(-0,0.35)+
  scale_color_manual(values = c("#8AB17D","#F6D27A"))+
  scale_fill_manual(values = c("#8AB17D","#F6D27A"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("CIN2", "CIN3")), scales = "free_y")+
  theme_bw()+
  ggtitle('E3 infected mice post infection - disease severity, top driving Taxa')


E3.Ist.Meta <- subset.data.frame(E3.Meta, Infected %in% c("Infected"))
E3.I.Meta <- subset.data.frame(E3.Ist.Meta, InfectStage2 =="Infected")
Micro.I.df <- subset.data.frame(Micro.df, row.names(Micro.df) %in% row.names(E3.I.Meta))

X.Ist <-Micro.I.df
Y.Ist <- E3.I.Meta$InfectStageDEM
Y.Ist<-factor(Y.Ist, levels = c("Infection Establishment (2wpi)", "Early Infection (4 & 6wpi)", "Mid Infection (8-10-12wpi)"))

plsda.Ist <- plsda(X.Ist, Y.Ist, ncomp = 3) # note no SCC in this group
plotIndiv(plsda.Ist, 
          ind.names = F, 
          ellipse = TRUE,
          legend=TRUE, cex=6,
          comp= c(1,2),
          #pch = c(3,4,5,2,1),
          col =c("#268A86","#13546F","#053061" ), #"#DE7350"), 
          #xlim = c(-5,5), ylim = c(-5,5),
          title = 'E3 infection stage - infected mice only')

cord1_Load <-plotLoadings(plsda.Ist, comp = 1, contrib = 'max', method = 'median')
cord1_Load <- cord1_Load %>% rename("Variate.1.importance" = "importance")
Cord_1_Loadings <- tibble::rownames_to_column(cord1_Load, "Variable")
TiesX3 <- filter(Cord_1_Loadings,  Contrib.Infection.Establishment..2wpi.=="TRUE"& Contrib.Early.Infection..4...6wpi. == "TRUE" & Contrib.Mid.Infection..8.10.12wpi.=="TRUE") # removing variables that contributed to ALL the groups
# "Normal/Hyperplasia", "CIN 3", "At least CIN3", "SCC"
Tie3Var<-pull(TiesX3, Variable)
Cord_1_filt.df<-Cord_1_Loadings[ !Cord_1_Loadings$Variable %in% Tie3Var,]
Cord_1_filt.df<- Cord_1_filt.df%>% dplyr::select(c("Variable","Contrib.Pre.Infection","Contrib.Post.Depo", "Contrib.Infection.Establishment..2wpi.", "Contrib.Early.Infection..4...6wpi.","Contrib.Mid.Infection..8.10.12wpi.",  "GroupContrib", "Variate.1.importance"))

cord2_Load<-plotLoadings(plsda.Ist, comp = 2, contrib = 'max', method = 'median')
cord2_Load <- cord2_Load %>% rename("Variate.2.importance" = "importance")
Cord_2_Loadings <- tibble::rownames_to_column(cord2_Load, "Variable")
TiesX3 <- filter(Cord_2_Loadings, Contrib.Infection.Establishment..2wpi.=="TRUE"& Contrib.Early.Infection..4...6wpi. == "TRUE" & Contrib.Mid.Infection..8.10.12wpi.=="TRUE") # removing variables that contributed to ALL the groups
# removing variables that contributed to ALL the groups
Tie3Var<-pull(TiesX3, Variable)
Cord_2_filt.df<-Cord_2_Loadings[ !Cord_2_Loadings$Variable %in% Tie3Var,]
Cord_2_filt.df<- Cord_2_filt.df%>% dplyr::select(c("Variable","Contrib.Pre.Infection","Contrib.Post.Depo", "Contrib.Infection.Establishment..2wpi.", "Contrib.Early.Infection..4...6wpi.","Contrib.Mid.Infection..8.10.12wpi.",  "GroupContrib", "Variate.2.importance"))

Cord_Loadings_dataframe <- merge(x=Cord_1_filt.df,y=Cord_2_filt.df, by=c("Variable","GroupContrib"))

Cord_Loadings_dataframe$TotalImportance = with(Cord_Loadings_dataframe, sqrt((Variate.1.importance^2)+(Variate.2.importance^2)))

TopCord <-subset.data.frame(Cord_Loadings_dataframe, TotalImportance > 0.05)
ggplot(TopCord, aes(x = TotalImportance, y = reorder(Variable, TotalImportance), color = GroupContrib, fill = GroupContrib))+
  geom_point(aes(size = 1.5*TotalImportance))+
  geom_col(width = .04)+
  xlim(0,0.3)+
  scale_color_manual(values = c("#13546F","#268A86","#053061", "#EFB366"))+
  scale_fill_manual(values = c("#13546F","#268A86","#053061", "#EFB366"))+
  #facet_wrap(~VariateAxis, scales = "free_y")+
  facet_wrap(~factor(GroupContrib, levels = c("Infection Establishment (2wpi)", "Early Infection (4 & 6wpi)", "Mid Infection (8-10-12wpi)")), 
             scales = "free_y", nrow =1)+
  theme_bw()+
  ggtitle('E3 Infection Stage- top driving Taxa')

# MAASLIN -Experiment 2 - for  figure 2F
# Running Maaslin2 for determining differential abundance - similar to above but theretically more accurate 
# https://huttenhower.sph.harvard.edu/maaslin/
E3BaseRelAbundance_Genus<-read.csv("./E3BaseRelAbundance_Genus.csv")
input_genus_data <- as.data.frame(E3BaseRelAbundance_Genus)
rownames(input_genus_data) <-input_genus_data[,1]

input_meta_file <- Base.meta
input_meta_file$TimeP<- as.numeric(input_meta_file$TimeP)
input_meta_file$ABX.Infect <- factor(input_meta_file$ABX.Infect, levels = c("None_Untouched", "None_Mock","None_Infected", "ABX_Untouched", "ABX_Mock", "ABX_Infected"))

ABXtreated <- subset.data.frame(input_meta_file, AntibioticTx == "ABX")
NoAntibiotics <- subset.data.frame(input_meta_file, AntibioticTx != "ABX")
Post_infection.MF <-subset.data.frame(NoAntibiotics, InfectStagePreUM %in% c("Untouched", "Mock", "Infected"))

# Evaluating the effects of the virus in the group NOT treated with antibiotics # Experiment 2 supplemental figure 1
E3.Natural.MF <- subset.data.frame(NoAntibiotics, NaturalMicro %in% c("Natural", "Post-Depo"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = NoAntibiotics, 
  output = "./MAASLIN/E3_InfectionALLtimepoints_Genus_maslin", 
  fixed_effects = c("Infected"),
  reference = c("Infected", "Infected"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = NoAntibiotics, 
  output = "./MAASLIN/E3_InfectionALLtimepoints_Untouched_Genus_maslin", 
  fixed_effects = c("Infected"),
  reference = c("Infected", "Untouched"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Post_infection.MF, 
  output = "./MAASLIN/E3_Infection_PostInfectionTimepoints_Untouched_Genus_maslin", 
  fixed_effects = c("Infected"),
  reference = c("Infected", "Untouched"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Post_infection.MF, 
  output = "./MAASLIN/E3_Infection_PostInfectionTimepoints_Mock_Genus_maslin", 
  fixed_effects = c("Infected"),
  reference = c("Infected", "Mock"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Post_infection.MF, 
  output = "./MAASLIN/PostinfectionTimepoints_IUM_Infected_Genus", 
  fixed_effects = c("Infected"),
  reference = c("Infected", "Infected"),
  random_effects = c("MP"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Post_infection.MF, 
  output = "./MAASLIN/PostinfectionTimepoints_IUM_Untouched_Genus", 
  fixed_effects = c("Infected"),
  reference = c("Infected", "Untouched"),
  random_effects = c("MP"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = Post_infection.MF, 
  output = "./MAASLIN/PostinfectionTimepoints_IUM_Mock_Genus", 
  fixed_effects = c("Infected"),
  reference = c("Infected", "Mock"),
  random_effects = c("MP"))

NoABX.infect.only <- subset.data.frame(NoAntibiotics, Infected == "Infected")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = NoABX.infect.only, 
  output = "./MAASLIN/E3_INfectedMice_InfectStageDEM_Pre-infection_Genus_maslin", 
  fixed_effects = c("InfectStageDEM"),
  reference = c("InfectStageDEM", "Pre-Infection"),
  random_effects = c("MP"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = NoABX.infect.only, 
  output = "./MAASLIN/E3_INfectedMice_InfectStageDEM_Post-Depo_Genus_maslin", 
  fixed_effects = c("InfectStageDEM"),
  reference = c("InfectStageDEM", "Post-Depo"),
  random_effects = c("MP"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = NoABX.infect.only, 
  output = "./MAASLIN/E3_INfectedMice_InfectStageDEM_Establisment_Genus_maslin", 
  fixed_effects = c("InfectStageDEM"),
  reference = c("InfectStageDEM", "Infection Establishment (2wpi)"),
  random_effects = c("MP"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = NoABX.infect.only, 
  output = "./MAASLIN/E3_INfectedMice_InfectStageDEM_EarlyInfection_Genus_maslin", 
  fixed_effects = c("InfectStageDEM"),
  reference = c("InfectStageDEM", "Early Infection (4 & 6wpi)"),
  random_effects = c("MP"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = NoABX.infect.only, 
  output = "./MAASLIN/E3_INfectedMice_InfectStageDEM_MidInf_Genus_maslin", 
  fixed_effects = c("InfectStageDEM"),
  reference = c("InfectStageDEM", "Mid Infection (8-10-12wpi)"),
  random_effects = c("MP"))

VC.E3 <-subset.data.frame(Post_infection.MF, ViralCopyNumberHighLow %in% c("Low", "High"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = VC.E3 , 
  output = "./MAASLIN/PostinfectionTimepoints_ViralLoad_Genus", 
  fixed_effects = c("ViralCopyNumberHighLow"),
  reference = c("ViralCopyNumberHighLow", "Low"),
  random_effects = c("MP"))

DisCIN.E3 <-subset.data.frame(Post_infection.MF, DiseaseCINEndpoint  %in% c("CIN2", "CIN3"))
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = DisCIN.E3, 
  output = "./MAASLIN/PostinfectionTimepoints_DiseaseCIN_Genus", 
  fixed_effects = c("DiseaseCINEndpoint"),
  reference = c("DiseaseCINEndpoint", "CIN2"),
  random_effects = c("MP"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = NoABX.infect.only, 
  output = "./MAASLIN/E3_INfectedMice_InfectStage2_Pre-infect_Genus_maslin", 
  fixed_effects = c("InfectStage2"),
  reference = c("InfectStage2", "Pre-Infection"),
  random_effects = c("MP"))

E3_PostInfect <- subset.data.frame(NoABX.infect.only, InfectStage2 == "Infected")
fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = E3_PostInfect, 
  output = "./MAASLIN/E3_POSTINfectedMice_InfectStageDEM_Establisment_Genus_maslin", 
  fixed_effects = c("InfectStageDEM"),
  reference = c("InfectStageDEM", "Infection Establishment (2wpi)"),
  random_effects = c("MP"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = E3_PostInfect, 
  output = "./MAASLIN/E3_POSTINfectedMice_InfectStageDEM_EarlyInfection_Genus_maslin", 
  fixed_effects = c("InfectStageDEM"),
  reference = c("InfectStageDEM", "Early Infection (4 & 6wpi)"),
  random_effects = c("MP"))

fit_data = Maaslin2(
  input_data = input_genus_data, 
  input_metadata = E3_PostInfect, 
  output = "./MAASLIN/E3_POSTINfectedMice_InfectStageDEM_MidInf_Genus_maslin", 
  fixed_effects = c("InfectStageDEM"),
  reference = c("InfectStageDEM", "Mid Infection (8-10-12wpi)"),
  random_effects = c("MP"))

