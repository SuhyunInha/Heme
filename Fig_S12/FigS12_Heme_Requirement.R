library(tidyverse)
library(openxlsx)
library(patchwork)

## KOtoEC -> convert to long format
KOtoEClong <- read.xlsx("FigS12_Heme_Requirement.xlsx", sheet="KOtoEC") %>% 	
  as_tibble() %>% 	
  pivot_longer(-KO, names_to="OrderinKO", values_to="EC", values_drop_na=T)	%>% 
  select("KO", "EC")
write.xlsx(KOtoEClong, file="KOtoEClong.xlsx")

## Matching EC numbers of heme/cytochrome-requiring enzymes to corresponding KOs
HemeECtoKO <- read.xlsx("FigS12_Heme_Requirement.xlsx", sheet="BRENDA.Curated") %>% 
  as_tibble() %>% 
  select("EC") %>%  
  inner_join(KOtoEClong, by="EC")
write.xlsx(HemeECtoKO, "HemeECtoKO.xlsx")

## Extract unique KOs and convert to the hal file format (for KofamScan)
HemeKO.hal <- HemeECtoKO %>% 
  distinct(KO) %>% 
  mutate(KO = str_replace(KO, "$", ".hmm")) 
write.table(HemeKO.hal, "HemeKO.hal", quote=FALSE, row.names=FALSE, col.names=FALSE)
  

### Scripts for Fig. S12
## Data load and edit
HemeKO_Genome_Count <- read.xlsx("FigS12_Heme_Requirement.xlsx", sheet="HemeKO.Genome.Count") %>% 
  as_tibble() 

hemecc <- read.xlsx("FigS12_Heme_Requirement.xlsx", sheet="Genome.HemeCompleteness") %>% 
  as_tibble() %>% 
  left_join(HemeKO_Genome_Count, by="Genome") %>% 
  mutate(HemeKO = replace_na(HemeKO, "0")) %>%
  # replace NA (genomes with no heme-requiring enzymes) into 0
  mutate(HemeKO = as.double(HemeKO)) %>% 
  mutate(NormHemeKO = HemeKO/Size*1000000) %>% 
  # Add a new column for normalized (by genome size) count of heme-requiring KO
  mutate(Nanopelagicaceae = if_else(Family=="f__Nanopelagicaceae",
                                    "Nanopelagicaceae", "Others", missing = NULL)) %>%
  # Add a new column: Nanopelagicaceae vs. Others
  mutate(Nanopelagicaceae = factor(Nanopelagicaceae)) %>% 
  mutate(HemeCompBin=cut_width(HemeCompleteness, width=50, boundary=0, closed="left")) 
  # Add a new column : Divide all genomes into two groups following heme completeness
  # with a threshold of 50%

# d__Bacteria only
Bacteria <- hemecc %>% 
  filter(Domain=="d__Bacteria")

# p__Actinobacteriota only
Actinobacteriota <- hemecc %>% 
  filter(Phylum=="p__Actinobacteriota")

# c__Actinobacteria only (for R89, maybe c__Actinomycetia from R95)
Actinobacteria <- hemecc %>% 
  filter(Class=="c__Actinobacteria")


## acI (=Nanopelagicaceae) vs. Others (Plot and Mann-Whitney test)

# against other c__Actinobacteria
# Heme requiring KO
ggplot(Actinobacteria, aes(x=Nanopelagicaceae, y=HemeKO, color=Nanopelagicaceae)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = "Class Actinobacteria", y = "Number of KOs requiring heme") +
  theme(axis.title = element_text(vjust=0, size=15)) +
  theme(axis.text = element_text(size = 12, face="bold")) +
  theme(legend.position = "none")

# different? (two-sided)
wilcox.test(HemeKO ~ Nanopelagicaceae, data=Actinobacteria)
# Check which is greater (one-sided)
levels(Actinobacteria$Nanopelagicaceae)
# Note the order -> Nanopelagicaceae (former) is compared to Others (latter)
# Designate alternative: "greater" or "less"
wilcox.test(HemeKO ~ Nanopelagicaceae, data=Actinobacteria, alternative="less")

# Normalized number of KO
ggplot(Actinobacteria, aes(x=Nanopelagicaceae, y=NormHemeKO, color=Nanopelagicaceae)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = "Class Actinobacteria", y = "Normalized number of KOs requiring heme") +
  theme(axis.title = element_text(vjust=0, size=15)) +
  theme(axis.text = element_text(size = 12, face="bold")) +
  theme(legend.position = "none")

wilcox.test(NormHemeKO ~ Nanopelagicaceae, data=Actinobacteria)
wilcox.test(NormHemeKO ~ Nanopelagicaceae, data=Actinobacteria, alternative="greater")

# against other p__Actinobacteriota
# Heme requiring KO
ggplot(Actinobacteriota, aes(x=Nanopelagicaceae, y=HemeKO, color=Nanopelagicaceae)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = "Phylum Actinobacteriota", y = "Number of KOs requiring heme") +
  theme(axis.title = element_text(vjust=0, size=15)) +
  theme(axis.text = element_text(size = 12, face="bold")) +
  theme(legend.position = "none")

wilcox.test(HemeKO ~ Nanopelagicaceae, data=Actinobacteriota)
levels(Actinobacteriota$Nanopelagicaceae)
wilcox.test(HemeKO ~ Nanopelagicaceae, data=Actinobacteriota, alternative="less")

# Normalized number of KO
ggplot(Actinobacteriota, aes(x=Nanopelagicaceae, y=NormHemeKO, color=Nanopelagicaceae)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = "Phylum Actinobacteriota", y = "Normalized number of KOs requiring heme") +
  theme(axis.title = element_text(vjust=0, size=15)) +
  theme(axis.text = element_text(size = 12, face="bold")) +
  theme(legend.position = "none")

wilcox.test(NormHemeKO ~ Nanopelagicaceae, data=Actinobacteriota)
wilcox.test(NormHemeKO ~ Nanopelagicaceae, data=Actinobacteriota, alternative="greater")

# against other d__Bacteria
# Heme requiring KO
ggplot(Bacteria, aes(x=Nanopelagicaceae, y=HemeKO, color=Nanopelagicaceae)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = "All bacteria", y = "Number of KOs requiring heme") +
  theme(axis.title = element_text(vjust=0, size=15)) +
  theme(axis.text = element_text(size = 12, face="bold")) +
  theme(legend.position = "none")

wilcox.test(HemeKO ~ Nanopelagicaceae, data=Bacteria)
levels(Bacteria$Nanopelagicaceae)
wilcox.test(HemeKO ~ Nanopelagicaceae, data=Bacteria, alternative="less")

# Normalized number of KO
ggplot(Bacteria, aes(x=Nanopelagicaceae, y=NormHemeKO, color=Nanopelagicaceae)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = "All bacteria", y = "Normalized number of KOs requiring heme") +
  theme(axis.title = element_text(vjust=0, size=15)) +
  theme(axis.text = element_text(size = 12, face="bold")) +
  theme(legend.position = "none")

wilcox.test(NormHemeKO ~ Nanopelagicaceae, data=Bacteria)
wilcox.test(NormHemeKO ~ Nanopelagicaceae, data=Bacteria, alternative="greater")


# Merge plots
p1 <- ggplot(Actinobacteria, aes(x=Nanopelagicaceae, y=HemeKO, color=Nanopelagicaceae)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

p2 <- ggplot(Actinobacteria, aes(x=Nanopelagicaceae, y=NormHemeKO, color=Nanopelagicaceae)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

p3 <- ggplot(Actinobacteriota, aes(x=Nanopelagicaceae, y=HemeKO, color=Nanopelagicaceae)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "Number of KOs requiring heme") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

p4 <- ggplot(Actinobacteriota, aes(x=Nanopelagicaceae, y=NormHemeKO, color=Nanopelagicaceae)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "Normalized number of KOs requiring heme") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

p5 <- ggplot(Bacteria, aes(x=Nanopelagicaceae, y=HemeKO, color=Nanopelagicaceae)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

p6 <- ggplot(Bacteria, aes(x=Nanopelagicaceae, y=NormHemeKO, color=Nanopelagicaceae)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

(p1 + p2) / (p3 + p4) / (p5 +  p6)


## Partitioning of bacterial and archaeal genomes into two groups according to heme completeness
# Cutoff = 50% (refer to the violin plot of heme completeness for isolates and MAG/SAG)

Archaea <- hemecc %>% 
  filter(Domain=="d__Archaea")

# Bacteria
# Heme requiring KO
ggplot(Bacteria, aes(x=HemeCompBin, y=HemeKO, color=HemeCompBin)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

wilcox.test(HemeKO ~ HemeCompBin, data=Bacteria)
levels(Bacteria$HemeCompBin)
wilcox.test(HemeKO ~ HemeCompBin, data=Bacteria, alternative="less")

# Normalized Heme requiring KO
ggplot(Bacteria, aes(x=HemeCompBin, y=NormHemeKO, color=HemeCompBin)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

wilcox.test(NormHemeKO ~ HemeCompBin, data=Bacteria)
wilcox.test(NormHemeKO ~ HemeCompBin, data=Bacteria, alternative="less")


# Archaea
# Heme requiring KO
ggplot(Archaea, aes(x=HemeCompBin, y=HemeKO, color=HemeCompBin)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

wilcox.test(HemeKO ~ HemeCompBin, data=Archaea)
levels(Archaea$HemeCompBin)
wilcox.test(HemeKO ~ HemeCompBin, data=Archaea, alternative="less")

# Normalized Heme requiring KO
ggplot(Archaea, aes(x=HemeCompBin, y=NormHemeKO, color=HemeCompBin)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

wilcox.test(NormHemeKO ~ HemeCompBin, data=Archaea)
wilcox.test(NormHemeKO ~ HemeCompBin, data=Archaea, alternative="less")


# Merge plots
p7 <- ggplot(Bacteria, aes(x=HemeCompBin, y=HemeKO, color=HemeCompBin)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

p8 <-ggplot(Bacteria, aes(x=HemeCompBin, y=NormHemeKO, color=HemeCompBin)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

p9 <- ggplot(Archaea, aes(x=HemeCompBin, y=HemeKO, color=HemeCompBin)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

p10 <- ggplot(Archaea, aes(x=HemeCompBin, y=NormHemeKO, color=HemeCompBin)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3, width=0.2) +
  labs(x = NULL, y = "") +
  theme(axis.title = element_text(vjust=0, size=12)) +
  theme(axis.text = element_text(size = 10, face="bold")) +
  theme(legend.position = "none")

(p7 + p8) / (p9 + p10)
