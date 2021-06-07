library(tidyverse)
library(openxlsx)
library(patchwork)

## Data load
# Heme completeness
heme <- read.xlsx("FigS11_Completeness_Heme_vs_Genome.xlsx", sheet="HemeCompleteness") %>% 
  as_tibble()

# Load GTDB data & 
# Rename column names
# Join with heme data
# Make "Category2" column: Isolate vs. MAG/SAG
# Mutate "Category" and "Category2" as factor
# Remove the 9 genomes with <50% completeness
genome_heme <- read.xlsx("FigS11_Completeness_Heme_vs_Genome.xlsx", sheet="GTDB_Completeness_Category") %>% 
  as_tibble() %>%
  rename(Genome = accession, GenomeCompleteness = checkm_completeness, Category = ncbi_genome_category) %>% 
  left_join(heme, by="Genome") %>% 
  mutate(Category2 = if_else(Category=="Isolate", "Isolate", "MAG/SAG", missing = NULL)) %>%
  mutate(Category = factor(Category)) %>% 
  mutate(Category2 = factor(Category2)) %>% 
  filter(GenomeCompleteness >= 50)

## Mann-Whitney test
#different? (two-sided)
wilcox.test(GenomeCompleteness ~ Category2, data=genome_heme)
wilcox.test(HemeCompleteness ~ Category2, data=genome_heme)
# Check which is greater (one-sided)
levels(genome_heme$Category2)
# Note the order -> Isolate (former) is compared to MAG/SAG (latter)
# Designate alternative: "greater"
wilcox.test(GenomeCompleteness ~ Category2, data=genome_heme, alternative="greater")
wilcox.test(HemeCompleteness ~ Category2, data=genome_heme, alternative="greater")


## Violin plots of genome/heme completeness: Isolate vs. MAG/SAG
# Genome
ggplot(genome_heme, aes(x=Category2, y=GenomeCompleteness, fill=Category2)) +
  geom_violin(scale="area")
# Heme  
ggplot(genome_heme, aes(x=Category2, y=HemeCompleteness, fill=Category2)) +
  geom_violin(scale="area")


## Heme completeness per genome completeness interval (MAG/SAG only)
## Not suitable for isolate genomes, because nearly all isolate genomes are >=90% complete
# Cut genome completeness with 10% interval and save into a new column (GenomeCompBin)
msag <- genome_heme %>% 
  filter(!Category == "Isolate") %>% 
  mutate(GenomeCompBin=cut_width(GenomeCompleteness, width=10, boundary=0, closed="left")) 

# Violin plot of heme completeness
# viridis color scale (discrete) for genome completeness
# Set scale="count" to indicate the number of genomes in each genome completeness bin
ggplot(msag, aes(x=GenomeCompBin, y=HemeCompleteness, fill=GenomeCompBin)) +
  geom_violin(scale = "count", alpha=0.7) +
  scale_fill_viridis_d(direction = 1) +
  geom_jitter(alpha=0.1, size=0.2, width=0.1)


## Merge plots using patchwork
p1 <- ggplot(genome_heme, aes(x=Category2, y=GenomeCompleteness, fill=Category2)) +
  geom_violin(scale="area") +
  theme(legend.position = "none")
  
p2 <- ggplot(genome_heme, aes(x=Category2, y=HemeCompleteness, fill=Category2)) +
  geom_violin(scale="area")

p3 <- ggplot(msag, aes(x=GenomeCompBin, y=HemeCompleteness, fill=GenomeCompBin)) +
  geom_violin(scale = "count", alpha=0.7) +
  scale_fill_viridis_d(direction = 1) +
  geom_jitter(alpha=0.1, size=0.2, width=0.1)


(p1 + p2) / p3

