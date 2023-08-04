.p <- c("tidyverse", "ggpubr", "gUtils")
suppressPackageStartupMessages(lapply(.p, require, character.only=T))

options(scipen = 999999)
## For CENBP (centromeric) regions;
kmer17 <- read.delim("Data/ChIP-seq/Mitotic_bookmarking/downstream/repeats/jellyfish_kmer/replicates/17-mer_counts.txt", header = F)

kmer17_tidy <- kmer17 %>%
  separate(V1, into= paste0("A",1:4,sep=""), sep= "[-_]") %>%
  separate(V2, into= paste0("B",1:2,sep=""), sep= " ") %>%
  separate(V3, into= paste0("C",1:2,sep=""), sep= " ") %>%
  rename(variant = A1, phase = A2, factor = A3, rep = A4, cenbp =B2, control = C2) %>%
  dplyr::select(-c(B1, C1)) %>%
  mutate(cenbp = as.numeric(cenbp), control = as.numeric(control)) %>%
  mutate(library_normalised_count = cenbp/control) %>% ## library normalisation
  dplyr::select(-c(cenbp, control)) %>%
  mutate(ID = paste0(variant,":",phase,":", rep)) %>%
  dplyr::select(-c(variant, phase, rep)) %>%
  pivot_wider(values_from = library_normalised_count, names_from = factor) %>%
  mutate(sample_library_normalised_count =  ChIP/Input) %>% ## sample normalisation
  separate(ID, into =paste0("A", 1:3, sep=""), sep=":") %>%
  rename(variant=A1, phase=A2, rep=A3) %>%
  dplyr::select(-c(ChIP, Input))

## Extract means from bl6 reads, and normalise to these:
bl6_interphase_mean <- kmer17_tidy %>% filter(variant == "Bl6" & phase == "Interphase") %>% summarise(mean = mean(sample_library_normalised_count)) %>% .$mean
bl6_mitosis_mean <- kmer17_tidy %>% filter(variant == "Bl6" & phase == "Mitotic") %>% summarise(mean = mean(sample_library_normalised_count)) %>% .$mean

kmer17_tidy_bl6_normalised  <- kmer17_tidy %>%
  filter(variant != "Bl6") %>%
  mutate(bl6_sample_library_normalised_count = ifelse(phase == "Interphase", sample_library_normalised_count/bl6_interphase_mean, sample_library_normalised_count/bl6_mitosis_mean))

## Check count distribution for all:
shapiro.test(kmer17_tidy_bl6_normalised$bl6_sample_library_normalised_count)
ggqqplot(kmer17_tidy_bl6_normalised$bl6_sample_library_normalised_count)

shapiro.test(kmer17_tidy_bl6_normalised$bl6_sample_library_normalised_count[kmer17_tidy_bl6_normalised$variant == "WT"])
ggqqplot(kmer17_tidy_bl6_normalised$bl6_sample_library_normalised_count[kmer17_tidy_bl6_normalised$variant == "WT"])
shapiro.test(kmer17_tidy_bl6_normalised$bl6_sample_library_normalised_count[kmer17_tidy_bl6_normalised$variant == "3A"])
ggqqplot(kmer17_tidy_bl6_normalised$bl6_sample_library_normalised_count[kmer17_tidy_bl6_normalised$variant == "3A"])
shapiro.test(kmer17_tidy_bl6_normalised$bl6_sample_library_normalised_count[kmer17_tidy_bl6_normalised$variant == "5MK"])
ggqqplot(kmer17_tidy_bl6_normalised$bl6_sample_library_normalised_count[kmer17_tidy_bl6_normalised$variant == "5MK"])


## Use t-test
x <- kmer17_tidy_bl6_normalised %>%
  mutate(phase = recode(phase, Mitotic = "Mitosis")) %>%
  ggplot(aes(x=phase, y=bl6_sample_library_normalised_count, fill=phase)) +
  geom_boxplot( width=0.4) +
  geom_point(aes(color=variant), size=3) +
  stat_compare_means(method = "t.test", label = "p.format") +
  scale_fill_brewer(name = "Phase", palette = "Dark2") +
  scale_color_brewer(name ="Variant", palette = "Set2") +
  # labs(x= "Phase", y="Normalised counts", title= "CENBP sequence (CTTCGTTGGAAACGGGA)") +
  labs(x= "Phase", y="Normalised counts") + 
  scale_x_discrete(expand=c(0.5,0)) +
  scale_y_continuous(expand =c(0.02,1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12),
        legend.title= element_text(size=13, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size=13))
x

## Split according to variants:
a <- kmer17_tidy_bl6_normalised %>%
  mutate(phase = recode(phase, Mitotic = "Mitosis")) %>%
  ggplot(aes(x=phase, y=bl6_sample_library_normalised_count, fill=phase)) +
  geom_boxplot(width=0.4) +
  geom_point() +
  stat_compare_means(method = "t.test", label = "p.format") +
  facet_grid(~variant) +
  scale_fill_brewer(name = "Phase", palette = "Dark2") +
  # labs(x= "Phase", y="Normalised counts", title= "CENBP sequence (CTTCGTTGGAAACGGGA)") +
  labs(x= "Phase", y="Normalised counts") + 
  scale_x_discrete(expand=c(0.5,0)) +
  scale_y_continuous(expand =c(0.02,1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12),
        legend.title= element_text(size=13, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size=13))
a

## Telomeric tandem repeats:
kmer12 <- read.delim("Data/ChIP-seq/Mitotic_bookmarking/downstream/repeats/jellyfish_kmer/replicates/12-mer_counts.txt", header = F)

kmer12_tidy <- kmer12 %>%
  separate(V1, into= paste0("A",1:4,sep=""), sep= "[-_]") %>%
  separate(V2, into= paste0("B",1:2,sep=""), sep= " ") %>%
  separate(V3, into= paste0("C",1:2,sep=""), sep= " ") %>%
  rename(variant = A1, phase = A2, factor = A3, rep = A4, tel_tandem =B2, control = C2) %>%
  dplyr::select(-c(B1, C1)) %>%
  mutate(tel_tandem = as.numeric(tel_tandem), control = as.numeric(control)) %>%
  mutate(library_normalised_count = tel_tandem/control) %>% ## library normalisation
  dplyr::select(-c(tel_tandem, control)) %>%
  mutate(ID = paste0(variant,":",phase,":", rep)) %>%
  dplyr::select(-c(variant, phase, rep)) %>%
  pivot_wider(values_from = library_normalised_count, names_from = factor) %>%
  mutate(sample_library_normalised_count =  ChIP/Input) %>% ## sample normalisation
  separate(ID, into =paste0("A", 1:3, sep=""), sep=":") %>%
  rename(variant=A1, phase=A2, rep=A3) %>%
  dplyr::select(-c(ChIP, Input))

## Extract means from bl6 reads, and normalise to these:
bl6_interphase_mean <- kmer12_tidy %>% filter(variant == "Bl6" & phase == "Interphase") %>% summarise(mean = mean(sample_library_normalised_count)) %>% .$mean
bl6_mitosis_mean <- kmer12_tidy %>% filter(variant == "Bl6" & phase == "Mitotic") %>% summarise(mean = mean(sample_library_normalised_count)) %>% .$mean

kmer12_tidy_bl6_normalised  <- kmer12_tidy %>%
  filter(variant != "Bl6") %>%
  mutate(bl6_sample_library_normalised_count = ifelse(phase == "Interphase", sample_library_normalised_count/bl6_interphase_mean, sample_library_normalised_count/bl6_mitosis_mean))

## Joint plots of all variants:
shapiro.test(kmer12_tidy_bl6_normalised$bl6_sample_library_normalised_count)
ggqqplot(kmer12_tidy_bl6_normalised$bl6_sample_library_normalised_count)

shapiro.test(kmer12_tidy_bl6_normalised$bl6_sample_library_normalised_count[kmer12_tidy_bl6_normalised$variant == "WT"])
ggqqplot(kmer12_tidy_bl6_normalised$bl6_sample_library_normalised_count[kmer12_tidy_bl6_normalised$variant == "WT"])
shapiro.test(kmer12_tidy_bl6_normalised$bl6_sample_library_normalised_count[kmer12_tidy_bl6_normalised$variant == "3A"])
ggqqplot(kmer12_tidy_bl6_normalised$bl6_sample_library_normalised_count[kmer12_tidy_bl6_normalised$variant == "3A"])
shapiro.test(kmer12_tidy_bl6_normalised$bl6_sample_library_normalised_count[kmer12_tidy_bl6_normalised$variant == "5MK"])
ggqqplot(kmer12_tidy_bl6_normalised$bl6_sample_library_normalised_count[kmer12_tidy_bl6_normalised$variant == "5MK"])

## All:
y <- kmer12_tidy_bl6_normalised %>%
  mutate(phase = recode(phase, Mitotic = "Mitosis")) %>%
  ggplot(aes(x=phase, y=bl6_sample_library_normalised_count, fill=phase)) +
  geom_boxplot( width=0.3) +
  geom_point(aes(color=variant), size=3) +
  stat_compare_means(method = "t.test", label = "p.format") +
  scale_fill_brewer(name = "Phase", palette = "Dark2") +
  scale_color_brewer(name ="Variant", palette = "Set2") +
  # labs(x= "Phase", y="Normalised counts", title= "Telomere tandem sequence (CCCTAACCCTAA)") +
  labs(x= "Phase", y="Normalised counts") + 
  scale_x_discrete(expand=c(0.5,0)) +
  scale_y_continuous(expand =c(0.02,1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12),
        legend.title= element_text(size=13, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size=13)) 
y
 
## Split
b <- kmer12_tidy_bl6_normalised %>%
  mutate(phase = recode(phase, Mitotic = "Mitosis")) %>%
  ggplot(aes(x=phase, y=bl6_sample_library_normalised_count, fill=phase)) +
  geom_boxplot(width=0.4) +
  geom_point() +
  stat_compare_means(method = "t.test", label = "p.format") +
  facet_grid(~variant) +
  scale_fill_brewer(name = "Phase", palette = "Dark2") +
  # labs(x= "Phase", y="Normalised counts", title= "Telomere tandem sequence (CCCTAACCCTAA)") +
  labs(x= "Phase", y="Normalised counts") + 
  scale_x_discrete(expand=c(0.5,0)) +
  scale_y_continuous(expand =c(0.02,1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=12),
        legend.title= element_text(size=13, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size=13))

d <- ggarrange(x,y, ncol=1, nrow=2, align = "hv", common.legend = T, legend = "right")
e <- ggarrange(a,b, ncol=1, nrow=2, align = "hv", common.legend = T, legend = "right")

pdf(paste0("Data/ChIP-seq/Mitotic_bookmarking/downstream/figures/centromeric_seq_jellyfish_all.pdf"), width = 4, height = 4.5, onefile=F)
d <- ggarrange(x,y, ncol=1, nrow=2, align = "hv", common.legend = T, legend = "right")
d
dev.off()

pdf(paste0("Data/ChIP-seq/Mitotic_bookmarking/downstream/figures/centromeric_seq_jellyfish_split_variants.pdf"), width = 7, height = 4.6, onefile=F)
e <- ggarrange(a,b, ncol=1, nrow=2, align = "hv", common.legend = T, legend = "right")
e
dev.off()