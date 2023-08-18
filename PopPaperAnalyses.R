library(tidyverse)
library(viridis)
library(tidytext)
library(data.table)
library(ggplot2)
library(plotly)

mydat <- fread("PopPaper/OutputTables.tsv", sep = "\t",
               header = T) 
mydat %>% mutate(SampleName=factor(SampleName), Population=factor(Population)) %>%
  arrange(Population, SampleName, Locus) %>%
  group_by(Population, SampleName, Locus) %>%
  summarise(`Avg. Read Depth (Sample)` = mean(HaplotypeSum)) %>% #->foo
  ggplot(aes(x=Locus,y=SampleName, fill=`Avg. Read Depth (Sample)`, height=1.5)) + geom_tile(
    #color = "white",
    #lwd = 0.05,
    #linetype = 1
  ) +
  scale_fill_viridis() + facet_grid(Population~., scales = "free") +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x="Locus",y="Sample name") +
  theme(axis.text.y = element_text(size = 7.5),
        axis.text.x = element_text(size = 9)
  )

# Avg Read Depth by locus
mydat %>% 
  group_by(SampleName, Locus) %>% 
  summarise(meanAD = sum(HaplotypeSum)) -> mydatSum
mydatSum %>% 
  ungroup() %>% 
  group_by(Locus) %>% rename("AD" = "meanAD") %>%
  summarise(mean.AD = mean(AD,na.rm = TRUE),
            sd.AD = sd(AD, na.rm = TRUE),
            n.AD = n()) %>%
  ungroup() %>%
  mutate(se.AD = sd.AD / sqrt(n.AD),
         lower.ci.AD = mean.AD - qt(1 - (0.05 / 2), n.AD - 1) * se.AD,
         upper.ci.AD = mean.AD + qt(1 - (0.05 / 2), n.AD - 1) * se.AD,
         cv.AD = sd.AD / mean.AD,
         Med = as.numeric(median(mean.AD)), 
         Mean = mean(mean.AD)) -> mydatAD  

mydatAD %>% summarise(median(mean.AD))
ggplot(mydatAD, aes(reorder(Locus, mean.AD), mean.AD)) + 
  geom_bar(stat = "identity", fill = "#3b528b", width = 0.7) + 
  geom_errorbar(stat = "identity", width=0.8,
                color = "#35b779", size=0.4,
                aes(ymin = mean.AD, 
                    ymax = upper.ci.AD)) + 
  scale_y_continuous(breaks = c(0,2000,4000,6000,8000,10000,12000,14000)) + 
  geom_hline(yintercept = 2664.813, size=0.3, linetype='dashed', color='red') +
  theme_classic() + 
  theme(panel.grid.major.y = element_line(color = "grey", size = 0.1)) +
  labs(#title = "Mean Read Depth", 
       #subtitle = "by Locus", 
       x = "Locus", 
       y = "Read Depth") + 
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))
scale_color_viridis(discrete=TRUE)



HetLoci <- popdf %>% 
  filter(AlleleRank == 1, HB >= 0)

hbWrapup <- HetLoci %>%
  group_by(Locus, Population) %>%
  summarise(mean.HB = mean(HB, na.rm = TRUE),
            sd.HB = sd(HB, na.rm = TRUE),
            n.HB = n()) %>%
  ungroup() %>%
  mutate(se.HB = sd.HB / sqrt(n.HB),
         lower.ci.HB = mean.HB - qt(1 - (0.05 / 2), n.HB - 1) * se.HB,
         upper.ci.HB = mean.HB + qt(1 - (0.05 / 2), n.HB - 1) * se.HB, 
         cv.HB = sd.HB / mean.HB)

#Calculate allele frequencies
alleleFreqbyPop <- haps %>% 
  group_by(Population, Locus, Allele) %>% 
  summarize(alleleCount = n()) %>% 
  ungroup() %>% 
  left_join(popSum1, by = c("Population", "Locus")) %>% 
  mutate(alleleFreq = alleleCount/total) %>% 
  select(-alleleCount)


# heterozygosity
myhet_oe <- read_delim("/home/snm0205/nimagen/PopPaper/Heterozygosity_obs_expec_05242023.txt", col_names = T, delim = "\t")
myhet_oe_long <- myhet_oe %>%
  pivot_longer(-c(Locus),names_to = "Population",values_to = "AlleleFreq")
myhet_oe_long %>% separate(Population, into = c("Population","ObsExp"), sep = "_") -> myhet_oe_long_1

myhet_oe_long_1 %>% 
  group_by(Population) %>% 
  arrange(AlleleFreq) %>% 
  #group_by(Locus) %>% 
  mutate(Heterozygosity=AlleleFreq,Locus=reorder_within(Locus, -AlleleFreq, Population),
         ObsExp=ifelse(ObsExp=="E","Expected","Observed"),
         `Observed / Expected` = ObsExp
  ) %>%
  ggplot(aes(fill=`Observed / Expected`)) + scale_fill_viridis_d(begin=0.8, end =0.3, option = "A",)+
  geom_bar(aes(y=Heterozygosity,x=Locus),stat = "identity", postion=position_dodge2(preserve = "single"), width=0.8) +
  scale_x_reordered() +
  coord_flip() +
  facet_wrap(~Population, scales = "free") + theme(legend.position="bottom") +
  theme_classic() + theme(legend.position="bottom")

# pivot_wider to get exp - obs
myhet_oe_long_1 %>%
  pivot_wider(names_from = ObsExp, values_from = AlleleFreq) -> myhet_oe_long_2
#myhet_oe_long_2 %>% ggplot(aes(O,E, color=Population)) + 
#  stat_qq()

# distance of obs -exp, arranged bar plots
myhet_oe_long_2 %>% mutate(Distance=(E-O)) %>% mutate(Locus=reorder_within(Locus, Distance,Population)) %>%
  ggplot(aes(x=Locus,y=Distance)) + geom_bar(stat="identity", width=0.7) +
  #scale_fill_viridis_d(begin=0.8, end =0.3, option = "A",) +
  facet_wrap(~Population, scales = "free") + coord_flip() + 
  theme_classic() + theme(legend.position = "none") +
  scale_fill_grey() +
  scale_x_reordered()

# Ae
myAe <- fread("/home/snm0205/nimagen/PopPaper/Ae_05242023.txt", sep = "\t",header = T)
myAe %>% 
  pivot_longer(cols = -c(Locus), names_to = "Pop", values_to = "Ae") %>%
  group_by(Locus) %>%
  arrange(desc(Ae), .by_group = TRUE) -> myAe_intermed

myAe_intermed %>% group_by(Pop) %>%
  #summarise(MeanAe = mean(Ae)) %>%
  count(Ae >=2.0)


myAe_intermed %>%
  ggplot(aes(x=reorder(Locus,desc(Ae)),y=Ae, group=Pop)) +
  geom_line(aes(color=Pop),linetype = "dashed", size=0.5) +
  geom_point(aes(color=Pop), size=3) + #ylim(0,3)+
  scale_y_continuous(n.breaks = 6, limits = c(0,3)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5)) + 
  scale_colour_manual(values = c("#f98e09", "#440154", "#5ec962")) +
  #scale_color_viridis(discrete = T, option = "D", begin=0, end = 0.92 ) +
  theme(panel.grid.major.y = element_line(color = "grey", size = 0.1)) +
  labs(x="Locus",y="Effective Number of Alleles (Ae)") +
  theme(legend.position="bottom") #-> p
ggplotly(p)

# AF
myAF <- fread("/home/snm0205/nimagen/PopPaper/alleleFreqbyPop_wide.tsv",
              sep="\t", header = T)
myAF %>%
  select(-c(starts_with("total") )) %>%
  pivot_longer(cols = -c(Locus,Allele), names_to = "Pop", 
               values_to = "AlleleFreq") %>%
  separate(Pop, into = c(NA,"Pop")) %>%
  group_by(Locus, Pop) %>% 
  arrange(desc(AlleleFreq), .by_group = T) %>%
  mutate(Allele_1 = row_number(), 
         Allele_2 = paste0("Allele",Allele_1),
         Allele_3 = ifelse(Allele_2 == "Allele2",AlleleFreq,0)
         ) -> myAF_1

myAF_1 %>%
  ungroup() %>%
  mutate(Locus = reorder_within(Locus, -Allele_3, Pop)) %>%
  ggplot(aes(x=Locus,y=AlleleFreq)) +
  geom_bar(stat = "identity", aes(fill=Allele_2), width = 0.8) + 
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5))  +
  geom_hline(yintercept = 0.1, color="red") +
  facet_wrap(~Pop, scales = "free") + coord_flip() +
  scale_fill_viridis(discrete = T, 
                     option = "D", begin=0.3, end = 0.8,
                     direction = 1) +
  theme_minimal() +
  labs(y='Allele Frequency', fill = "Alleles") +
  scale_x_reordered() +
  theme(legend.position="bottom") + 
  geom_vline(xintercept = 0.05)# -> q
ggplotly(q)

# PD
myPD <- fread("/home/snm0205/nimagen/PopPaper/PD.txt", sep = "\t", header=T)
myPD %>% 
  pivot_longer(cols = -c(Locus), names_to = "Pop", values_to = "PD") %>%
  separate(Pop, into = c(NA,"Pop")) %>% 
  # mutate(Locus = reorder_within(Locus, -PD, Pop)) %>%
  # ggplot(aes(x=Locus, y= PD)) +
  # geom_bar(stat="identity", width = 0.8) + 
  # coord_flip() + 
  # facet_wrap(~Pop, scales = "free") +
  # scale_x_reordered() + 
  ggplot(aes(x=reorder(Locus, desc(PD)),y=PD, group=Pop)) +
  geom_line(aes(color=Pop),linetype = "dashed", size=0.5) +
  geom_point(aes(color=Pop), size=3) + #ylim(0,3)+
  scale_y_continuous(n.breaks = 8, limits = c(0,1)) +
  # theme_minimal() + 
  # #scale_fill_viridis(discrete = T, begin=0, end = 0) +
  # theme(legend.position="none")
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5)) + 
  scale_colour_manual(values = c("#f98e09", "#440154", "#5ec962")) +
  #scale_color_viridis(discrete = T, option = "D", begin=0, end = 0.92 ) +
  theme(panel.grid.major.y = element_line(color = "grey", size = 0.1)) +
  labs(x="Locus",y="Power of Discrimination (PD)") +
  theme(legend.position="bottom") #-> p
ggplotly(p)

# Forensic Param, PIC and PD
myfp_pic <- fread("/home/snm0205/nimagen/PopPaper/ForensicParam.txt", sep="\t")
  
myfp_pic %>%
  pivot_longer(cols = -c(Locus), names_to = "Names", values_to = "PIC") %>%
  separate(Names, into = c("Pop","ForensicParameter")) %>%
  filter(ForensicParameter == "PIC") %>% select(-ForensicParameter) %>%
  ggplot(aes(x=reorder(Locus, desc(PIC)),y=PIC, group=Pop)) +
  geom_line(aes(color=Pop),linetype = "dashed", size=0.5) +
  geom_point(aes(color=Pop), size=3) + #ylim(0,3)+
  scale_y_continuous(n.breaks = 8, limits = c(0,1)) +
  # theme_minimal() + 
  # #scale_fill_viridis(discrete = T, begin=0, end = 0) +
  # theme(legend.position="none")
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5)) + 
  scale_colour_manual(values = c("#f98e09", "#440154", "#5ec962")) +
  #scale_color_viridis(discrete = T, option = "D", begin=0, end = 0.92 ) +
  theme(panel.grid.major.y = element_line(color = "grey", size = 0.1)) +
  labs(x="Locus",y="Polymorphism Information Content (PIC)") +
  theme(legend.position="bottom")
  
# HB
HetLoci <- popdf %>% 
  filter(AlleleRank == 1, HB >= 0)

hbWrapup <- HetLoci %>%
  group_by(Locus, Population) %>% # to get locus and pop level info
  #group_by(Locus) %>% #to get Locus level info
  summarise(mean.HB = mean(HB, na.rm = TRUE),
            sd.HB = sd(HB, na.rm = TRUE),
            n.HB = n()) %>%
  ungroup() %>%
  mutate(se.HB = sd.HB / sqrt(n.HB),
         lower.ci.HB = mean.HB - qt(1 - (0.05 / 2), n.HB - 1) * se.HB,
         upper.ci.HB = mean.HB + qt(1 - (0.05 / 2), n.HB - 1) * se.HB, 
         cv.HB = sd.HB / mean.HB)

here <- hbWrapup %>% filter(Locus == "rs430046"| Locus == "rs1490413" | Locus == "rs1015250"| Locus == "rs907100")

# mean
ggplot(here, aes(x=reorder(Locus, mean.HB), y=mean.HB)) + 
  geom_bar(stat = "identity", fill = "#440154") + 
  geom_errorbar(color = "#5ec962", stat = "identity",position = "identity",
                aes(ymin = mean.HB, ymax = mean.HB+se.HB), size=1) + 
  theme_classic() + 
  labs(title = "Mean Heterozygote Balance", 
       subtitle = "by Locus", x = "Locus", y = "Heterozygote Balance") + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) #+ facet_grid(Population~.)

# by pop
ggplot(here, aes(x=reorder(Locus, mean.HB), y=mean.HB)) + 
  geom_bar(stat = "identity", aes(fill = Population), position = "dodge2") + 
  geom_errorbar(stat = "identity",position = "dodge2",
                aes(ymin = mean.HB, ymax = mean.HB+se.HB)) + 
  scale_fill_manual(values = c("#f98e09", "#3b528b", "#5ec962")) +
  #geom_hline(yintercept = 0.5, color="red", size=0.3, linetype='dashed') +
  #geom_hline(yintercept = 2.0, color="red", size=0.3, linetype='dashed') +
  theme_classic() + 
  labs(title = "Mean Heterozygote Balance", 
       subtitle = "by Locus", x = "Locus", y = "Heterozygote Balance") + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
