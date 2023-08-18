library(tidyverse)
library(viridis)
library(tidytext)
mydata_old <-read_delim("OutputTables_OCTOBER_16_Magda_analyses.txt", 
           col_names = T, delim = "\t")
mydata_old %>% #distinct(SampleID) %>% filter(is.na(SampleID))
  mutate(SampleName=factor(SampleName), Population=factor(Population)) %>%
  filter(!is.na(HaplotypeSum)) %>% 
  arrange(Population, SampleName, Locus) %>%
  group_by(Population, SampleName, Locus) %>% 
  summarise(`Avg. Read Depth` = mean(HaplotypeSum), Population=Population)%>%
  ggplot(aes(x=Locus,y=SampleName, fill=`Avg. Read Depth`)) + geom_tile() +
  scale_fill_viridis() + facet_grid(Population~., scales = "free") +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x="Locus",y="Sample name") +
  theme(axis.text.y = element_text(size = 7.5),
        axis.text.x = element_text(size = 9)
        )


mydata_old %>% #distinct(SampleID) %>% filter(is.na(SampleID))
  mutate(SampleName=factor(SampleName), Population=factor(Population)) %>%
  filter(!is.na(HaplotypeSum)) %>% 
  arrange(Population, SampleName, Locus) %>%
  group_by(Population, SampleName, Locus) %>% 
  summarise(`Avg. Read Depth` = mean(HaplotypeSum), Population=Population)%>%
  filter(Population=="AFA") %>%
  ggplot(aes(x=Locus,y=SampleName, fill=`Avg. Read Depth`)) + geom_tile(aes(height=1.5)) +
  scale_fill_viridis() + facet_grid(Population~., scales = "free", space="free") +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x="Locus",y="Sample name")


(mydata_old %>% group_by(SampleName, Locus) %>% 
  summarise(n=n()))

# HB across all loci boxplot
mydata_old %>% 
  mutate(SampleName=factor(SampleName), Population=factor(Population)) %>%
  filter(!is.na(HB),HB!="-") %>%
  group_by(Population, SampleName,Locus) %>%
  distinct(HB) %>% 
  ggplot(aes(x=factor(Locus), y=as.numeric(HB), fill=Locus)) + 
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot() + coord_flip() + theme(legend.position = "none") +
  labs(x="Locus",y="Heterozygote balance")
  #geom_point(size=1, alpha=0.19, position = position_jitter(seed=1,width=0.2)) +
  theme(axis.text.x = element_text(angle = 90))


# HB across all loci and population boxplot
  mydata_old %>% 
    mutate(SampleName=factor(SampleName), Population=factor(Population)) %>%
    filter(!is.na(HB),HB!="-") %>%
    group_by(Population, SampleName,Locus) %>%
    distinct(HB) %>% 
    ggplot(aes(x=factor(Locus), y=as.numeric(HB), fill=Population)) + 
    stat_boxplot(geom = "errorbar", width = 0.5) +
    geom_boxplot() + theme(legend.position = "none") +
    labs(x="Locus",y="Heterozygote balance (log10 scale)") + 
    geom_hline(yintercept = c(0.5,2), linetype="dashed", color="red",size=0.75) +
    scale_y_continuous(trans = "log10") +
    facet_wrap(~Population) + coord_flip()  #-> g
#ggplotly(g)
  
    
#working on loci that have off HB with jitter points
  mydata_old %>% 
    mutate(SampleName=factor(SampleName), Population=factor(Population)) %>%
    filter(!is.na(HB),HB!="-") %>%
    group_by(Population, SampleName,Locus) %>%
    distinct(HB) %>% 
    filter(Locus=="rs907100"|Locus=="rs430046"|
             Locus=="rs740598"|Locus=="rs445251"|
             Locus=="rs1490413"|Locus=="rs1015250") %>%
    ggplot(aes(x=factor(Locus), y=as.numeric(HB)#, fill=Population
               )) + 
    #stat_boxplot(geom = "errorbar", width = 0.5) +
    geom_boxplot()  + #theme(legend.position = "none")
    geom_point(aes(color=Population), size=2, alpha=0.5, 
               position = position_jitter(seed=1,width=0.2), 
               ) +
    labs(x="Locus",y="Heterozygote balance") + coord_flip()
  
#working on loci that have off HB, grouped boxplot..dodge population
  mydata_old %>% 
    mutate(SampleName=factor(SampleName), Population=factor(Population)) %>%
    filter(!is.na(HB),HB!="-") %>%
    group_by(Population, SampleName,Locus) %>%
    distinct(HB) %>% 
    filter(Locus=="rs907100"|Locus=="rs430046"|
             #Locus=="rs740598"|Locus=="rs445251"|
             Locus=="rs1490413"|Locus=="rs1015250") %>%
    ggplot(aes(x=factor(Locus), y=as.numeric(HB), fill=Population
    )) + scale_y_continuous(breaks = seq(0,10,1)) + 
    geom_hline(yintercept = c(0.5,2), linetype="dashed", color="red", size=0.75) +
    stat_boxplot(geom = "errorbar") + #geom_text(aes())
    geom_boxplot() + coord_flip() + #theme(legend.position = "none")
    #geom_point(aes(color=Population), size=2, alpha=0.5, 
    #           position = position_jitter(seed=1,width=0.2), 
    #) +
    labs(x="Locus",y="Heterozygote balance")
  

  
# mydata_old %>% group_by(Population, SampleName, Locus) %>%
#   distinct(HB) %>%
#   summarise(n())

# HB across 3 populations boxplot
  mydata_old %>% 
  mutate(SampleName=factor(SampleName), Population=factor(Population)) %>%
  filter(!is.na(HB),HB!="-") %>%
  group_by(Population, SampleName,Locus) %>%
  distinct(HB) %>% 
  ggplot(aes(x=factor(Population), y=as.numeric(HB), color=Population)) + 
  geom_boxplot() +
  geom_point(size=1, alpha=0.19, 
             position = position_jitter(seed=1,width=0.2)
             ) +
  labs(x="Population",y="Heterozygote balance")
  #+ theme(legend.position = "none")
  #theme(axis.text.x = element_text(angle = 90))


# Allele freq
myaf <- read_delim("AlleleFreq.txt", col_names = T, delim = "\t")
myaf_w <- myaf %>%
  pivot_longer(-c(Locus,Allele),
                       names_to = "Population",
                       values_to = "AlleleFrequency")
myaf_w %>% #mutate()
  ggplot(aes(x=AlleleFrequency,y=Locus, fill=Allele)) + 
  geom_bar(stat="identity",posiion=position_dodge2(preserve = "single")) + 
  facet_wrap(~Population) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x="Allele Frequency",y="Locus") + scale_colour_hue(l = 45)#scale_fill_viridis(discrete = T)#-> g
ggplotly(g)


# heterozygosity
myhet_oe <- read_delim("Heterozygosity_Obs_Expec.txt", col_names = T, delim = "\t")
myhet_oe_long <- myhet_oe %>%
  pivot_longer(-c(Locus),names_to = "Population",values_to = "AlleleFreq")

myhet_oe_long %>% separate(Population, into = c("Population","ObsExp"), sep = "_") -> myhet_oe_long_1
myhet_oe_long_1 %>% group_by(Population) %>% arrange(AlleleFreq) %>% 
  #group_by(Locus) %>% 
  mutate(Heterozygosity=AlleleFreq,Locus=reorder_within(Locus, -AlleleFreq, Population),
         ObsExp=ifelse(ObsExp=="E","Expected","Observed"),
         `Observed / Expected` = ObsExp
  ) %>%
  ggplot(aes(fill=`Observed / Expected`)) + 
  geom_bar(aes(x=Heterozygosity,y=Locus),stat = "identity", postion=position_dodge2(preserve = "single")) +
  facet_wrap(~Population, scales = "free") 

# pivot_wider to get exp - obs
myhet_oe_long_1 %>%
  pivot_wider(names_from = ObsExp, values_from = AlleleFreq) -> myhet_oe_long_2
myhet_oe_long_2 %>% ggplot(aes(O,E, color=Population)) + 
  stat_qq()

# distance of obs - exp bar plots
myhet_oe_long_2 %>% mutate(Distance=(E-O)) %>%
  ggplot(aes(x=Locus,y=Distance, fill=Population)) + geom_bar(stat="identity") +
  facet_wrap(~Population, scales = "free") + coord_flip() + theme(legend.position = "none")

# distance of obs -exp, arranged bar plots
myhet_oe_long_2 %>% mutate(Distance=(E-O)) %>% mutate(Locus=reorder_within(Locus, Distance,Population)) %>%
  ggplot(aes(x=Locus,y=Distance, fill=Population)) + geom_bar(stat="identity") +
  facet_wrap(~Population, scales = "free") + coord_flip() + theme(legend.position = "none")



# myaf_w %>% #mutate()
#   ggplot(aes(x=AlleleFrequency,y=Locus, fill=Allele)) + 
#   geom_point(aes(color=Allele))
#   
#   geom_bar(stat="identity",posiion=position_dodge2(preserve = "single")) + 
#   facet_wrap(~Population) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   labs(x="Allele Frequency",y="Locus")

library(gapminder)
df = gapminder %>%
  filter(year %in% c(1952,2007)) %>%
  filter(continent %in% c("Asia")) %>%
  select(country,year,lifeExp, gdpPercap)%>%
  mutate(paired = rep(1:(n()/2),each=2),
         year=factor(year))

# dumbell plot
myhet_oe_long_1 %>% group_by(Population) %>% mutate(paired=rep(1:(n()/2),each=2),
                                                    ObsExp=ifelse(ObsExp=="E","Expected","Observed"),
                                                    `Observed / Expected` = ObsExp,
                                                    Heterozygosity = AlleleFreq,
                                                    Locus=reorder_within(Locus, Heterozygosity, Population)
                                                    ) %>%
  ggplot(aes(x=Heterozygosity,y=Locus)) + 
  geom_line(aes(group=paired)) + geom_point(aes(color=`Observed / Expected`), size=2) +
  facet_wrap(~Population, scales = "free") #-> g
ggplotly(g)
