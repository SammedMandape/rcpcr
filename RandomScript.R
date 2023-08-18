popdf <-  read_delim("PopPaper/OutputTables.tsv", "\t", 
                     escape_double = FALSE, 
                     col_types = cols(
                       Population = col_character(), 
                       SampleName = col_character(), 
                       SampleID = col_character(), 
                       Locus = col_character(), 
                       Nomenclature = col_character(), 
                       Allele = col_character(), 
                       Haplotype = col_character(), 
                       HaplotypeSum = col_number(), 
                       RAP = col_number(), 
                       HB = col_number()), 
                     trim_ws = TRUE)

popdfbackup <- popdf

CombinedAlleleDepth <- popdf %>% 
  group_by(SampleName, Locus) %>% 
  summarise(AD = sum(HaplotypeSum)) %>% 
  ungroup()

popdf <- left_join(popdf, CombinedAlleleDepth, by = c("SampleName", "Locus")) %>% 
  group_by(SampleName, Locus) %>% 
  mutate(
    AlleleRank = rank(Allele, ties.method = "first"),
    Interp = if_else(is.na(HB), "Hom", "Het")
  ) %>% 
  ungroup()


popdf %>%
  filter(AlleleRank == 1, AD >= 1) -> foobar

foobar %>% group_by(Locus) %>%
  summarise(mean.AD = mean(AD, na.rm = TRUE),
            sd.AD = sd(AD, na.rm = TRUE),
            n.AD = n()) %>%
  ungroup() %>%
  mutate(se.AD = sd.AD / sqrt(n.AD),
         lower.ci.AD = mean.AD - qt(1 - (0.05 / 2), n.AD - 1) * se.AD,
         upper.ci.AD = mean.AD + qt(1 - (0.05 / 2), n.AD - 1) * se.AD,
         cv.AD = sd.AD / mean.AD) -> foobar1


library(tidytext)
library(babynames)
install.packages("babynames")
library(tidyverse)

top_names <- babynames %>%
  filter(year >= 1950,
         year < 1990) %>%
  mutate(decade = (year %/% 10) * 10) %>%
  group_by(decade) %>%
  count(name, wt = n, sort = TRUE) %>%
  ungroup() %>%
  mutate(n = ifelse(decade == "1950" & name == "Michael", 846042, n))

top_names %>%
  group_by(decade) %>%
  top_n(15) %>%
  ungroup %>%
  mutate(decade = as.factor(decade),
         name = reorder_within(name, -n, decade)) %>%
  arrange(decade, desc(n), desc(name)) %>% 
  mutate(name = fct_inorder(name)) %>% 
  ggplot(aes(name, n, fill = decade)) +
  geom_col(show.legend = FALSE) +
  #geom_text(aes(label = n), hjust = 1) +
  facet_wrap(~decade, scales = "free_y") +
  coord_flip() +
  #scale_x_reordered() +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "Number of babies per decade",
       x = NULL,
       title = "What were the most common baby names in each decade?",
       subtitle = "Via US Social Security Administration")
  
ggplot(mpg,aes(class, hwy)) + 
  geom_dotplot(binaxis = "y",stackdir = "center", width = 0.1)
#geom_point()
geom_dotp


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
  theme(legend.position="bottom")

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

ggplot(here, aes(x=reorder(Locus, mean.HB), y=mean.HB)) + 
  geom_bar(stat = "identity", fill = "#440154") + 
  geom_errorbar(color = "#5ec962", stat = "identity",position = "identity",
                aes(ymin = mean.HB, ymax = mean.HB+se.HB), size=1) + 
  theme_classic() + 
  labs(title = "Mean Heterozygote Balance", 
       subtitle = "by Locus", x = "Locus", y = "Heterozygote Balance") + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) #+ facet_grid(Population~.)
