library(tidyverse)
library(data.table)
library(plotly)
library(tidytext)
library(treemapify)
library(waffle)

readFile <- function(x){
  tib <- fread(x, sep = "\t", header = F, skip = 2)
  colnames(tib) <- c("leftp", "leftp_seq",
                     "leftp_ham", "tseq",
                     "rightp", "rightp_seq", 
                     "rightp_ham", "readcount")
  
  tib$Filename <- x
  return(tib)
}

sen00<-"Sensitivity/00/*00_results.txt"
sen11<-"Sensitivity/11/*11_results.txt"
sen22<-"Sensitivity/22/*22_results.txt"
sen33<-"Sensitivity/33/*33_results.txt"

files00<-Sys.glob(sen00)
files11<-Sys.glob(sen11)
files22<-Sys.glob(sen22)
files33<-Sys.glob(sen33)

allsen00<-rbindlist(map(files00, readFile))
allsen11<-rbindlist(map(files11, readFile))
allsen22<-rbindlist(map(files22, readFile))
allsen33<-rbindlist(map(files33, readFile))


allsen <- bind_rows(allsen00, allsen11, allsen22, allsen33)

allsen %>%
  separate(Filename, into = c(NA, "Fuzz", "Sample"), sep="/",remove = F) -> bar1

bar1 %>%
  separate(Sample, into = c("Sample1","DNAInput"), sep = "-", remove = F) -> bar2

bar2 %>%
  separate(DNAInput, into = c("DNAInput"), sep = "_") %>%
  mutate(DNAInput = ifelse(is.na(DNAInput), "1ng", DNAInput))-> bar3


bar3 %>% separate(Sample1, into = c("Pop","Sample2"), 
                  sep ="(?<=[a-zA-Z])(?=[0-9])", remove = F) %>%
  separate(Sample2, into = c("Sample2"), sep = "_") %>%
  mutate(Sample3 = ifelse(is.na(Sample2),Pop,Sample2),
         Pop2 = ifelse(is.na(Sample2), "Control", Pop)
  ) %>%
  separate(Sample3, into = c("Sample3"), sep = "_") %>%
  mutate(DNAInput = ifelse(is.na(DNAInput),"1ng",DNAInput)) %>% select(-c(Sample1, Sample2)) %>%
  rename(Sample2=Sample3) %>%
  mutate(Pop2 = case_when(
    Sample2 == "22768" ~ "HIS",
    Sample2 == "21832" ~ "HIS",
    Sample2 == "5544" ~ "AFA",
    Sample2 == "41131" ~ "CAU",
    TRUE ~ Pop2
  )
  )  -> bar4

# allsen %>% 
#   separate(Filename, into = c(NA, "Fuzz", "Sample"), remove = F) -> bar1 #%>%
#   separate(Sample, into = c("Pop","Sample1"), 
#            sep ="(?<=[a-zA-Z])(?=[0-9])", remove = F) %>%
#     mutate(Sample2 = ifelse(is.na(Sample1),Pop,Sample1),
#            Pop2 = ifelse(is.na(Sample1), "Control", Pop)
#            ) %>%
#     select(
#       -Sample1
#     ) %>%
bar4 %>%
  group_by(Fuzz, Pop2, Sample2, DNAInput, leftp, tseq) %>%
    summarise(ReadCTotal = sum(readcount)) %>%
  arrange(desc(ReadCTotal), .by_group = TRUE)-> allsen.1

allsen.1 %>% 
  #group_by(Fuzz, Pop2, DNAInput, Sample2, leftp, tseq) %>%
  top_n(2, ReadCTotal) %>%
  # top_n gives all counts for ties, the following takes care of it or just 
  # use rank instead
  group_by(Fuzz, Pop2, Sample2, DNAInput, leftp, ReadCTotal) %>% summarise(n())-> foo.3
  #filter( leftp =="AAAACCGGAGAGCTGGCGCTGAA") %>%

foo.3 %>% filter(DNAInput=="31pg") -> bar5

primer_in <- fread("Updated_IDseek_SNP85_Primer_Seq_hg38_03232021_SM04082022_primerMod.txt")
primer_in %>% mutate(leftp = ifelse(New_switch==1,`Left primer`,`Right primer`))-> foo

allsen.1 %>%
  left_join(foo, by = "leftp") -> allsen.2

allsen.2 %>% 
  mutate(rank = rank(-ReadCTotal, ties.method = "first")) %>%
  top_n(-2, rank) %>%
  mutate(logReadCTotal = log10(ReadCTotal))-> foo.3


 foo.3 %>% filter((Fuzz=="00" | Fuzz=="11")) %>% 
#                  DNAInput == "31pg", 
#                  (Sample2=="05643" | Sample2=="34125")) %>% 
filter(!(Pop2=="Control" | Pop2=="NTC" )) %>% #-> foo.4 
ggplot(aes(x=reorder(RefSNP, ReadCTotal),y=logReadCTotal)) +
  geom_bar(aes(fill=factor(rank)),stat = "identity") +
  theme(axis.text.x = element_text(angle = 90,vjust=0.45)) +
  #facet_wrap(Sample2~Fuzz)#-> q
  facet_grid(Pop2~Fuzz)#-> q

ggplotly(q)

# take care of drop-out problem by changing fuzzes
foo.3 %>% filter((Fuzz=="00" |Fuzz=="11"), Pop2=="AFA") %>% 
  mutate(DNAInput = factor(DNAInput, levels = c("1ng","500pg","250pg","125pg","62pg","31pg"))) %>%
  #                  DNAInput == "31pg", 
  #                  (Sample2=="05643" | Sample2=="34125")) %>% 
  filter(!(Pop2=="Control" | Pop2=="NTC" )) %>% #-> foo.4 
  ggplot(aes(x=reorder(RefSNP, ReadCTotal),y=logReadCTotal)) +
  geom_bar(aes(fill=factor(rank)),stat = "identity") +
  theme(axis.text.x = element_text(angle = 90,vjust=0.45)) +
  #facet_wrap(Sample2~Fuzz)#-> q
  facet_grid(DNAInput~Fuzz)

bar4 %>% 
  left_join(foo, by = "leftp") -> foo2

foo2 %>%  
filter(Fuzz == "00" | Fuzz == "11") %>%
  filter(RefSNP == "rs735155", DNAInput == "31pg") -> foo4

foo.3 %>% ungroup() %>% filter(Fuzz == "00" | Fuzz=="11") %>%
  group_by(Fuzz, Pop2, Sample2, DNAInput) %>% 
  distinct(RefSNP) %>% 
  summarise(n()) -> foo.4

foo.3 %>%
  filter(Fuzz=="11", Pop2=="AFA", ) %>%
  ggplot(aes(x=reorder(RefSNP, ReadCTotal),y=logReadCTotal)) +
  geom_bar(aes(fill=factor(rank)), stat = "identity", position = "dodge2") +
  theme(axis.text.x = element_text(angle = 90,vjust=0.45)) +
  facet_grid(DNAInput~Fuzz)

# qqplot
foo.3 %>% filter((Fuzz=="00" |Fuzz=="11"), Pop2=="AFA") %>%
  filter(!(Pop2=="Control" | Pop2=="NTC" )) %>%
  ggplot(aes(x=))

foo.4 %>% filter(!(Pop2=="Control" | Pop2=="NTC" )) %>%
  ggplot(aes(x=`n()`, y=Pop2)) +
  geom_bar(stat = "identity", position = "dodge2", aes(fill=Fuzz)) +
  facet_wrap(~DNAInput)

# Recovery of snps at different DNAInput
foo.4 %>% ungroup() %>% filter(!(Pop2=="Control" | Pop2=="NTC" )) %>% rename(count=`n()`) %>%
  group_by(Fuzz, Pop2, DNAInput) %>% summarise(maxi = max(count), mini= min(count)) %>%
  mutate(Count1 = ifelse((Fuzz=="00" & DNAInput=="31pg"), mini,maxi)) %>%
  mutate(DNAInput = factor(DNAInput, levels = c("1ng","500pg","250pg","125pg","62pg","31pg"))) %>%
  ggplot(aes(x=Count1, y=DNAInput)) +
  geom_bar(stat = "identity", position = "dodge2", aes(fill=Fuzz)) +
  facet_wrap(~Pop2)

allsen.2 %>% ungroup() %>% filter(Fuzz=="11") %>%
  group_by(Sample2, DNAInput,RefSNP) %>% 
  mutate(y3 = quantile(ReadCTotal, 0.03),y10 = quantile(ReadCTotal,0.1),
                       y15=quantile(ReadCTotal,0.15), y20 = quantile(ReadCTotal, 0.2)
                       ) -> allsen.3

# Locus depth for all loci
allsen.2 %>%  filter(Fuzz=="11") %>% filter(!(Pop2=="Control" | Pop2=="NTC" )) %>%
  mutate(rank = rank(-ReadCTotal, ties.method = "first")) %>%
  top_n(-2, rank) %>% rename(Population = Pop2) %>%
  mutate(DNAInput = factor(DNAInput, levels = c("31pg", "62pg", "125pg", "250pg", "500pg", "1ng"))) %>%
  ggplot(aes(x=DNAInput, y=log10(ReadCTotal), fill=Population)) +
  geom_boxplot() + scale_fill_viridis_d() + 
  labs(y = "Log based locus depth", x= "DNA input", title = "Locus depth across DNA inputs") +
  theme(text=element_text(size=15))#+
    facet_wrap(~Population)

# compute SD
allsen.2 %>%  filter(Fuzz=="11") %>% filter(!(Pop2=="Control" | Pop2=="NTC" )) %>%
  mutate(rank = rank(-ReadCTotal, ties.method = "first")) %>%
  top_n(-2, rank) %>% rename(Population = Pop2) %>% mutate(logReadCTotal = log10(ReadCTotal))-> foofoo #%>% group_by(Population, DNAInput) %>%
  
View(foofoo %>% filter(DNAInput == "1ng") %>% arrange(desc(ReadCTotal)))

foofoo %>% ungroup() %>% group_by(DNAInput) %>% summarise(SD = sd(log10(ReadCTotal),na.rm=TRUE)) -> SD
    
allsen.2 %>% ungroup() %>% filter(Fuzz=="11") %>% group_by(Sample2, DNAInput, RefSNP) %>%
  summarise(y3 = quantile(ReadCTotal, 0.03),y10 = quantile(ReadCTotal,0.1),
            y15=quantile(ReadCTotal,0.15), y20 = quantile(ReadCTotal, 0.2)) -> bar10

allsen.3 %>% filter(RefSNP == "rs907100" | 
                      RefSNP == "rs430046" | 
                      RefSNP == "rs1490413" | RefSNP == "rs1015250", Fuzz == "11") %>%
  #filter(DNAInput == "1ng" | DNAInput == "250pg" | DNAInput == "31pg") %>%
  filter(!(Pop2=="Control" | Pop2=="NTC" )) -> foobar1

foobar1 %>% mutate(rank = rank(-ReadCTotal, ties.method = "first")) %>%
  top_n(-2, rank) %>% mutate(row_id=row_number()) %>%
  mutate(logReadCTotal = log10(ReadCTotal)) -> foobar2

# locus depth for 4 loci (logreaddepth)
foobar2 %>%
ggplot(aes(x=RefSNP, y=logReadCTotal, fill = Pop2)) +
  #geom_bar(position = "dodge2", stat = "identity") +
  geom_boxplot() + scale_fill_viridis_d() +
  facet_wrap(~DNAInput)
  
foobar2 %>% mutate(y10log = log10(y10), y15log = log10(y15)) %>% 
  summarise(y10uniq = unique(y10log), y15uniq = unique(y15log)) %>%
  filter(!(Sample2 == "21832" |
      Sample2 == "22768" | 
      Sample2== "41131" |
      Sample2 == "5544")) -> foobar3 #%>%
  
foobar2 %>% filter(Fuzz=="11", 
                   !(Sample2 == "21832" |
                      Sample2 == "22768" | 
                      Sample2== "41131" |
                      Sample2 == "5544")) %>%
ggplot(aes(x=RefSNP, y = logReadCTotal)) +
  geom_bar(aes(fill=factor(rank)),stat = "identity", position = "dodge2", size = 1) + 
  facet_grid(Sample2~DNAInput, scales = "free") +
  #geom_segment(aes(x=row_id-0.5,xend=row_id+0.5,y=y10uniq,yend=y10uniq),
  #             inherit.aes = FALSE, color="orange")
  geom_errorbar(data= foobar3,aes(y=y10uniq,ymax=y10uniq,ymin=y10uniq), color="orange") +
  geom_errorbar(data= foobar3,aes(y=y15uniq,ymax=y15uniq,ymin=y15uniq), color="red") +
  geom_hline(yintercept = log10(50)) + scale_fill_manual(values = c("#73D055FF", "#F89441FF"))

allsen.2 %>% ungroup() %>% 
  group_by(Sample2, DNAInput,RefSNP) %>% 
  mutate(y3 = quantile(ReadCTotal, 0.03),y10 = quantile(ReadCTotal,0.1),
         y15=quantile(ReadCTotal,0.15), y20 = quantile(ReadCTotal, 0.2)
  ) 

amplicon_length <- c("42-49 bp","50-59 bp","61-69 bp","70-79 bp","80-89 bp","90-99 bp","100-128 bp")
Percentage <- c(10,24,14,14,8,15,15)
df <- data.frame(amplicon_length,Percentage)
ggplot(df, aes(area = Percentage,fill=amplicon_length, 
               label = paste(amplicon_length, paste(Percentage,"%"), sep = "\n"))) +
  geom_treemap() +
  geom_treemap_text(colour = c(rep("white", 5),1,"white"),
                    place = "centre", size = 19) +
  theme(legend.position = "none") + scale_fill_viridis_d(option="inferno")

SNP_amplicons <- data.frame(Length = c("101-128 bp", "42-100 bp"),
                            Percentage = c(13, 87))
ggplot(SNP_amplicons, aes(fill=Length, values = Percentage)) +
  geom_waffle(color="white") + geom_text(aes(x = 6, y = 5.5), label = "87%", size = 18) +
  geom_text(aes(x = 1.5, y = 2), label = "13%", size = 12) +
  coord_equal() + theme(legend.position = "bottom", legend.text=element_text(size=15)) +
  theme_void() + theme( 
                       legend.text=element_text(size=14), 
                       legend.title=element_text(size=14)) #+ theme_enhance_waffle()
